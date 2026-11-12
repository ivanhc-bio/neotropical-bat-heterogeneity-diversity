# 17_gdm_beta.R
# Fase 5 del plan: diversidad beta (recambio de Simpson) por pares de celdas con Modelos de Disimilitud Generalizados (GDM).
# Pregunta: ¿cuánto del recambio entre celdas de un género se explica por diferencias de clima actual, de topografía, de historia
# climática y por la distancia geográfica, y cuánto es aportado de forma única por cada grupo?
#
# Importante: esto NO es la misma pregunta que los scripts 12-13. Allí, el recambio medio de cada celda se modelaba contra la
# heterogeneidad DENTRO de la celda (hipótesis de heterogeneidad). Aquí el recambio entre dos celdas se modela contra las
# DIFERENCIAS entre ellas (gradientes). Los resultados son exploratorios (decisión de la fase 0: beta no es confirmatoria).
#
# Respuestas (por género, celdas del criterio principal, base_modelos):
#   tax    recambio taxonómico de Simpson entre pares de celdas
#   filo   recambio filogenético de Simpson (media de n_arboles_gdm árboles de la posterior)
# Grupos de predictores (diferencias entre celdas, transformadas por splines I de GDM):
#   clima_actual   clima_bio01_media, clima_bio12_media
#   topografia     topo_elev_media, topo_sd_elev
#   historia       hist_anom_temp, hist_sd_temp
#   esfuerzo       ctrl_esfuerzo_log (control: parte del recambio puede ser diferencia de muestreo)
#   geografia      distancia geográfica entre centros (geo = TRUE en GDM)
# Para cada género y respuesta: deviance explicada por el modelo completo y la contribución ÚNICA de cada grupo (completo menos
# el modelo sin ese grupo). La contribución única subestima a los grupos correlacionados (elevación y temperatura media, p. ej.).
# También se guarda la suma de los coeficientes de los splines por predictor (altura máxima de la curva de recambio) del completo.
#
# Cuidado:
#   - Se submuestrean hasta max_celdas celdas por género (n_pares crece con el cuadrado). Pon max_celdas <- Inf para usar todas.
#   - Los pares de GDM no son independientes; aquí no se prueban p-valores. Para eso: gdm.varImp (permutaciones, lento) con
#     usar_varImp <- TRUE (usa sampleSites = 0.5 y nPerm = n_perm).
#   - El recambio filogenético está correlacionado con el taxonómico; se imprime su correlación por género.
#
# Salidas (datos/modelos): gdm_resumen_{tam}km.csv, gdm_coeficientes_{tam}km.csv (y gdm_varimp_{género}_{resp}_{tam}km.csv si se pide)
# No probado: no tengo R en mi entorno. Prueba con tamanos_km <- 100, generos <- "Carollia", n_arboles_gdm <- 3, max_celdas <- 150.
# Requiere: gdm, ape, dplyr, readr, tidyr   (install.packages("gdm"))

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(gdm)
library(ape)
library(dplyr)
library(readr)
library(tidyr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_arb  <- file.path(dir_base, "datos", "filogenia")

etiqueta        <- "completo"
archivo_arboles <- file.path(dir_arb, sprintf("arboles_%s_500.rds", etiqueta))
tamanos_km      <- c(100)
generos         <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
n_arboles_gdm   <- 20
max_celdas      <- 400
min_celdas      <- 30
# Géneros con pocas especies (Carollia: 8, Anoura: 10) tienen casi todo el recambio en 0 o 1 y GDM puede no ajustar. Prueba una segunda
# corrida con celdas de al menos 2 especies: s_min_por_genero <- c(Carollia = 2, Anoura = 2) y sufijo_salida <- "_s2" (no pisa la primera).
s_min_por_genero <- c(Anoura = 2, Artibeus = 2, Carollia = 2, Platyrrhinus = 2, Sturnira = 2)
sufijo_salida   <- "_s2"
usar_varImp     <- FALSE
n_perm          <- 50
semilla         <- 1

grupos <- list(
  clima_actual = c("clima_bio01_media", "clima_bio12_media"),
  topografia   = c("topo_elev_media", "topo_sd_elev"),
  historia     = c("hist_anom_temp", "hist_sd_temp"),
  esfuerzo     = "ctrl_esfuerzo_log")
todas_vars <- unlist(grupos, use.names = FALSE)

# ---- 1. Funciones de recambio por pares (las mismas que en 12_diversidad_beta.R) -----------
incidencia <- function(tr) {
  n <- Ntip(tr); pp <- prop.part(tr); hijo <- tr$edge[, 2]
  A <- matrix(0, n, length(hijo))
  for (e in seq_along(hijo)) { h <- hijo[e]; A[if (h <= n) h else pp[[h - n]], e] <- 1 }
  list(A = A, len = tr$edge.length)
}
beta_simpson <- function(a, S) {
  b <- pmax(outer(S, rep(1, length(S))) - a, 0); c <- t(b)
  m <- pmin(b, c); den <- m + a
  out <- m / den; out[den <= 0] <- NA
  out
}
beta_taxonomico <- function(X) beta_simpson(X %*% t(X), rowSums(X))
beta_filogenetico <- function(X, A, len) {
  P <- ((X %*% A) > 0) * 1
  a <- (P * rep(len, each = nrow(P))) %*% t(P)
  beta_simpson(a, as.vector(P %*% len))
}

# ---- 2. GDM ------------------------------------------------------------------------------
# M: matriz de disimilitud con nombres de fila y columna = celda; pred: data.frame con celda, x, y y predictores
ajustar_gdm <- function(M, pred, vars, geo) {
  if (!length(vars) && !geo) return(NULL)
  # formatsitepair con matriz de disimilitud (bioFormat = 3) exige identificadores NUMÉRICOS: se usa el índice de la fila de M
  ids <- match(rownames(M), pred$celda)
  stopifnot(!anyNA(ids))
  M2 <- M; dimnames(M2) <- list(NULL, ids)
  bio <- data.frame(site = ids, M2, check.names = FALSE)
  # gdm::formatsitepair compara class(predData) con un texto; un tibble tiene varias clases y falla ("the condition has length > 1")
  pd <- as.data.frame(pred[ids, c("x", "y", vars), drop = FALSE])
  pd <- cbind(site = ids, pd)
  sp <- tryCatch(formatsitepair(bioData = bio, bioFormat = 3, siteColumn = "site", XColumn = "x", YColumn = "y", predData = pd),
                 error = function(e) { message("  formatsitepair: ", conditionMessage(e)); NULL })
  if (is.null(sp)) return(NULL)
  m <- tryCatch(gdm(sp, geo = geo), error = function(e) { message("  gdm: ", conditionMessage(e)); NULL })
  list(m = m, sp = sp)
}
explicada <- function(r) if (is.null(r) || is.null(r$m)) NA_real_ else r$m$explained

# Suma de coeficientes (3 splines I por predictor) = altura máxima de la curva de recambio
coef_por_pred <- function(m) {
  if (is.null(m)) return(NULL)
  # m$predictors tiene un nombre por predictor y m$coefficients un coeficiente por spline (m$splines por predictor)
  tibble(predictor = rep(m$predictors, times = m$splines), coef = m$coefficients) %>%
    group_by(predictor) %>% summarise(suma_coef = sum(coef), .groups = "drop") %>% arrange(desc(suma_coef))
}

# ---- 3. Datos comunes -------------------------------------------------------------------------
arboles <- readRDS(archivo_arboles)
arboles <- arboles[seq_len(min(n_arboles_gdm, length(arboles)))]
tips_arbol <- arboles[[1]]$tip.label
set.seed(semilla)

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  falta <- setdiff(todas_vars, names(base))
  if (length(falta)) stop("Faltan columnas en base_modelos_", tam_km, "km.csv: ", paste(falta, collapse = ", "),
                          ". Corre 05b, 05 y 09 de nuevo (versiones con hist_anom_temp y topo_elev_media).")
  res <- list(); coefs <- list()

  for (g in generos) {
    d <- base %>% filter(genero == g) %>% filter(if_all(all_of(todas_vars), ~ !is.na(.x)))
    p <- pres %>% filter(genero_analisis == g, celda %in% d$celda, especie %in% tips_arbol)
    X <- p %>% distinct(celda, especie) %>% mutate(v = 1) %>% pivot_wider(names_from = especie, values_from = v, values_fill = 0)
    cel <- X$celda; X <- as.matrix(X[, -1]); rownames(X) <- cel
    s_min <- if (g %in% names(s_min_por_genero)) s_min_por_genero[[g]] else 1
    X <- X[rowSums(X) >= s_min, , drop = FALSE]
    if (ncol(X) < 3 || nrow(X) < min_celdas) { cat(sprintf("  [%s] pocas especies o celdas (%d); se omite\n", g, nrow(X))); next }
    if (nrow(X) > max_celdas) X <- X[sort(sample.int(nrow(X), max_celdas)), , drop = FALSE]
    esp <- colnames(X)
    d <- d %>% filter(celda %in% rownames(X)) %>% arrange(match(celda, rownames(X)))
    stopifnot(identical(d$celda, rownames(X)))
    cat(sprintf("  [%s] celdas: %d | especies: %d | pares: %d\n", g, nrow(X), length(esp), choose(nrow(X), 2)))

    pred <- d %>% transmute(celda, x = x_centro_km, y = y_centro_km, across(all_of(todas_vars)))

    Mt <- beta_taxonomico(X); dimnames(Mt) <- list(rownames(X), rownames(X))
    Mf <- Reduce(`+`, lapply(arboles, function(tr) {
      tr <- keep.tip(tr, esp); inc <- incidencia(tr)
      beta_filogenetico(X[, tr$tip.label, drop = FALSE], inc$A, inc$len)
    })) / length(arboles)
    dimnames(Mf) <- list(rownames(X), rownames(X))
    diag(Mt) <- 0; diag(Mf) <- 0
    Mt[is.na(Mt)] <- 0; Mf[is.na(Mf)] <- 0     # pares sin información (ambas celdas vacías): no ocurre con S >= 1
    cat(sprintf("    correlación recambio taxonómico-filogenético (pares): %.2f\n", cor(Mt[lower.tri(Mt)], Mf[lower.tri(Mf)])))

    for (rn in c("tax", "filo")) {
      M <- if (rn == "tax") Mt else Mf
      full <- ajustar_gdm(M, pred, todas_vars, geo = TRUE)
      if (is.null(full) || is.null(full$m)) { cat(sprintf("    [%s] el modelo completo no se ajustó\n", rn)); next }
      fila <- tibble(tam_km = tam_km, genero = g, respuesta = rn, n_celdas = nrow(X), n_pares = nrow(full$sp),
                     explicada_completo = explicada(full))
      # contribución única: completo menos el modelo sin el grupo
      for (gr in names(grupos)) {
        sin <- ajustar_gdm(M, pred, setdiff(todas_vars, grupos[[gr]]), geo = TRUE)
        fila[[paste0("unica_", gr)]] <- explicada(full) - explicada(sin)
      }
      sin_geo <- ajustar_gdm(M, pred, todas_vars, geo = FALSE)
      fila$unica_geografia <- explicada(full) - explicada(sin_geo)
      res[[paste(g, rn)]] <- fila
      coefs[[paste(g, rn)]] <- coef_por_pred(full$m) %>% mutate(tam_km = tam_km, genero = g, respuesta = rn)
      cat(sprintf("    [%s] deviance explicada (completo) = %.1f %%\n", rn, explicada(full)))

      if (usar_varImp) {
        vi <- tryCatch(gdm.varImp(full$sp, geo = TRUE, nPerm = n_perm, sampleSites = 0.5, parallel = FALSE),
                       error = function(e) { message("  gdm.varImp: ", conditionMessage(e)); NULL })
        if (!is.null(vi)) saveRDS(vi, file.path(dir_mod, sprintf("gdm_varimp_%s_%s_%dkm%s.rds", g, rn, tam_km, sufijo_salida)))
      }
    }
  }

  res <- bind_rows(res); coefs <- bind_rows(coefs)
  if (!nrow(res)) { cat("\nNingún modelo se ajustó; revisa los mensajes de arriba.\n"); next }
  write_csv(res, file.path(dir_mod, sprintf("gdm_resumen_%dkm%s.csv", tam_km, sufijo_salida)))
  write_csv(coefs, file.path(dir_mod, sprintf("gdm_coeficientes_%dkm%s.csv", tam_km, sufijo_salida)))

  # ---- 4. Resumen legible ----
  r1 <- function(x) round(x, 1)
  cat("\n--- Deviance explicada (%) por el modelo completo y contribución única de cada grupo (puntos porcentuales) ---\n")
  print(as.data.frame(res %>% mutate(across(where(is.numeric) & !c(tam_km, n_celdas, n_pares), r1)) %>%
                        select(genero, respuesta, n_celdas, n_pares, explicada_completo, starts_with("unica_"))), row.names = FALSE)
  cat("\n--- Altura máxima de la curva de recambio por predictor (suma de coeficientes; mayor = más recambio con esa diferencia) ---\n")
  print(as.data.frame(coefs %>% mutate(suma_coef = round(suma_coef, 3)) %>% select(genero, respuesta, predictor, suma_coef) %>%
                        pivot_wider(names_from = predictor, values_from = suma_coef)), row.names = FALSE)
  cat("\nGuardado:", file.path(dir_mod, sprintf("gdm_resumen_%dkm%s.csv", tam_km, sufijo_salida)), "\n")
}
