# 21_beta_pruebas.R
# Diseño definitivo, parte beta. Pruebas formales de la DIFERENCIACIÓN ambiental (presente e histórica) entre pares de celdas
# frente al recambio taxonómico y filogenético de Simpson, en los 4 géneros frugívoros.
#
# Unidad: par de celdas del mismo género. Respuestas: beta_tax y beta_filo (recambio de Simpson; el filogenético se promedia
# sobre n_arboles_beta árboles). Ambas están en [0, 1] y se modelan con un GLM cuasibinomial (enlace logit).
# Predictores de hipótesis (diferencia absoluta entre las dos celdas, estandarizada dentro de género):
#   presente  = |Δ topo_elev_media|      histórica = |Δ hist_sd_temp|
# Controles (pendiente común): |Δ clima_bio01_media|, |Δ clima_bio12_media|, |Δ ctrl_esfuerzo_log|, |Δ log S|, log S mínimo del par;
#   más un efecto de la distancia geográfica propio de cada género (spline natural de 3 gl).
#
# Pruebas (2 predictores x 2 respuestas x 2 = 8 por tamaño de celda):
#   (a) pendiente COMÚN: coeficiente del predictor en M0 (ambos predictores con pendiente común), z con errores robustos.
#   (b) INTERACCIÓN género x predictor: Wald de 3 gl de las desviaciones por género en M1 (común + desviación de ambos).
# Inferencia: los pares comparten celdas, así que no son independientes. Se usa la matriz de covarianza robusta "diádica"
#   (suma de las puntuaciones por celda: cada par aporta a sus dos celdas; las celdas se agrupan entre géneros por su identificador),
#   con corrección G/(G-1); p de pendiente con t de G-1 gl, p de interacción con chi-cuadrado de 3 gl.
#   Es un cambio respecto a lo anunciado (permutaciones): da las dos pruebas con un solo ajuste y evita miles de ajustes con ~300 mil pares.
#   Las permutaciones por bloques quedan como robustez.
#
# Celdas por género: las del criterio principal con al menos s_min especies (2) y presentes en el árbol; hasta max_celdas por género al azar.
# Salidas (datos/modelos): beta_pruebas_{tam}km.csv, beta_pendientes_genero_{tam}km.csv
# No probado: no hay R en mi entorno. Prueba con tamanos_km <- 200 (más rápido) y max_celdas <- 150, n_arboles_beta <- 3.
# Requiere: dplyr, readr, tidyr, ape, splines
# Para reutilizar las funciones: solo_funciones <- TRUE; source("21_beta_pruebas.R")

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr); library(readr); library(tidyr); library(ape); library(splines)

# terra/raster enmascaran select(); forzamos el de dplyr
select <- dplyr::select

dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_arb  <- file.path(dir_base, "datos", "filogenia")

archivo_arboles <- file.path(dir_arb, "arboles_completo_500.rds")
tamanos_km     <- c(100, 50, 200)
generos        <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
n_arboles_beta <- 20
max_celdas     <- 400
min_celdas     <- 30
s_min          <- 2
semilla        <- 1
bloque_km <- 500
# Predictor topografico = diferencia en la SD de elevacion, la misma medida de
# heterogeneidad que usan los modelos alfa. La diferencia en elevacion media
# pasa a control, porque dos celdas lejanas en el gradiente comparten menos
# especies por razones ajenas a la heterogeneidad.
hip   <- c(presente = "d_topo_sd_elev", historica = "d_hist_sd_temp")
ctrl  <- c("d_bio01", "d_bio12", "d_esf", "d_logS", "logS_min", "d_topo_elev_media")
vars_celda <- c("topo_elev_media", "topo_sd_elev", "hist_sd_temp", "clima_bio01_media", "clima_bio12_media", "ctrl_esfuerzo_log")

# ---- recambio de Simpson por pares (mismas funciones que 12 y 17) -----------------------------------------------
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
beta_taxonomico   <- function(X) beta_simpson(X %*% t(X), rowSums(X))
beta_filogenetico <- function(X, A, len) {
  P <- ((X %*% A) > 0) * 1
  a <- (P * rep(len, each = nrow(P))) %*% t(P)
  beta_simpson(a, as.vector(P %*% len))
}

# ---- pares de un género -------------------------------------------------------------------------------------------
construir_pares <- function(g, d, X, arboles) {
  n <- nrow(X); esp <- colnames(X)
  Mt <- beta_taxonomico(X)
  Mf <- Reduce(`+`, lapply(arboles, function(tr) {
    tr <- keep.tip(tr, esp); inc <- incidencia(tr)
    beta_filogenetico(X[, tr$tip.label, drop = FALSE], inc$A, inc$len)
  })) / length(arboles)
  ij <- which(lower.tri(Mt), arr.ind = TRUE); i <- ij[, 1]; j <- ij[, 2]
  dif <- function(x) abs(x[i] - x[j])
  lS <- log(d$S)
  tibble(genero = g, ci = as.character(d$celda[i]), cj = as.character(d$celda[j]),
         beta_tax = Mt[ij], beta_filo = Mf[ij],
         dist_km = sqrt((d$x_centro_km[i] - d$x_centro_km[j])^2 + (d$y_centro_km[i] - d$y_centro_km[j])^2),
         d_topo_elev_media = dif(d$topo_elev_media), d_topo_sd_elev = dif(d$topo_sd_elev),
         d_hist_sd_temp = dif(d$hist_sd_temp),
         d_bio01 = dif(d$clima_bio01_media), d_bio12 = dif(d$clima_bio12_media), d_esf = dif(d$ctrl_esfuerzo_log),
         d_logS = dif(lS), logS_min = pmin(lS[i], lS[j])) %>%
    filter(!is.na(beta_tax), !is.na(beta_filo))
}

# ---- GLM cuasibinomial con covarianza robusta diádica ---------------------------------------------------------------
ajuste_diadico <- function(X, y, ci, cj) {
  f <- suppressWarnings(glm.fit(X, y, family = quasibinomial()))
  b <- f$coefficients
  if (anyNA(b)) stop("columnas colineales en el diseño: ", paste(colnames(X)[is.na(b)], collapse = ", "))
  mu <- f$fitted.values; w <- mu * (1 - mu)
  bread <- solve(crossprod(X * sqrt(w)))
  U <- X * (y - mu)
  ids <- as.character(union(ci, cj))
  S <- matrix(0, length(ids), ncol(X), dimnames = list(ids, colnames(X)))
  for (idx in list(ci, cj)) { a <- rowsum(U, idx); S[rownames(a), ] <- S[rownames(a), ] + a }
  G <- length(ids)
  V <- bread %*% crossprod(S) %*% bread * G / (G - 1)
  dimnames(V) <- list(colnames(X), colnames(X))
  list(b = b, V = V, G = G)
}
pendiente <- function(a, col) {
  est <- a$b[[col]]; se <- sqrt(a$V[col, col]); t <- est / se
  list(est = est, se = se, ic = est + c(-1, 1) * qt(0.975, a$G - 1) * se, p = 2 * pt(-abs(t), a$G - 1))
}
wald_interaccion <- function(a, cols) {
  b <- a$b[cols]; V <- a$V[cols, cols, drop = FALSE]
  w <- as.numeric(t(b) %*% solve(V) %*% b)
  c(gl = length(cols), p = pchisq(w, length(cols), lower.tail = FALSE))
}

# ---- una respuesta ---------------------------------------------------------------------------------------------------
beta_respuesta <- function(P, y, etiqueta, tam_km, hipotesis = hip) {
  h <- unname(hipotesis)
  P$genero <- factor(P$genero)
  Xd <- function(f) model.matrix(f, P)
  nsd <- "genero * ns(dist_km, 3)"
  ctr <- paste(ctrl, collapse = " + ")
  f0 <- as.formula(paste("~", nsd, "+", ctr, "+", paste(h, collapse = " + ")))
  f1 <- as.formula(paste("~", nsd, "+", ctr, "+", paste(h, collapse = " + "), "+", paste0("genero:", h, collapse = " + ")))
  fs <- as.formula(paste("~", nsd, "+", ctr, "+", paste0("genero:", h, collapse = " + ")))
  a0 <- ajuste_diadico(Xd(f0), P[[y]], P$ci, P$cj)
  a1 <- ajuste_diadico(Xd(f1), P[[y]], P$ci, P$cj)
  as_ <- ajuste_diadico(Xd(fs), P[[y]], P$ci, P$cj)
  cat(sprintf("  %s: %d pares, %d celdas\n", etiqueta, nrow(P), a0$G))
  pruebas <- bind_rows(lapply(names(hipotesis), function(tipo) {
    p <- hipotesis[[tipo]]; s <- pendiente(a0, p)
    cols <- grep(paste0("^genero[^:]+:", p, "$"), colnames(a1$V), value = TRUE)
    wi <- wald_interaccion(a1, cols)
    tibble(tam_km = tam_km, respuesta = etiqueta, tipo = tipo, predictor = p, n_pares = nrow(P), n_celdas = a0$G,
           generos = paste(levels(P$genero), collapse = ","), est_comun = s$est, se_comun = s$se,
           ic_inf = s$ic[1], ic_sup = s$ic[2], p_comun = s$p, gl_int = wi[["gl"]], p_int = wi[["p"]])
  }))
  pend <- bind_rows(lapply(h, function(p) {
    cols <- grep(paste0("^genero[^:]+:", p, "$"), colnames(as_$V), value = TRUE)
    bind_rows(lapply(cols, function(cn) { s <- pendiente(as_, cn)
      tibble(tam_km = tam_km, respuesta = etiqueta, genero = sub("^genero([^:]+):.*$", "\\1", cn), predictor = p,
             estimate = s$est, se = s$se, p = s$p, ic_inf = s$ic[1], ic_sup = s$ic[2]) }))
  }))
  list(pruebas = pruebas, pendientes = pend)
}

# ---- un tamaño de celda ------------------------------------------------------------------------------------------------
beta_tamano <- function(base, pres, arboles, tam_km, generos, hipotesis = hip) {
  P <- pares_tamano(base, pres, arboles, tam_km, generos, hipotesis)
  pruebas_pares(P, tam_km, hipotesis)
}
pruebas_pares <- function(P, tam_km, hipotesis = hip) {
  r1 <- beta_respuesta(P, "beta_tax", "beta_tax", tam_km, hipotesis)
  r2 <- beta_respuesta(P, "beta_filo", "beta_filo", tam_km, hipotesis)
  list(pruebas = bind_rows(r1$pruebas, r2$pruebas), pendientes = bind_rows(r1$pendientes, r2$pendientes))
}
# Tabla de pares (estandarizada dentro de género) de un tamaño de celda
pares_tamano <- function(base, pres, arboles, tam_km, generos, hipotesis = hip) {
  tips_arbol <- arboles[[1]]$tip.label
  set.seed(semilla)
  lista <- list()
  for (g in generos) {
    d <- base %>% filter(genero == g) %>% filter(if_all(all_of(vars_celda), ~ !is.na(.x)))
    p <- pres %>% filter(genero_analisis == g, celda %in% d$celda, especie %in% tips_arbol)
    X <- p %>% distinct(celda, especie) %>% mutate(v = 1) %>% pivot_wider(names_from = especie, values_from = v, values_fill = 0)
    cel <- X$celda; X <- as.matrix(X[, -1]); rownames(X) <- cel
    X <- X[rowSums(X) >= s_min, , drop = FALSE]
    if (ncol(X) < 3 || nrow(X) < min_celdas) { cat(sprintf("  [%s] pocas especies o celdas (%d); se omite\n", g, nrow(X))); next }
    if (nrow(X) > max_celdas) X <- X[sort(sample.int(nrow(X), max_celdas)), , drop = FALSE]
    d <- d %>% filter(celda %in% rownames(X)) %>% arrange(match(celda, rownames(X)))
    stopifnot(identical(as.character(d$celda), rownames(X)))
    cat(sprintf("  [%s] celdas: %d | especies: %d | pares: %d\n", g, nrow(X), ncol(X), choose(nrow(X), 2)))
    pr <- construir_pares(g, d, X, arboles)
    est <- c(unname(hipotesis), ctrl)                     # estandarizar dentro de género (a nivel de par)
    pr[est] <- lapply(pr[est], function(x) as.numeric(scale(x)))
    lista[[g]] <- pr
  }
  bind_rows(lista)
}
# Agrupa las celdas de un tabla de pares en bloques espaciales (para la covarianza robusta por bloques)
a_bloques <- function(P, base, bloque_km = 500) {
  cc <- base %>% distinct(celda, x_centro_km, y_centro_km) %>%
    mutate(b = paste(floor(x_centro_km / bloque_km), floor(y_centro_km / bloque_km)))
  m <- setNames(cc$b, as.character(cc$celda))
  P$ci <- unname(m[P$ci]); P$cj <- unname(m[P$cj]); P
}

if (!exists("solo_funciones") || !solo_funciones) {
  arboles <- readRDS(archivo_arboles); arboles <- arboles[seq_len(min(n_arboles_beta, length(arboles)))]
  for (tam_km in tamanos_km) {
    cat("\n=== Celdas de", tam_km, "km ===\n")
    base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
    pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
    P <- pares_tamano(base, pres, arboles, tam_km, generos)
    # INFERENCIA PRIMARIA: errores agrupados por bloques de 500 km
    r <- pruebas_pares(a_bloques(P, base, bloque_km), tam_km)
    write_csv(r$pruebas,    file.path(dir_mod, sprintf("beta_pruebas_%dkm.csv", tam_km)))
    write_csv(r$pendientes, file.path(dir_mod, sprintf("beta_pendientes_genero_%dkm.csv", tam_km)))
    # SUPLEMENTARIO: errores agrupados por celda
    rc <- pruebas_pares(P, tam_km)
    write_csv(rc$pruebas,    file.path(dir_mod, sprintf("beta_pruebas_%dkm_celda.csv", tam_km)))
    write_csv(rc$pendientes, file.path(dir_mod, sprintf("beta_pendientes_genero_%dkm_celda.csv", tam_km)))
    cat("\nPruebas (est_comun = coeficiente logit por DE de la diferencia; p_comun con t de G-1 gl; p_int = Wald de 3 gl):\n")
    print(as.data.frame(r$pruebas %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
                          select(respuesta, tipo, predictor, n_pares, n_celdas, est_comun, ic_inf, ic_sup, p_comun, gl_int, p_int)), row.names = FALSE)
  }
}
