# 11_modelo_conjunto.R
# Modelo conjunto con interacción por género (GAM), a partir de base_modelos_{tam}km.csv (09) y presencia_{tam}km.csv.
# Responde: ¿el efecto de cada predictor es común a los géneros o específico de alguno?
#
#   1. Pruebas de interacción género x predictor. Con controles_comunes = TRUE (por defecto, igual que 15 y 18):
#        los CONTROLES (esfuerzo, área, clima medio, EVI, log_S) tienen una pendiente común a los géneros y solo los
#        predictores de hipótesis pueden variar por género.
#        Prueba global: nulo = todos los predictores de hipótesis comunes; completo = todos por género.
#        Prueba por predictor: completo = solo ese predictor por género (los demás, comunes); nulo = todos comunes.
#        Se reportan la razón de verosimilitud (ML; no válida si la devianza <= 0), el dAIC y la prueba de Wald (REML).
#        La prueba de Wald (p_wald, gl_wald = n_géneros - 1) contrasta la interacción: se ajusta un modelo REML con una
#        pendiente común del predictor MÁS su desviación por género (genero:predictor) y se prueba esa desviación.
#        (Versiones anteriores de este script probaban genero:predictor sin pendiente común = 'todas las pendientes
#        son cero', con n_géneros gl; esa prueba no mide si las pendientes difieren entre géneros y fue corregida.)
#      Con controles_comunes = FALSE se reproduce la versión anterior, donde todo predictor no común tenía pendiente por
#      género (unas 50 pendientes por género con 5 géneros y 10 predictores): ahí H2 se diluye (ver script 18).
#      Como los predictores se estandarizan dentro de cada género, una pendiente es "por DE dentro del género".
#      Esto evita la dependencia del género de referencia que tiene la prueba con factor ordenado.
#   2. Respuestas: riqueza (binomial, cbind(S, S_fallo)) y SES-MPD (gaussiano, modelo reducido, log_S como covariable).
#   3. Partición de devianza por género (A = control, B = A + espacio, C = A + ambiente, F = todo):
#        ambiente único = F - B;  espacio único = F - C;  compartido = B + C - F - A;  residual = 1 - F.
#      Con concurvidad alta del suavizado espacial (~0.9) el reparto entre "espacio único" y "compartido" es poco estable:
#      léelo como orden de magnitud, no como cifra exacta.
#   4. Sensibilidad no lineal: hist_sd_temp como suavizado s(x, k = 5) por género vs lineal (AIC y p del suavizado).
#
# Salidas (datos/modelos): prueba_global_conjunto_{tam}km{sufijo}.csv, coef_conjunto_{tam}km{sufijo}.csv
#                          (sufijo = "_cc" con controles comunes, "" con la versión anterior, para no sobrescribir),
#                          y, si correr_particion / correr_no_lineal = TRUE, particion_devianza_{tam}km.csv y
#                          no_lineal_hist_{tam}km.csv (no dependen de controles_comunes).
#
# No probado: no tengo R en mi entorno. Corre con tamanos_km <- 100 y revisa lo que imprime.
# Requiere: mgcv, ape, dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv)
library(ape)
library(dplyr)
library(readr)
library(tidyr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")

todas_las_escalas <- TRUE   # TRUE = corre 100, 50 y 200 km (en ese orden); FALSE = solo 100 km
tamanos_km <- if (todas_las_escalas) c(100, 50, 200) else c(100)
generos      <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
esfuerzo_var <- "ctrl_esfuerzo_log"
controles_comunes <- FALSE   # TRUE = controles con pendiente común en las pruebas de interacción (como 15 y 18)
correr_particion  <- FALSE  # la partición de devianza no cambia con controles_comunes; TRUE solo para repetirla
correr_no_lineal  <- FALSE  # ídem para la sensibilidad no lineal de hist_sd_temp
sufijo_salida     <- if (controles_comunes) "_cc" else ""
n_min_modelo <- 25          # mínimo de celdas por género para entrar al modelo conjunto
k_por_celdas <- 8           # k del suavizado espacial por género = n / 8 (entre 5 y 40)

pred_completo  <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec",
                    "clima_bio01_media", "clima_bio12_media", "hab_evi_media", "ctrl_esfuerzo_log", "ctrl_area_tierra_km2")
pred_reducido  <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp", "clima_bio01_media", "ctrl_esfuerzo_log")
hip_riqueza    <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec")
hip_ses        <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp")

# ---- 1. Funciones -------------------------------------------------------------------
zscore <- function(x) as.numeric(scale(x))

preparar <- function(d, vars) {
  d %>% group_by(genero) %>% mutate(across(all_of(vars), zscore)) %>% ungroup() %>%
    mutate(genero = factor(genero)) %>%
    bind_cols(setNames(lapply(levels(factor(d$genero)), function(g) as.numeric(d$genero == g)),
                       paste0("I_", levels(factor(d$genero)))))
}

# k espacial propio por género (n / 8, entre 5 y 40); se aplica con una variable indicadora I_<género>
k_esp_por_genero <- function(d) {
  n <- table(d$genero)
  setNames(pmax(5, pmin(40, floor(as.numeric(n) / k_por_celdas))), names(n))
}

# fórmula: predictores "comunes" con una sola pendiente, el resto con pendiente por género
# sep_extra: predictores que, además de su pendiente común (si está en 'comunes'), llevan la desviación por género
formula_conjunta <- function(y, vars, comunes, k_esp, extra_suave = NULL, sep_extra = NULL) {
  sep <- union(setdiff(vars, comunes), sep_extra)
  partes <- c("genero",
              if (length(comunes)) comunes,
              if (length(sep)) paste0("genero:", sep),
              sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k_esp), as.integer(k_esp)),
              extra_suave)
  as.formula(paste0(y, " ~ ", paste(partes, collapse = " + ")))
}

ajustar_gam <- function(f, d, fam, metodo) {
  tryCatch(gam(f, family = fam, data = d, method = metodo),
           error = function(e) gam(f, family = fam, data = d, method = "GCV.Cp"))
}

moran_por_genero <- function(m, d) {
  res <- residuals(m, type = "pearson")
  sapply(levels(d$genero), function(g) {
    i <- d$genero == g
    w <- 1 / as.matrix(dist(d[i, c("x_centro_km", "y_centro_km")])); diag(w) <- 0
    mi <- Moran.I(res[i], w); mi$p.value
  })
}

# prueba de interacción: global (todos los predictores de hipótesis comunes frente a por género) y uno por uno
pruebas_interaccion <- function(d, y, vars, hip, fam, k_esp, etiqueta, ctrl_comunes = TRUE) {
  ctrl <- setdiff(vars, hip)
  # predictores con pendiente común en cada modelo
  com_full_global <- if (ctrl_comunes) ctrl else character()
  com_nulo_global <- if (ctrl_comunes) c(ctrl, hip) else hip
  com_full_p <- function(p) if (ctrl_comunes) c(ctrl, setdiff(hip, p)) else character()
  com_nulo_p <- function(p) if (ctrl_comunes) c(ctrl, hip) else p
  f <- function(com, metodo) ajustar_gam(formula_conjunta(y, vars, com, k_esp), d, fam, metodo)
  comp <- function(full, nulo) {
    a <- tryCatch(anova(nulo, full, test = "Chisq"), error = function(e) NULL)
    dev <- if (is.null(a)) NA_real_ else a$Deviance[2]
    p <- if (is.null(a)) NA_real_ else a[["Pr(>Chi)"]][2]
    ok <- !is.na(dev) && dev > 0                   # devianza <= 0 = el ajuste ML falló: se usa solo el AIC
    tibble(gl = if (is.null(a)) NA_real_ else a$Df[2], devianza = dev, p = if (ok) p else NA_real_,
           dAIC = AIC(full) - AIC(nulo),           # <0 = el modelo con interacción es mejor
           nota = if (ok) "" else "prueba LR no válida; usar p_wald y dAIC")
  }
  wald <- function(p) {                            # Wald de la interacción genero:p (pendiente común + desviación por género)
    com <- c(com_full_p(p), p)                     # p con pendiente común; el resto como en el modelo 'completo' del predictor
    m <- ajustar_gam(formula_conjunta(y, vars, com, k_esp, sep_extra = p), d, fam, "REML")
    at <- tryCatch(anova(m)$pTerms.table, error = function(e) NULL)
    fila <- if (is.null(at)) integer() else which(rownames(at) == paste0("genero:", p))
    if (length(fila) == 1) tibble(gl_wald = at[fila, "df"], p_wald = at[fila, "p-value"])
    else tibble(gl_wald = NA_real_, p_wald = NA_real_)
  }
  nulo_g <- f(com_nulo_global, "ML")
  out <- list(bind_cols(tibble(prueba = "global (todos los predictores de hipótesis)"),
                        comp(f(com_full_global, "ML"), nulo_g), tibble(gl_wald = NA_real_, p_wald = NA_real_)))
  for (p in hip) {
    nulo_p <- if (ctrl_comunes) nulo_g else f(com_nulo_p(p), "ML")
    full_p <- f(com_full_p(p), "ML")
    out[[p]] <- bind_cols(tibble(prueba = p), comp(full_p, nulo_p), wald(p))
  }
  bind_rows(out) %>% mutate(respuesta = etiqueta, n = nrow(d), controles_comunes = ctrl_comunes,
                            generos = paste(levels(d$genero), collapse = ","))
}

# coeficientes por género del modelo completo (REML)
coef_conjunto <- function(d, y, vars, fam, k_esp, etiqueta, hip, ctrl_comunes = TRUE) {
  comunes <- if (ctrl_comunes) setdiff(vars, hip) else character()
  m <- ajustar_gam(formula_conjunta(y, vars, comunes, k_esp), d, fam, "REML")
  tab <- as.data.frame(summary(m)$p.table); names(tab) <- c("estimate", "se", "stat", "p")
  tab$term <- rownames(tab)
  tab <- tab %>% filter(grepl("^genero.*:", term)) %>%
    mutate(genero = sub("^genero([^:]+):.*$", "\\1", term), predictor = sub("^genero[^:]+:", "", term),
           ic_inf = estimate - 1.96 * se, ic_sup = estimate + 1.96 * se, respuesta = etiqueta) %>%
    select(respuesta, genero, predictor, estimate, se, stat, p, ic_inf, ic_sup)
  mp <- moran_por_genero(m, d)
  cat(sprintf("  [%s] dev.expl = %.2f; Moran p por género: %s\n", etiqueta, summary(m)$dev.expl,
              paste(names(mp), sprintf("%.3f", mp), collapse = ", ")))
  tab
}

# partición de devianza (dev. explicada) para un género
particion <- function(d, y, ctrl, env, fam, k_esp) {
  ajusta <- function(terms, espacio) {
    t <- c(terms, if (espacio) sprintf("s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", k_esp))
    if (!length(t)) t <- "1"
    m <- ajustar_gam(as.formula(paste0(y, " ~ ", paste(t, collapse = " + "))), d, fam, "REML")
    summary(m)$dev.expl
  }
  A <- ajusta(ctrl, FALSE); B <- ajusta(ctrl, TRUE); C <- ajusta(c(ctrl, env), FALSE); Fm <- ajusta(c(ctrl, env), TRUE)
  tibble(A_control = A, B_control_espacio = B, C_control_ambiente = C, F_todo = Fm,
         ambiente_unico = Fm - B, espacio_unico = Fm - C, compartido = B + C - Fm - A,
         control = A, residual = 1 - Fm)
}

# ---- 2. Por tamaño de celda ------------------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)

  # grupo de especies (pool) por género y S_fallo
  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g)
    if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    if (any(d$S_fallo < 0)) stop("Hay celdas con más especies que el grupo de especies (", g, ")")
    d$pool <- pool
    d
  }))

  pruebas <- list(); coefs <- list(); parts <- list(); nolin <- list()

  # ---- 2a. Riqueza: modelo conjunto -----------------------------------------------------
  vars_r <- sub("^ctrl_esfuerzo_log$", esfuerzo_var, pred_completo)
  hip_r  <- hip_riqueza
  gr <- names(which(table(base$genero) >= n_min_modelo))
  dr <- base %>% filter(genero %in% gr) %>% preparar(vars_r)
  k_r <- k_esp_por_genero(dr)
  cat(sprintf("\nRiqueza: géneros = %s; n = %d; k espacial por género: %s\n", paste(gr, collapse = ", "), nrow(dr), paste(names(k_r), k_r, collapse = ", ")))
  pruebas[["riqueza"]] <- pruebas_interaccion(dr, "cbind(S, S_fallo)", vars_r, hip_r, binomial(), k_r, "riqueza", controles_comunes)
  coefs[["riqueza"]]   <- coef_conjunto(dr, "cbind(S, S_fallo)", vars_r, binomial(), k_r, "riqueza", hip_r, controles_comunes)

  # ---- 2b. SES-MPD: modelo conjunto (predictores reducidos + log_S) ------------------------
  vars_s <- c(sub("^ctrl_esfuerzo_log$", esfuerzo_var, pred_reducido), "log_S")
  ds <- base %>% filter(filo_valida, !is.na(ses_mpd))
  gs <- names(which(table(ds$genero) >= n_min_modelo))
  ds <- ds %>% filter(genero %in% gs) %>% preparar(vars_s)
  k_s <- k_esp_por_genero(ds)
  cat(sprintf("\nSES-MPD: géneros = %s; n = %d; k espacial por género: %s\n", paste(gs, collapse = ", "), nrow(ds), paste(names(k_s), k_s, collapse = ", ")))
  if (length(gs) >= 2) {
    pruebas[["ses"]] <- pruebas_interaccion(ds, "ses_mpd", vars_s, hip_ses, gaussian(), k_s, "ses_mpd", controles_comunes)
    coefs[["ses"]]   <- coef_conjunto(ds, "ses_mpd", vars_s, gaussian(), k_s, "ses_mpd", hip_ses, controles_comunes)
  }

  # ---- 2c. Partición de devianza por género --------------------------------------------------
  ctrl_r <- c(esfuerzo_var, "ctrl_area_tierra_km2")
  env_r  <- setdiff(vars_r, ctrl_r)
  ctrl_s <- c(esfuerzo_var, "log_S")
  env_s  <- setdiff(vars_s, ctrl_s)
  if (correr_particion) for (g in gr) {
    dg <- dr %>% filter(genero == g)
    k <- max(5, min(60, floor(nrow(dg) / k_por_celdas)))
    parts[[paste(g, "riqueza")]] <- particion(dg, "cbind(S, S_fallo)", ctrl_r, env_r, binomial(), k) %>%
      mutate(genero = g, respuesta = "riqueza", n = nrow(dg))
  }
  if (correr_particion) for (g in gs) {
    dg <- ds %>% filter(genero == g)
    k <- max(5, min(60, floor(nrow(dg) / k_por_celdas)))
    parts[[paste(g, "ses")]] <- particion(dg, "ses_mpd", ctrl_s, env_s, gaussian(), k) %>%
      mutate(genero = g, respuesta = "ses_mpd", n = nrow(dg))
  }

  # ---- 2d. Sensibilidad no lineal de hist_sd_temp (por género, riqueza) -------------------------
  if (correr_no_lineal) for (g in gr) {
    dg <- dr %>% filter(genero == g)
    k <- max(5, min(60, floor(nrow(dg) / k_por_celdas)))
    base_t <- c(setdiff(vars_r, "hist_sd_temp"), sprintf("s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", k))
    f_lin <- as.formula(paste("cbind(S, S_fallo) ~", paste(c(base_t, "hist_sd_temp"), collapse = " + ")))
    f_sua <- as.formula(paste("cbind(S, S_fallo) ~", paste(c(base_t, "s(hist_sd_temp, k = 5)"), collapse = " + ")))
    m_l <- ajustar_gam(f_lin, dg, binomial(), "REML"); m_s <- ajustar_gam(f_sua, dg, binomial(), "REML")
    st <- summary(m_s)$s.table
    fila <- st[grep("hist_sd_temp", rownames(st)), , drop = FALSE]
    nolin[[g]] <- tibble(genero = g, respuesta = "riqueza", edf = fila[1, "edf"], p_suavizado = fila[1, "p-value"],
                         AIC_lineal = AIC(m_l), AIC_suavizado = AIC(m_s), dAIC = AIC(m_s) - AIC(m_l))
  }

  # ---- 2e. Guardar e imprimir ----------------------------------------------------------------------
  pr <- bind_rows(pruebas) %>% mutate(tam_km = tam_km) %>% select(tam_km, respuesta, prueba, everything())
  co <- bind_rows(coefs) %>% mutate(tam_km = tam_km) %>% select(tam_km, everything())
  write_csv(pr, file.path(dir_mod, sprintf("prueba_global_conjunto_%dkm%s.csv", tam_km, sufijo_salida)))
  write_csv(co, file.path(dir_mod, sprintf("coef_conjunto_%dkm%s.csv", tam_km, sufijo_salida)))
  redondea <- function(x) x %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
  cat("\nControles comunes:", controles_comunes, "\n")
  cat("Pruebas de interacción (p_wald = Wald de la interacción genero x predictor, gl = n_géneros - 1; p_wald pequeña = el efecto difiere entre géneros; dAIC < 0 favorece la interacción):\n")
  print(as.data.frame(redondea(pr %>% select(respuesta, prueba, gl, devianza, p, gl_wald, p_wald, dAIC, nota))), row.names = FALSE)
  cat("\nPendientes por género (por DE dentro del género), predictores de hipótesis:\n")
  print(as.data.frame(redondea(co %>% filter(predictor %in% c(hip_riqueza, hip_ses)) %>%
                                 mutate(estimado = sprintf("%.2f (%.2f)", estimate, se)) %>%
                                 select(respuesta, genero, predictor, estimado, p))), row.names = FALSE)
  if (correr_particion) {
    pa <- bind_rows(parts) %>% mutate(tam_km = tam_km) %>% select(tam_km, respuesta, genero, n, everything())
    write_csv(pa, file.path(dir_mod, sprintf("particion_devianza_%dkm.csv", tam_km)))
    cat("\nPartición de devianza (fracción de la devianza total):\n")
    print(as.data.frame(redondea(pa %>% select(respuesta, genero, n, control, ambiente_unico, espacio_unico, compartido, residual))), row.names = FALSE)
  }
  if (correr_no_lineal) {
    nl <- bind_rows(nolin) %>% mutate(tam_km = tam_km) %>% select(tam_km, everything())
    write_csv(nl, file.path(dir_mod, sprintf("no_lineal_hist_%dkm.csv", tam_km)))
    cat("\nhist_sd_temp: ¿no lineal? (edf ~ 1 = lineal; dAIC < 0 favorece el suavizado)\n")
    print(as.data.frame(redondea(nl)), row.names = FALSE)
  }
}
