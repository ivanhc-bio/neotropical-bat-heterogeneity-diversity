# 18_h2_especificacion.R
# Diagnóstico: ¿por qué H2 (el efecto de topo_sd_elev sobre la riqueza difiere entre géneros) sale clara en el script 15 y no en el 11?
#   15: 4 géneros frugívoros, 7 predictores (3 de hipótesis + 4 controles)                 -> Wald de la interacción p ~ .001 (100 km)
#   11: 5 géneros (con Anoura), 10 predictores (5 de hipótesis + 5 controles)              -> prueba de interacción p ~ .5-.8 (50 y 200 km)
# Diferencias entre ambos: (a) Anoura incluida o no, (b) conjunto de predictores (7 o 10), (c) la prueba (Wald en REML frente a razón de
# verosimilitud en ML). Aquí se cruzan (a) y (b) y se calculan las dos pruebas en cada celda del cruce.
#
# Especificaciones (todas con el suavizado espacial "actual" de 15: s(x, y) por género con k = n / 8):
#   A_7pred_4gen     como 15
#   B_7pred_5gen     como 15 más Anoura
#   C_10pred_4gen    predictores del 11, sin Anoura
#   D_10pred_5gen    como 11
# Respuesta: riqueza (binomial). H1 = pendiente común de topo_sd_elev; H2 = interacción género x topo_sd_elev.
#
# Salida (datos/modelos): h2_especificacion_{tam}km.csv
# No probado: no tengo R en mi entorno. Debe reproducir H1 y H2 de 15 en A_7pred_4gen (100 km: H1 = 0.115, H2 Wald p = .001).
# Requiere: mgcv, dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv)
library(dplyr)
library(readr)
library(tidyr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")

tamanos_km   <- c(50, 100, 200)
frugivoros   <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
con_anoura   <- c("Anoura", frugivoros)
n_min_modelo <- 25
k_por_celdas <- 8
h_conf       <- "topo_sd_elev"

hip_7   <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp")
ctrl_7  <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2", "clima_bio01_media", "clima_bio12_media")
hip_10  <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec")
ctrl_10 <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2", "clima_bio01_media", "clima_bio12_media", "hab_evi_media")
especificaciones <- list(
  A_7pred_4gen  = list(hip = hip_7,  ctrl = ctrl_7,  gens = frugivoros),
  B_7pred_5gen  = list(hip = hip_7,  ctrl = ctrl_7,  gens = con_anoura),
  C_10pred_4gen = list(hip = hip_10, ctrl = ctrl_10, gens = frugivoros),
  D_10pred_5gen = list(hip = hip_10, ctrl = ctrl_10, gens = con_anoura))

# ---- 1. Funciones (las de 15 y 16) ----------------------------------------------------
zscore <- function(x) as.numeric(scale(x))
preparar <- function(d, vars) {
  gens <- sort(unique(d$genero))
  d <- d %>% group_by(genero) %>% mutate(across(all_of(vars), zscore)) %>% ungroup()
  for (g in gens) d[[paste0("I_", g)]] <- as.numeric(d$genero == g)
  d$genero <- factor(d$genero)
  d
}
k_por_genero <- function(d) { n <- table(d$genero); setNames(pmax(5, pmin(40, floor(as.numeric(n) / k_por_celdas))), names(n)) }
formula_conjunta <- function(y, d, ctrl, separados, comunes_hip) {
  k <- k_por_genero(d)
  esp <- sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k), as.integer(k))
  partes <- c("genero", ctrl, comunes_hip, if (length(separados)) paste0("genero:", separados), esp)
  as.formula(paste(y, "~", paste(partes, collapse = " + ")))
}
ajusta <- function(f, d, metodo) tryCatch(gam(f, family = binomial(), data = d, method = metodo), error = function(e) NULL)

ajustar_especificacion <- function(base, id, sp) {
  vars <- c(sp$hip, sp$ctrl)
  d0 <- base %>% filter(genero %in% sp$gens) %>% filter(if_all(all_of(vars), ~ !is.na(.x)))
  ok <- names(which(table(d0$genero) >= n_min_modelo))
  d <- preparar(d0 %>% filter(genero %in% ok), vars)
  y <- "cbind(S, S_fallo)"
  m_com <- ajusta(formula_conjunta(y, d, sp$ctrl, character(), sp$hip), d, "REML")
  m_gen <- ajusta(formula_conjunta(y, d, sp$ctrl, h_conf, sp$hip), d, "REML")             # común + interacción (Wald)
  m_sep <- ajusta(formula_conjunta(y, d, sp$ctrl, h_conf, setdiff(sp$hip, h_conf)), d, "REML")  # pendiente por género (para reportar)
  mf <- ajusta(formula_conjunta(y, d, sp$ctrl, h_conf, setdiff(sp$hip, h_conf)), d, "ML")
  mn <- ajusta(formula_conjunta(y, d, sp$ctrl, character(), sp$hip), d, "ML")
  if (is.null(m_com)) return(NULL)

  tc <- as.data.frame(summary(m_com)$p.table)[h_conf, ]
  p_wald <- NA_real_
  if (!is.null(m_gen)) {
    at <- tryCatch(anova(m_gen)$pTerms.table, error = function(e) NULL)
    fila <- if (is.null(at)) integer() else which(rownames(at) == paste0("genero:", h_conf))
    if (length(fila) == 1) p_wald <- at[fila, "p-value"]
  }
  p_lr <- NA_real_; d_aic <- NA_real_
  if (!is.null(mf) && !is.null(mn)) {
    a <- tryCatch(anova(mn, mf, test = "Chisq"), error = function(e) NULL)
    d_aic <- AIC(mf) - AIC(mn)
    if (!is.null(a) && !is.na(a$Deviance[2]) && a$Deviance[2] > 0) p_lr <- a[["Pr(>Chi)"]][2]
  }
  por_gen <- NA_character_
  if (!is.null(m_sep)) {
    tg <- as.data.frame(summary(m_sep)$p.table); tg$nombre <- rownames(tg)
    tg <- tg[grepl(paste0("^genero.*:", h_conf, "$"), tg$nombre), ]
    por_gen <- paste(sprintf("%s %.2f (p %.3f)", sub("^genero([^:]+):.*$", "\\1", tg$nombre), tg$Estimate, tg[[4]]), collapse = "; ")
  }
  tibble(especificacion = id, n_predictores = length(vars), generos = paste(ok, collapse = ","), n = nrow(d),
         H1_estimado = tc[[1]], H1_se = tc[[2]], H1_p = tc[[4]],
         H2_p_wald = p_wald, H2_p_LR = p_lr, H2_dAIC = d_aic, por_genero = por_gen)
}

# ---- 2. Por tamaño de celda --------------------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  base <- bind_rows(lapply(con_anoura, function(g) {
    d <- base %>% filter(genero == g)
    if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    d
  }))
  falta <- setdiff(unique(c(unlist(lapply(especificaciones, `[[`, "hip")), unlist(lapply(especificaciones, `[[`, "ctrl")))), names(base))
  if (length(falta)) stop("Faltan columnas en base_modelos_", tam_km, "km.csv: ", paste(falta, collapse = ", "))

  res <- bind_rows(lapply(names(especificaciones), function(id) ajustar_especificacion(base, id, especificaciones[[id]]))) %>%
    mutate(tam_km = tam_km)
  write_csv(res, file.path(dir_mod, sprintf("h2_especificacion_%dkm.csv", tam_km)))

  r3 <- function(x) round(x, 3)
  cat("\n--- H1 (pendiente común de topo_sd_elev) y H2 (interacción género x topo_sd_elev) según la especificación ---\n")
  print(as.data.frame(res %>% mutate(H1 = sprintf("%.3f (%.3f)", H1_estimado, H1_se), H1_p = r3(H1_p), H2_p_wald = r3(H2_p_wald),
                                     H2_p_LR = r3(H2_p_LR), H2_dAIC = round(H2_dAIC, 1)) %>%
                        select(especificacion, n_predictores, n, H1, H1_p, H2_p_wald, H2_p_LR, H2_dAIC)), row.names = FALSE)
  cat("\n--- Pendientes por género de topo_sd_elev ---\n")
  print(as.data.frame(res %>% select(especificacion, por_genero)), row.names = FALSE)
  cat("\nGuardado:", file.path(dir_mod, sprintf("h2_especificacion_%dkm.csv", tam_km)), "\n")
}
