# 18b_h2_wald_controles_por_genero.R
# Prueba de interacción genero x topo_sd_elev (Wald, 3 gl) cuando TODOS los controles y los demás predictores de hipótesis
# tienen pendiente propia por género (la versión sin _cc del script 11). Riqueza, 4 frugívoros, 10 predictores.
# Corrige el p_wald del script 11, que prueba "las 4 pendientes = 0" (4 gl), no la igualdad entre géneros.
# Mismo modelo que el script 11 sin _cc; la diferencia es que aquí topo_sd_elev tiene una pendiente común MÁS la
# desviación por género (genero:topo_sd_elev), de modo que el término de interacción tiene 3 gl.
# Salida: datos/modelos/h2_wald_controles_por_genero_{tam}km.csv
# No probado (no hay R en mi entorno). Requiere: mgcv, dplyr, readr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv); library(dplyr); library(readr)

dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
tamanos_km <- c(50, 100, 200)
frugivoros <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
n_min_modelo <- 25; k_por_celdas <- 8; h_conf <- "topo_sd_elev"
hip  <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec")
ctrl <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2", "clima_bio01_media", "clima_bio12_media", "hab_evi_media")
vars <- c(hip, ctrl)

zscore <- function(x) as.numeric(scale(x))
preparar <- function(d, vars) {
  gens <- sort(unique(d$genero))
  d <- d %>% group_by(genero) %>% mutate(across(all_of(vars), zscore)) %>% ungroup()
  for (g in gens) d[[paste0("I_", g)]] <- as.numeric(d$genero == g)
  d$genero <- factor(d$genero); d
}
k_por_genero <- function(d) { n <- table(d$genero); setNames(pmax(5, pmin(40, floor(as.numeric(n) / k_por_celdas))), names(n)) }
# comunes = predictores con una sola pendiente; separados = predictores con pendiente por género
fml <- function(d, comunes, separados) {
  k <- k_por_genero(d)
  esp <- sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k), as.integer(k))
  as.formula(paste("cbind(S, S_fallo) ~", paste(c("genero", comunes, paste0("genero:", separados), esp), collapse = " + ")))
}
ajusta <- function(f, d, metodo) tryCatch(gam(f, family = binomial(), data = d, method = metodo), error = function(e) NULL)

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  base <- bind_rows(lapply(frugivoros, function(g) {
    d <- base %>% filter(genero == g)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S; d
  }))
  d0 <- base %>% filter(if_all(all_of(vars), ~ !is.na(.x)))
  ok <- names(which(table(d0$genero) >= n_min_modelo))
  d <- preparar(d0 %>% filter(genero %in% ok), vars)

  todo_gen  <- c(ctrl, hip)                         # todo con pendiente por género
  sin_topo  <- c(ctrl, setdiff(hip, h_conf))
  m_wald <- ajusta(fml(d, h_conf, todo_gen), d, "REML")   # topo común + genero:topo (3 gl) y el resto por género
  mf <- ajusta(fml(d, character(), todo_gen), d, "ML")    # topo por género
  mn <- ajusta(fml(d, h_conf, sin_topo), d, "ML")         # topo común
  p_wald <- NA_real_; gl_wald <- NA_real_
  at <- tryCatch(anova(m_wald)$pTerms.table, error = function(e) NULL)
  fila <- if (is.null(at)) integer() else which(rownames(at) == paste0("genero:", h_conf))
  if (length(fila) == 1) { p_wald <- at[fila, "p-value"]; gl_wald <- at[fila, "df"] }
  a <- tryCatch(anova(mn, mf, test = "Chisq"), error = function(e) NULL)
  p_lr <- if (!is.null(a) && !is.na(a$Deviance[2]) && a$Deviance[2] > 0) a[["Pr(>Chi)"]][2] else NA_real_
  res <- tibble(tam_km = tam_km, n = nrow(d), generos = paste(ok, collapse = ","),
                gl_wald = gl_wald, p_wald = p_wald, p_LR = p_lr, dAIC = AIC(mf) - AIC(mn))
  write_csv(res, file.path(dir_mod, sprintf("h2_wald_controles_por_genero_%dkm.csv", tam_km)))
  print(as.data.frame(res), row.names = FALSE)
}
