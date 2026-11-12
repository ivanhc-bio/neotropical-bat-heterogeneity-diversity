# 10_modelos.R
# Modelos aditivos (GAM) por género y tamaño de celda, a partir de base_modelos_{tam}km.csv (script 09):
#   - Riqueza observada del género (binomial: S de las especies del grupo de especies del género), modelo completo.
#   - SES-MPD (gaussiano): modelo completo (Artibeus, Carollia) o reducido (Platyrrhinus, Sturnira, Anoura), con log(S)
#     como covariable en Carollia y Sturnira (el SES depende algo de la riqueza).
#   - Sensibilidades filogenéticas: SES-MPD sin hermanas, SES-MPD solo ADN y SES-PD (comprobación).
# Todos incluyen un suavizado espacial s(x, y) de dimensión baja para absorber la autocorrelación espacial residual y los
# predictores estandarizados (media 0, DE 1) como términos lineales, así que los coeficientes son comparables.
#
# Salidas (datos/modelos): coeficientes_modelos_{tam}km{sufijo}.csv y diagnosticos_modelos_{tam}km{sufijo}.csv,
# donde sufijo = "_esfuerzo_genero" si esfuerzo_var = "esfuerzo_genero_log" y "_sin_espacial" si con_espacial = FALSE.
#
# Qué mirar antes de interpretar:
#   - dispersion (riqueza): si es mucho mayor que 1, la binomial subestima la varianza; prueba familia_riqueza = "quasibinomial".
#   - moran_I y moran_p de los residuos: si hay autocorrelación espacial residual, aumenta k_espacial o revisa el modelo.
#   - conc_espacial: concurvidad del suavizado espacial con el resto (cercana a 1 = el suavizado compite con los predictores).
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

tamanos_km      <- c(50, 100, 200)                     # luego 50 y 200 como sensibilidad
generos         <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
familia_riqueza <- "binomial"                 # "binomial", "quasibinomial" o "nb"
esfuerzo_var    <- "esfuerzo_genero_log"        # sensibilidad: "esfuerzo_genero_log" (registros del género)
genero_con_logS <- c("Carollia", "Sturnira")
n_min_modelo    <- 25                         # no se ajusta un modelo con menos celdas
k_por_celdas    <- 8                          # k del suavizado espacial = n / 8, entre 5 y 60
con_espacial    <- TRUE                       # FALSE = sin suavizado espacial (para comparar coeficientes y ver la autocorrelación)
# El nombre de los archivos indica el esfuerzo usado: sin marca = esfuerzo general (ctrl_esfuerzo_log);
# "_esfuerzo_genero" = registros del género. Así una corrida no sobrescribe a la otra.
sufijo          <- paste0(if (esfuerzo_var == "ctrl_esfuerzo_log") "" else "_esfuerzo_genero",
                          if (con_espacial) "" else "_sin_espacial")

pred_completo <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec",
                   "clima_bio01_media", "clima_bio12_media", "hab_evi_media", "ctrl_esfuerzo_log", "ctrl_area_tierra_km2")
pred_reducido <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp", "clima_bio01_media", "ctrl_esfuerzo_log")
pred_hipotesis <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec")

respuestas_filo <- c(ses_mpd = "principal", ses_mpd_sin_hermanas = "sens_sin_hermanas",
                     ses_mpd_solo_adn = "sens_solo_adn", ses_pd = "comprobacion_pd")

# ---- 1. Funciones -------------------------------------------------------------------
zscore <- function(x) as.numeric(scale(x))

ajustar <- function(d, y, preds, familia, extra = character()) {
  vars <- c(preds, extra)
  vars <- sub("^ctrl_esfuerzo_log$", esfuerzo_var, vars)
  d[vars] <- lapply(d[vars], zscore)
  k_esp <- max(5, min(60, floor(nrow(d) / k_por_celdas)))
  termino_esp <- if (con_espacial) sprintf(" + s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", k_esp) else ""
  f <- as.formula(paste0(y, " ~ ", paste(vars, collapse = " + "), termino_esp))
  fam <- switch(familia, binomial = binomial(), quasibinomial = quasibinomial(), nb = nb(), gaussian = gaussian())
  m <- tryCatch(gam(f, family = fam, data = d, method = "REML"),
                error = function(e) gam(f, family = fam, data = d, method = "GCV.Cp"))
  res <- residuals(m, type = "pearson")
  w <- 1 / as.matrix(dist(d[, c("x_centro_km", "y_centro_km")])); diag(w) <- 0
  mi <- Moran.I(res, w)
  cc <- tryCatch(concurvity(m, full = TRUE)["worst", ], error = function(e) NA)
  conc <- if (all(is.na(cc))) NA_real_ else unname(cc[grepl("^s\\(", names(cc))])[1]
  tab <- as.data.frame(summary(m)$p.table); names(tab) <- c("estimate", "se", "stat", "p")
  tab$term <- rownames(tab); tab <- tab %>% filter(term != "(Intercept)") %>%
    mutate(ic_inf = estimate - 1.96 * se, ic_sup = estimate + 1.96 * se)
  list(coef = tab,
       diag = tibble(n = nrow(d), k_espacial = k_esp, n_terminos = length(vars), dev_expl = summary(m)$dev.expl,
                     AIC = AIC(m), dispersion = sum(res^2) / df.residual(m), moran_I = mi$observed,
                     moran_p = mi$p.value, conc_espacial = conc))
}

# ---- 2. Modelos por tamaño de celda y género ----------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  coefs <- list(); diags <- list()

  for (g in generos) {
    d <- base %>% filter(genero == g)
    # grupo de especies del género: las que aparecen en alguna celda incluida
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    if (any(d$S_fallo < 0)) stop("Hay celdas con más especies que el grupo de especies (", g, ")")

    # ---- riqueza
    if (nrow(d) >= n_min_modelo) {
      r <- ajustar(d, "cbind(S, S_fallo)", pred_completo, familia_riqueza)
      coefs[[paste(g, "riqueza")]] <- r$coef %>% mutate(genero = g, respuesta = "riqueza", modelo = "completo", sensibilidad = "principal")
      diags[[paste(g, "riqueza")]] <- r$diag %>% mutate(genero = g, respuesta = "riqueza", modelo = "completo", sensibilidad = "principal", pool = pool)
    }

    # ---- SES (celdas con S >= 3)
    modelo_f <- unique(d$modelo_filo)
    preds_f <- if (modelo_f == "completo") pred_completo else pred_reducido
    extra_f <- if (g %in% genero_con_logS) "log_S" else character()
    for (y in names(respuestas_filo)) {
      dy <- d %>% filter(filo_valida, !is.na(.data[[y]]))
      if (nrow(dy) < n_min_modelo) { cat(sprintf("  [%s] %s: solo %d celdas; se omite\n", g, y, nrow(dy))); next }
      r <- ajustar(dy, y, preds_f, "gaussian", extra_f)
      coefs[[paste(g, y)]] <- r$coef %>% mutate(genero = g, respuesta = "ses_mpd_o_pd", modelo = modelo_f, sensibilidad = respuestas_filo[[y]])
      diags[[paste(g, y)]] <- r$diag %>% mutate(genero = g, respuesta = y, modelo = modelo_f, sensibilidad = respuestas_filo[[y]], pool = pool)
    }
  }

  coefs <- bind_rows(coefs) %>% mutate(tam_km = tam_km, esfuerzo = esfuerzo_var) %>%
    dplyr::select(tam_km, genero, respuesta, modelo, sensibilidad, esfuerzo, term, estimate, se, stat, p, ic_inf, ic_sup)
  diags <- bind_rows(diags) %>% mutate(tam_km = tam_km, esfuerzo = esfuerzo_var) %>%
    dplyr::select(tam_km, genero, respuesta, modelo, sensibilidad, esfuerzo, everything())
  write_csv(coefs, file.path(dir_mod, sprintf("coeficientes_modelos_%dkm%s.csv", tam_km, sufijo)))
  write_csv(diags, file.path(dir_mod, sprintf("diagnosticos_modelos_%dkm%s.csv", tam_km, sufijo)))

  cat("\nDiagnósticos:\n")
  print(as.data.frame(diags %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
                        dplyr::select(genero, respuesta, sensibilidad, n, k_espacial, dev_expl, dispersion, moran_I, moran_p, conc_espacial)),
        row.names = FALSE)
  cat("\nPredictores de hipótesis (modelo principal): estimado (error estándar) y p\n")
  print(as.data.frame(coefs %>% filter(sensibilidad == "principal", term %in% pred_hipotesis) %>%
                        mutate(estimado = sprintf("%.2f (%.2f)", estimate, se), p = round(p, 3)) %>%
                        dplyr::select(genero, respuesta, term, estimado, p)), row.names = FALSE)
}
