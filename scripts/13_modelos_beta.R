# 13_modelos_beta.R
# Modelos aditivos (GAM) de la diversidad beta por género (salida de 12_diversidad_beta.R) contra los predictores de
# heterogeneidad ambiental e inestabilidad climática histórica (base_modelos_{tam}km.csv, script 09).
#   - Recambio taxonómico (beta_tax, proporción 0-1 con muchos ceros a 500 km): quasibinomial, enlace logit.
#   - Recambio filogenético (beta_filo_ses, SES): gaussiano.
#   Predictores: modelo reducido (topo_sd_elev, clima_bio12_sd, hist_sd_temp, clima_bio01_media, ctrl_esfuerzo_log),
#   estandarizados dentro de cada género, + log_S (el recambio depende de la riqueza) + suavizado espacial s(x, y).
#   Radio principal: radio_principal_km (500); el radio Inf se corre como comparación con el borrador.
#
# Salidas (datos/modelos): coeficientes_beta_{tam}km.csv, diagnosticos_beta_{tam}km.csv
# Qué mirar: moran_p (autocorrelación residual; las medias contra vecinas se traslapan en el espacio, así que es esperable
#   que haya algo), dispersion, y que el signo de topo_sd_elev y hist_sd_temp sea parecido entre radios y géneros.
#
# No probado: no tengo R en mi entorno. Requiere: mgcv, ape, dplyr, readr
source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv); library(ape); library(dplyr); library(readr)

dir_base <- RUTA_BASE
dir_mod  <- file.path(dir_base, "datos", "modelos")
etiqueta <- "completo"
todas_las_escalas <- TRUE    # TRUE = corre 100, 50 y 200 km (en ese orden); FALSE = solo 100 km
tamanos_km <- if (todas_las_escalas) c(100, 50, 200) else c(100)
generos <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
radios_modelo_km <- c(500, Inf)
n_min_modelo <- 25; k_por_celdas <- 8
esfuerzo_var <- "ctrl_esfuerzo_log"
preds <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp", "clima_bio01_media", esfuerzo_var, "log_S")
pred_hipotesis <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp")
respuestas <- list(beta_tax = quasibinomial(), beta_filo_ses = gaussian())

zscore <- function(x) as.numeric(scale(x))
ajustar <- function(d, y, fam) {
  d[preds] <- lapply(d[preds], zscore)
  k_esp <- max(5, min(60, floor(nrow(d) / k_por_celdas)))
  f <- as.formula(paste0(y, " ~ ", paste(preds, collapse = " + "),
                         sprintf(" + s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", k_esp)))
  m <- tryCatch(gam(f, family = fam, data = d, method = "REML"),
                error = function(e) gam(f, family = fam, data = d, method = "GCV.Cp"))
  res <- residuals(m, type = "pearson")
  w <- 1 / as.matrix(dist(d[, c("x_centro_km", "y_centro_km")])); diag(w) <- 0
  mi <- Moran.I(res, w)
  cc <- tryCatch(concurvity(m, full = TRUE)["worst", ], error = function(e) NA)
  conc <- if (all(is.na(cc))) NA_real_ else unname(cc[grepl("^s\\(", names(cc))])[1]
  tab <- as.data.frame(summary(m)$p.table); names(tab) <- c("estimate", "se", "stat", "p")
  tab$term <- rownames(tab)
  tab <- tab %>% filter(term != "(Intercept)") %>% mutate(ic_inf = estimate - 1.96 * se, ic_sup = estimate + 1.96 * se)
  list(coef = tab, diag = tibble(n = nrow(d), k_espacial = k_esp, dev_expl = summary(m)$dev.expl,
         dispersion = sum(res^2) / df.residual(m), moran_I = mi$observed, moran_p = mi$p.value, conc_espacial = conc))
}

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  beta <- read_csv(file.path(dir_mod, sprintf("beta_%s_%dkm.csv", etiqueta, tam_km)), show_col_types = FALSE)
  coefs <- list(); diags <- list()
  for (g in generos) for (r in radios_modelo_km) for (y in names(respuestas)) {
    d <- beta %>% filter(genero == g, radio_km == r) %>% select(celda, all_of(y)) %>%
      inner_join(base %>% filter(genero == g) %>% select(celda, x_centro_km, y_centro_km, all_of(setdiff(preds, "log_S")), log_S),
                 by = "celda") %>% filter(!is.na(.data[[y]]))
    if (nrow(d) < n_min_modelo) { cat(sprintf("  [%s, %s km, %s]: solo %d celdas; se omite\n", g, r, y, nrow(d))); next }
    # beta_tax en [0, 1]; se asegura el rango
    if (y == "beta_tax") d[[y]] <- pmin(pmax(d[[y]], 0), 1)
    if (sd(d[[y]]) == 0) { cat(sprintf("  [%s, %s km, %s]: respuesta constante; se omite\n", g, r, y)); next }
    a <- ajustar(d, y, respuestas[[y]]); k <- paste(g, r, y)
    coefs[[k]] <- a$coef %>% mutate(genero = g, radio_km = r, respuesta = y)
    diags[[k]] <- a$diag %>% mutate(genero = g, radio_km = r, respuesta = y)
  }
  coefs <- bind_rows(coefs) %>% mutate(tam_km = tam_km) %>%
    select(tam_km, genero, radio_km, respuesta, term, estimate, se, stat, p, ic_inf, ic_sup)
  diags <- bind_rows(diags) %>% mutate(tam_km = tam_km) %>% select(tam_km, genero, radio_km, respuesta, everything())
  write_csv(coefs, file.path(dir_mod, sprintf("coeficientes_beta_%dkm.csv", tam_km)))
  write_csv(diags, file.path(dir_mod, sprintf("diagnosticos_beta_%dkm.csv", tam_km)))
  cat("\nDiagnósticos:\n")
  print(as.data.frame(diags %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
  cat("\nPredictores de hipótesis: estimado (error estándar) y p\n")
  print(as.data.frame(coefs %>% filter(term %in% pred_hipotesis) %>%
          mutate(estimado = sprintf("%.2f (%.2f)", estimate, se), p = round(p, 3)) %>%
          select(genero, radio_km, respuesta, term, estimado, p)), row.names = FALSE)
}
