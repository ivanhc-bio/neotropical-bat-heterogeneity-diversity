# 23_anoura.R
# Contraste exploratorio con Anoura (género nectarívoro): las mismas 8 combinaciones de predictor y respuesta, pero solo con Anoura,
# comparadas con la pendiente común de los 4 frugívoros. Queda FUERA de la corrección de Holm (exploratorio).
#   Alfa : riqueza (binomial) y SES-MPD (gaussiano) ~ topo_sd_elev + hist_sd_temp + controles + s(x, y)   (GAM, REML)
#   Beta : recambio taxonómico y filogenético (cuasibinomial) entre pares de celdas de Anoura ~ |Δ elevación media| + |Δ hist_sd_temp|
#          + controles + distancia (spline), con covarianza robusta agrupada por BLOQUES de 500 km (igual que beta principal).
# Contraste: z = (b_Anoura - b_frugívoros) / sqrt(se_Anoura^2 + se_frugívoros^2); supone que las dos estimaciones son independientes (Anoura
#   no entra en la estimación de los frugívoros, y las celdas comparten solo en parte): es una aproximación, no una prueba exacta.
# Lee: alfa_pruebas_{tam}km.csv (script 20) y beta_robustez_{tam}km.csv, variante bloque_s1 (script 24).
# Salida (datos/modelos): anoura_contraste_{tam}km.csv
# No probado: no hay R en mi entorno. Requiere los scripts 20 y 21 en la misma carpeta.
# Con pocas celdas (SES-MPD con >= 3 especies y beta con >= 2 especies), algunos modelos pueden no ajustarse: el script lo avisa y sigue.

solo_funciones <- TRUE
source("20_alfa_pruebas.R")                       # ajusta, preparar, k_por_genero, ...
hip_alfa <- c(presente = "topo_sd_elev", historica = "hist_sd_temp")
ctrl_riq_a <- ctrl_riq; ctrl_ses_a <- ctrl_ses
source("21_beta_pruebas.R")                       # pares_tamano, a_bloques, ajuste_diadico, pendiente, ...
hip_beta <- c(presente = "d_topo_elev_media", historica = "d_hist_sd_temp")
ctrl_beta <- c("d_bio01", "d_bio12", "d_esf", "d_logS", "logS_min")

tamanos_km <- c(100, 50, 200)
n_min_anoura <- 25
bloque_km <- 500

# ---- alfa: Anoura solo ------------------------------------------------------------------------------------------------------
alfa_anoura <- function(d, y, fam, ctrl, etiqueta, tam_km) {
  h <- unname(hip_alfa)
  d[c(h, ctrl)] <- lapply(d[c(h, ctrl)], function(x) as.numeric(scale(x)))
  k <- max(5, min(40, floor(nrow(d) / k_por_celdas)))
  f <- as.formula(paste(y, "~", paste(c(ctrl, h, sprintf("s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", k)), collapse = " + ")))
  m <- ajusta(f, d, fam)
  tab <- as.data.frame(summary(m)$p.table)
  bind_rows(lapply(names(hip_alfa), function(tipo) { p <- hip_alfa[[tipo]]
    tibble(faceta = "alfa", respuesta = etiqueta, tipo = tipo, predictor = p, n = nrow(d), est_anoura = tab[p, 1], se_anoura = tab[p, 2],
           ic_inf = tab[p, 1] - 1.96 * tab[p, 2], ic_sup = tab[p, 1] + 1.96 * tab[p, 2], p_anoura = tab[p, 4]) }))
}

# ---- beta: Anoura solo ------------------------------------------------------------------------------------------------------
beta_anoura <- function(P, y, tam_km) {
  h <- unname(hip_beta)
  f <- as.formula(paste("~ ns(dist_km, 3) +", paste(c(ctrl_beta, h), collapse = " + ")))
  a <- ajuste_diadico(model.matrix(f, P), P[[y]], P$ci, P$cj)
  bind_rows(lapply(names(hip_beta), function(tipo) { p <- hip_beta[[tipo]]; s <- pendiente(a, p)
    tibble(faceta = "beta", respuesta = y, tipo = tipo, predictor = p, n = nrow(P), est_anoura = s$est, se_anoura = s$se,
           ic_inf = s$ic[1], ic_sup = s$ic[2], p_anoura = s$p) }))
}

# ---- contraste con los frugívoros ------------------------------------------------------------------------------------------
contrastar <- function(anoura, ref) {
  anoura %>% left_join(ref, by = c("faceta", "respuesta", "tipo")) %>%
    mutate(z_dif = (est_anoura - est_frug) / sqrt(se_anoura^2 + se_frug^2), p_dif = 2 * pnorm(-abs(z_dif)))
}

arboles <- readRDS(archivo_arboles); arboles <- arboles[seq_len(min(n_arboles_beta, length(arboles)))]
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  res <- list()

  # alfa
  pool <- pres %>% filter(genero_analisis == "Anoura", celda %in% base$celda[base$genero == "Anoura"]) %>% summarise(n_distinct(especie)) %>% pull()
  d <- base %>% filter(genero == "Anoura"); d$S_fallo <- pool - d$S
  vr <- c(unname(hip_alfa), ctrl_riq_a); dr <- d %>% filter(if_all(all_of(vr), ~ !is.na(.x)))
  if (nrow(dr) >= n_min_anoura) res$riq <- tryCatch(alfa_anoura(dr, "cbind(S, S_fallo)", binomial(), ctrl_riq_a, "riqueza", tam_km),
                                                    error = function(e) { cat("  riqueza: falló:", conditionMessage(e), "\n"); NULL })
  vs <- c(unname(hip_alfa), ctrl_ses_a); ds <- d %>% filter(filo_valida, !is.na(ses_mpd), if_all(all_of(vs), ~ !is.na(.x)))
  cat(sprintf("  Anoura: riqueza n = %d; SES-MPD n = %d\n", nrow(dr), nrow(ds)))
  if (nrow(ds) >= n_min_anoura) res$ses <- tryCatch(alfa_anoura(ds, "ses_mpd", gaussian(), ctrl_ses_a, "ses_mpd", tam_km),
                                                    error = function(e) { cat("  SES-MPD: falló:", conditionMessage(e), "\n"); NULL })
  else cat("  SES-MPD: muy pocas celdas (< ", n_min_anoura, "); se omite\n", sep = "")

  # beta
  P <- pares_tamano(base, pres, arboles, tam_km, "Anoura", hip_beta)
  if (nrow(P) > 0) {
    P <- a_bloques(P, base, bloque_km)
    for (y in c("beta_tax", "beta_filo")) res[[y]] <- tryCatch(beta_anoura(P, y, tam_km),
                                                               error = function(e) { cat("  ", y, ": falló:", conditionMessage(e), "\n"); NULL })
  } else cat("  beta: Anoura no tiene celdas suficientes; se omite\n")

  anoura <- bind_rows(res)
  if (!nrow(anoura)) next
  # referencia: pendiente común de los frugívoros
  ra <- read_csv(file.path(dir_mod, sprintf("alfa_pruebas_%dkm.csv", tam_km)), show_col_types = FALSE) %>%
    transmute(faceta = "alfa", respuesta, tipo, est_frug = est_comun, se_frug = se_comun)
  rb <- read_csv(file.path(dir_mod, sprintf("beta_robustez_%dkm.csv", tam_km)), show_col_types = FALSE) %>% filter(variante == "bloque_s1") %>%
    transmute(faceta = "beta", respuesta, tipo, est_frug = est_comun, se_frug = (ic_sup - ic_inf) / (2 * 1.96))
  out <- contrastar(anoura, bind_rows(ra, rb)) %>% mutate(tam_km = tam_km) %>% select(tam_km, everything())
  write_csv(out, file.path(dir_mod, sprintf("anoura_contraste_%dkm.csv", tam_km)))
  cat("\nAnoura frente a los frugívoros (est = pendiente por DE; p_dif = prueba de la diferencia; exploratorio):\n")
  print(as.data.frame(out %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          select(respuesta, tipo, n, est_anoura, ic_inf, ic_sup, p_anoura, est_frug, p_dif)), row.names = FALSE)
}
