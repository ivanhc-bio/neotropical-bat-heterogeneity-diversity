# =============================================================================
# 28_precipitacion.R
#
# Dos cosas:
#
# (a) Correlaciones entre las medidas de heterogeneidad, a nivel de CELDA
#     (no de celda-género). Sirve para confirmar si la desviación estándar de
#     la temperatura dentro de la celda es, en la práctica, topografía.
#
# (b) Análisis EXPLORATORIO del efecto de la precipitación sobre la riqueza y
#     el SES-MPD, con la misma estructura de modelo que el script 20.
#     Dos predictores de lluvia:
#        actual    = coeficiente de variación de bio12 dentro de la celda
#                    (clima_bio12_sd / clima_bio12_media; relativo, no absoluto)
#        historica = hist_sd_prec_log  (SD temporal de log(P+1) en 22 ka)
#     Cada modelo se corre dos veces:
#        sin_topografia = controles estándar
#        con_topografia = controles estándar + topo_sd_elev
#     La segunda dice qué añade la lluvia por encima del bloque topografía-
#     temperatura, que es la pregunta que interesa.
#
# ESTO ES EXPLORATORIO. No entra en la familia de 16 pruebas con Holm, que
# quedó fijada de antemano. Se reporta con tamaños de efecto e intervalos.
# Con 24 contrastes sin corregir, se esperan ~1.2 significativos por azar.
#
# Salidas (datos/modelos):
#   correlaciones_heterogeneidad.csv
#   precipitacion_exploratorio.csv
#   precipitacion_pendientes_genero.csv
#
# Ponlo junto a 20_alfa_pruebas.R. Requiere: mgcv, dplyr, readr, tidyr.
# =============================================================================

solo_funciones <- TRUE
source("20_alfa_pruebas.R")

suppressPackageStartupMessages({library(tidyr)})

tamanos_km <- c(50, 100, 200)

hip_prec <- c(actual = "clima_bio12_cv", historica = "hist_sd_prec_log")

vars_cor <- c("topo_sd_elev", "topo_rango_elev", "topo_tri_media", "topo_elev_media",
              "clima_bio01_sd", "clima_bio01_sd_res", "clima_bio01_media",
              "clima_bio12_sd", "clima_bio12_cv", "clima_bio12_media",
              "hist_sd_temp", "hist_anom_temp",
              "hist_anom_prec", "hist_sd_prec_log", "hist_sd_prec_mm")

anade_cv <- function(d) {
  if (all(c("clima_bio12_sd", "clima_bio12_media") %in% names(d)))
    d$clima_bio12_cv <- d$clima_bio12_sd / d$clima_bio12_media
  d
}

cor_todo <- list(); res_todo <- list(); pend_todo <- list()

for (tam_km in tamanos_km) {

  cat("\n\n=====================  ", tam_km, " km  =====================\n", sep = "")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)),
                   show_col_types = FALSE) %>% anade_cv()
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)),
                   show_col_types = FALSE)

  # --------------------------------------------------------------------------
  # (a) correlaciones, una fila por celda
  # --------------------------------------------------------------------------
  vv <- intersect(vars_cor, names(base))
  celdas <- base %>% distinct(celda, .keep_all = TRUE) %>% select(all_of(vv))

  cat("\nCorrelaciones sobre", nrow(celdas), "celdas distintas.\n")
  cm <- cor(celdas, use = "pairwise.complete.obs")

  if (all(c("topo_sd_elev", "clima_bio01_sd") %in% vv))
    cat("  *** topo_sd_elev vs clima_bio01_sd (SD de temperatura):  r =",
        round(cm["topo_sd_elev", "clima_bio01_sd"], 3), "***\n")
  for (v in intersect(c("clima_bio12_sd", "clima_bio12_cv", "hist_sd_prec_log",
                        "hist_anom_prec", "hist_sd_temp"), vv))
    cat("      topo_sd_elev vs ", v, ":  r = ",
        round(cm["topo_sd_elev", v], 3), "\n", sep = "")

  cor_todo[[as.character(tam_km)]] <- as.data.frame(cm) %>%
    tibble::rownames_to_column("variable") %>%
    pivot_longer(-variable, names_to = "con", values_to = "r") %>%
    mutate(tam_km = tam_km, n_celdas = nrow(celdas)) %>%
    filter(variable < con)

  # --------------------------------------------------------------------------
  # (b) modelos exploratorios de precipitación
  # --------------------------------------------------------------------------
  if (!all(unname(hip_prec) %in% names(base))) {
    message("Faltan ", paste(setdiff(unname(hip_prec), names(base)), collapse = ", "),
            "; se omiten los modelos de ", tam_km, " km.")
    next
  }

  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g); if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>%
      summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    d
  }))

  for (modo in c("sin_topografia", "con_topografia")) {

    extra <- if (modo == "con_topografia") "topo_sd_elev" else character()
    cat("\n--", modo, "--\n")

    # riqueza
    vr <- c(unname(hip_prec), ctrl_riq, extra)
    d0 <- base %>% filter(if_all(all_of(vr), ~ !is.na(.x)))
    ok <- names(which(table(d0$genero) >= n_min_modelo))
    dr <- preparar(d0 %>% filter(genero %in% ok), vr)
    cat("   riqueza: n =", nrow(dr), "\n")
    rr <- alfa_respuesta(dr, "cbind(S, S_fallo)", binomial(), c(ctrl_riq, extra),
                         "riqueza", tam_km, hip_prec)

    # SES-MPD
    vs <- c(unname(hip_prec), ctrl_ses, extra)
    d1 <- base %>% filter(filo_valida, !is.na(ses_mpd), if_all(all_of(vs), ~ !is.na(.x)))
    ok <- names(which(table(d1$genero) >= n_min_modelo))
    ds <- preparar(d1 %>% filter(genero %in% ok), vs)
    cat("   SES-MPD: n =", nrow(ds), "\n")
    rs <- alfa_respuesta(ds, "ses_mpd", gaussian(), c(ctrl_ses, extra),
                         "ses_mpd", tam_km, hip_prec)

    res_todo[[paste(tam_km, modo)]] <-
      bind_rows(rr$pruebas, rs$pruebas) %>% mutate(modelo = modo)
    pend_todo[[paste(tam_km, modo)]] <-
      bind_rows(rr$pendientes, rs$pendientes) %>% mutate(modelo = modo)
  }
}

cor_todo  <- bind_rows(cor_todo)
res_todo  <- bind_rows(res_todo)
pend_todo <- bind_rows(pend_todo)

write_csv(cor_todo,  file.path(dir_mod, "correlaciones_heterogeneidad.csv"))
write_csv(res_todo,  file.path(dir_mod, "precipitacion_exploratorio.csv"))
write_csv(pend_todo, file.path(dir_mod, "precipitacion_pendientes_genero.csv"))

cat("\n\n=== Pendiente común [IC] (p) | interacción p ===\n")
print(as.data.frame(
  res_todo %>%
    mutate(resumen = sprintf("%.3f [%.3f, %.3f] (p=%.3f) | int p=%.4f",
                             est_comun, ic_inf, ic_sup, p_comun, p_int)) %>%
    select(tam_km, modelo, respuesta, tipo, predictor, n, resumen) %>%
    arrange(respuesta, tipo, tam_km, modelo)), row.names = FALSE)

cat("\nRecuerda: exploratorio, sin corrección por multiplicidad.\n")
