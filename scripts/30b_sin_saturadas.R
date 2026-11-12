# =============================================================================
# 30b_sin_saturadas.R
# Prueba precisa de la saturación: reajusta los modelos de precipitación
# excluyendo exactamente las celdas con píxeles topados que identificó el
# script 30a, en lugar de aproximarlas por la media de lluvia.
#
# Tres escenarios:
#   completo        todas las celdas
#   sin_topadas     fuera las celdas con algún píxel en el tope
#   sin_topadas_10  fuera las celdas con más del 10 % de píxeles topados
#
# Salidas (datos/modelos):
#   precipitacion_sin_saturadas.csv
#   precipitacion_sin_saturadas_generos.csv
#
# Ponlo junto a 20_alfa_pruebas.R. Requiere haber corrido 30a.
# =============================================================================

solo_funciones <- TRUE
source("20_alfa_pruebas.R")

tamanos_km <- c(50, 100, 200)
hip_prec   <- c(actual = "clima_bio12_cv", historica = "hist_sd_prec_log")

res_todo <- list(); pend_todo <- list()

for (tam_km in tamanos_km) {

  cat("\n\n=====================  ", tam_km, " km  =====================\n", sep = "")

  ruta_sat <- file.path(dir_mod, sprintf("celdas_saturadas_%dkm.csv", tam_km))
  if (!file.exists(ruta_sat)) {
    message("No existe ", basename(ruta_sat), "; corre antes el script 30a. Se omite.")
    next
  }
  sat <- read_csv(ruta_sat, show_col_types = FALSE)
  cat("Celdas con algún píxel topado:", nrow(sat),
      "| con más del 10 %:", sum(sat$frac > 0.1), "\n")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)),
                   show_col_types = FALSE) %>%
    mutate(clima_bio12_cv = clima_bio12_sd / clima_bio12_media)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)),
                   show_col_types = FALSE)

  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g); if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>%
      summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    d
  }))

  # cuántas celdas-género se pierden en cada escenario
  for (nm in c("completo", "sin_topadas", "sin_topadas_10")) {
    fuera <- switch(nm,
                    completo       = character(0),
                    sin_topadas    = sat$celda,
                    sin_topadas_10 = sat$celda[sat$frac > 0.1])
    cat("\n--", nm, "-- celdas excluidas:", length(fuera), "\n")

    bb <- base %>% filter(!celda %in% fuera)
    extra <- "topo_sd_elev"

    vr <- c(unname(hip_prec), ctrl_riq, extra)
    d0 <- bb %>% filter(if_all(all_of(vr), ~ !is.na(.x)))
    ok <- names(which(table(d0$genero) >= n_min_modelo))
    if (length(ok) < 2) { cat("   pocos géneros; se omite\n"); next }
    dr <- preparar(d0 %>% filter(genero %in% ok), vr)
    cat("   riqueza: n =", nrow(dr), "\n")
    rr <- alfa_respuesta(dr, "cbind(S, S_fallo)", binomial(), c(ctrl_riq, extra),
                         "riqueza", tam_km, hip_prec)

    vs <- c(unname(hip_prec), ctrl_ses, extra)
    d1 <- bb %>% filter(filo_valida, !is.na(ses_mpd), if_all(all_of(vs), ~ !is.na(.x)))
    ok <- names(which(table(d1$genero) >= n_min_modelo))
    rs <- NULL
    if (length(ok) >= 2) {
      ds <- preparar(d1 %>% filter(genero %in% ok), vs)
      cat("   SES-MPD: n =", nrow(ds), "\n")
      rs <- alfa_respuesta(ds, "ses_mpd", gaussian(), c(ctrl_ses, extra),
                           "ses_mpd", tam_km, hip_prec)
    } else cat("   SES-MPD: pocos géneros; se omite\n")

    res_todo[[paste(tam_km, nm)]]  <- bind_rows(rr$pruebas, rs$pruebas) %>%
      mutate(escenario = nm, n_excluidas = length(fuera))
    pend_todo[[paste(tam_km, nm)]] <- bind_rows(rr$pendientes, rs$pendientes) %>%
      mutate(escenario = nm)
  }
}

res_todo  <- bind_rows(res_todo)
pend_todo <- bind_rows(pend_todo)
write_csv(res_todo,  file.path(dir_mod, "precipitacion_sin_saturadas.csv"))
write_csv(pend_todo, file.path(dir_mod, "precipitacion_sin_saturadas_generos.csv"))

cat("\n\n=== Interacción género x lluvia actual ===\n")
print(as.data.frame(
  res_todo %>% filter(predictor == "clima_bio12_cv") %>%
    mutate(resumen = sprintf("comun %+.3f (p=%.3f) | int p=%.4f",
                             est_comun, p_comun, p_int)) %>%
    select(tam_km, respuesta, escenario, n, resumen) %>%
    arrange(respuesta, tam_km, escenario)), row.names = FALSE)

cat("\n=== Pendientes por género, lluvia actual ===\n")
print(as.data.frame(
  pend_todo %>% filter(predictor == "clima_bio12_cv") %>%
    mutate(pend = sprintf("%+.3f (p=%.3f)", estimate, p)) %>%
    select(tam_km, respuesta, genero, escenario, pend) %>%
    arrange(respuesta, tam_km, genero, escenario)), row.names = FALSE)
