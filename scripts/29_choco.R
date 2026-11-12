# =============================================================================
# 29_choco.R
#
# ¿La señal de variabilidad de lluvia depende de la saturación de CHELSA?
#
# CHELSA v2.1 trunca bio12 en 6553.5 mm ("6553.5 o más", no dato faltante), así
# que en las celdas muy lluviosas (Chocó, Caribe de Costa Rica y Nicaragua) la
# media y sobre todo la SD de precipitación están subestimadas. Esas celdas caen
# justo en el territorio de Platyrrhinus y Sturnira, que son los dos géneros que
# muestran el efecto, así que hay que descartar que el patrón venga del artefacto.
#
# Se reajusta el modelo "con_topografia" del script 28 quitando las celdas más
# lluviosas con cuatro umbrales, y se compara la interacción género x lluvia y
# las pendientes por género con las de la corrida completa.
#
# Salidas (datos/modelos):
#   precipitacion_choco.csv           pruebas por umbral
#   precipitacion_choco_generos.csv   pendientes por género y umbral
#
# Ponlo junto a 20_alfa_pruebas.R.
# =============================================================================

solo_funciones <- TRUE
source("20_alfa_pruebas.R")

tamanos_km <- c(50, 100, 200)
hip_prec   <- c(actual = "clima_bio12_cv", historica = "hist_sd_prec_log")
TOPE       <- 6553.5          # tope de CHELSA v2.1 para bio12

res_todo <- list(); pend_todo <- list()

for (tam_km in tamanos_km) {

  cat("\n\n=====================  ", tam_km, " km  =====================\n", sep = "")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)),
                   show_col_types = FALSE) %>%
    mutate(clima_bio12_cv = clima_bio12_sd / clima_bio12_media)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)),
                   show_col_types = FALSE)

  celdas <- base %>% distinct(celda, .keep_all = TRUE)
  cat("Celdas distintas:", nrow(celdas),
      "\n  clima_bio12_media: máximo", round(max(celdas$clima_bio12_media, na.rm = TRUE)),
      "mm; tope de CHELSA", TOPE, "mm\n")
  for (q in c(0.90, 0.95, 0.99))
    cat("  percentil", q * 100, ":", round(quantile(celdas$clima_bio12_media, q, na.rm = TRUE)), "mm\n")
  cat("  celdas con media > 3000 mm:", sum(celdas$clima_bio12_media > 3000, na.rm = TRUE),
      "; > 4000 mm:", sum(celdas$clima_bio12_media > 4000, na.rm = TRUE), "\n")

  umbrales <- list(
    completo      = Inf,
    sin_p95       = unname(quantile(celdas$clima_bio12_media, 0.95, na.rm = TRUE)),
    sin_p90       = unname(quantile(celdas$clima_bio12_media, 0.90, na.rm = TRUE)),
    sin_mas_4000  = 4000,
    sin_mas_3000  = 3000)

  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g); if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>%
      summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    d
  }))

  for (nm in names(umbrales)) {

    lim <- umbrales[[nm]]
    bb  <- base %>% filter(is.na(clima_bio12_media) | clima_bio12_media <= lim)
    cat("\n--", nm, "(<=", round(lim), "mm ) --\n")

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

    res_todo[[paste(tam_km, nm)]] <-
      bind_rows(rr$pruebas, rs$pruebas) %>% mutate(umbral = nm, lim_mm = lim)
    pend_todo[[paste(tam_km, nm)]] <-
      bind_rows(rr$pendientes, rs$pendientes) %>% mutate(umbral = nm, lim_mm = lim)
  }
}

res_todo  <- bind_rows(res_todo)
pend_todo <- bind_rows(pend_todo)
write_csv(res_todo,  file.path(dir_mod, "precipitacion_choco.csv"))
write_csv(pend_todo, file.path(dir_mod, "precipitacion_choco_generos.csv"))

cat("\n\n=== Interacción género x lluvia actual, por umbral ===\n")
print(as.data.frame(
  res_todo %>% filter(predictor == "clima_bio12_cv") %>%
    mutate(resumen = sprintf("comun %.3f (p=%.3f) | int p=%.4f", est_comun, p_comun, p_int)) %>%
    select(tam_km, respuesta, umbral, n, resumen) %>%
    arrange(respuesta, tam_km, umbral)), row.names = FALSE)

cat("\n=== Pendientes por género, lluvia actual ===\n")
print(as.data.frame(
  pend_todo %>% filter(predictor == "clima_bio12_cv") %>%
    mutate(pend = sprintf("%+.3f (p=%.3f)", estimate, p)) %>%
    select(tam_km, respuesta, genero, umbral, pend) %>%
    arrange(respuesta, tam_km, genero, umbral)), row.names = FALSE)
