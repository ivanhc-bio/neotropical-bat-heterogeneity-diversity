# =============================================================================
# 31b_resumen.R
# Solo imprime los resultados que el script 31 ya escribió en disco.
# Usa dplyr:: explícito porque terra enmascara select().
# =============================================================================

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(readr); library(dplyr); library(tidyr)

dir_mod <- file.path(RUTA_BASE, "datos", "modelos")

res  <- read_csv(file.path(dir_mod, "beta_precipitacion.csv"), show_col_types = FALSE)
pend <- read_csv(file.path(dir_mod, "beta_precipitacion_generos.csv"), show_col_types = FALSE)

cat("\n=== PRIMARIA (bloques de 500 km) ===\n")
print(as.data.frame(
  res |>
    dplyr::filter(agrupamiento == "bloque") |>
    dplyr::mutate(resumen = sprintf("%+.3f [%+.3f, %+.3f] (p=%.3f) | int p=%.4f",
                                    est_comun, ic_inf, ic_sup, p_comun, p_int)) |>
    dplyr::select(tam_km, escenario, respuesta, tipo, resumen) |>
    dplyr::arrange(escenario, respuesta, tipo, tam_km)), row.names = FALSE)

cat("\n\n=== Bloques frente a celda ===\n")
print(as.data.frame(
  res |>
    dplyr::select(tam_km, escenario, respuesta, tipo, agrupamiento, est_comun, p_comun) |>
    tidyr::pivot_wider(names_from = agrupamiento,
                       values_from = c(est_comun, p_comun)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4))) |>
    dplyr::arrange(escenario, respuesta, tipo, tam_km)), row.names = FALSE)

cat("\n\n=== Pendientes por género (bloques) ===\n")
print(as.data.frame(
  pend |>
    dplyr::filter(agrupamiento == "bloque") |>
    dplyr::mutate(pend = sprintf("%+.3f (p=%.3f)", estimate, p)) |>
    dplyr::select(tam_km, escenario, respuesta, predictor, genero, pend) |>
    dplyr::arrange(escenario, respuesta, predictor, tam_km, genero)), row.names = FALSE)
