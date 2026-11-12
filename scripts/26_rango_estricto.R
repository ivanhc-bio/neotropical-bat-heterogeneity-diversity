# =============================================================================
# 26_rango_estricto.R
# Arma los archivos de entrada de la variante "rango estricto": sin ningún
# registro a más de 100 km del polígono de rango de su especie.
#
# Salidas:
#   datos/rejilla/presencia_{tam}km_rango_estricto.csv
#   datos/modelos/base_modelos_{tam}km_rango_estricto.csv
#
# Después de correr esto, el script 25 toma la variante solo (ya la tiene
# programada como "rango_estricto") y la añade a alfa_robustez_resumen.csv.
#
# Importante: esta variante solo vale para la RIQUEZA. El SES-MPD tendría que
# recalcularse con el modelo nulo y los árboles, así que las columnas de SES se
# dejan en NA y filo_valida en FALSE, para que el script 25 no las modele.
#
# Las celdas incluidas no cambian: se mantiene la misma lista de celdas y solo
# se recalculan S y los registros por género.
# =============================================================================

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
})
if (!requireNamespace("sf", quietly = TRUE))
  stop("Falta el paquete sf. Instálalo con install.packages('sf') y vuelve a correr.")

base    <- RUTA_BASE
dir_rej <- file.path(base, "datos", "rejilla")
dir_mod <- file.path(base, "datos", "modelos")
tamanos_km <- c(50, 100, 200)

# misma proyección que la rejilla
crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# ---- 1. Registros que la variante estricta elimina --------------------------

fuera <- read_csv(file.path(dir_rej, "registros_fuera_rango.csv"), show_col_types = FALSE)
cat("Registros más allá del margen de 100 km:", nrow(fuera), "\n\n")
print(as.data.frame(count(fuera, estado_rango, decision_rango)), row.names = FALSE)

# los que ya se habían eliminado no se vuelven a restar
ya_fuera <- str_detect(tolower(fuera$decision_rango), "excluir")
cat("\nYa eliminados en la curación:", sum(ya_fuera),
    "\nQue la variante estricta elimina de más:", sum(!ya_fuera), "\n")

quitar <- fuera[!ya_fuera, ]

p  <- sf::st_as_sf(quitar, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
xy <- sf::st_coordinates(sf::st_transform(p, crs_laea)) / 1000
quitar <- quitar %>% mutate(x_km = xy[, 1], y_km = xy[, 2])

# ---- 2. Un tamaño de celda a la vez -----------------------------------------

for (tam in tamanos_km) {

  cat("\n\n=====================  ", tam, " km  =====================\n", sep = "")

  ruta_pres <- file.path(dir_rej, sprintf("presencia_%dkm.csv", tam))
  ruta_bm   <- file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam))
  if (!file.exists(ruta_pres) || !file.exists(ruta_bm)) {
    message("Falta ", basename(ruta_pres), " o ", basename(ruta_bm), "; se omite.")
    next
  }

  pres <- read_csv(ruta_pres, show_col_types = FALSE)
  bm   <- read_csv(ruta_bm,   show_col_types = FALSE)

  q <- quitar %>%
    mutate(celda = sprintf("c%d_%d_%d", tam, floor(x_km / tam), floor(y_km / tam))) %>%
    count(celda, especie, name = "n_quitar")

  # comprobación de la decodificación de la rejilla
  emparejadas <- sum(q$n_quitar[paste(q$celda, q$especie) %in%
                                  paste(pres$celda, pres$especie)])
  cat("Registros a quitar que caen en una fila de presencia: ",
      emparejadas, " de ", sum(q$n_quitar),
      " (", round(100 * emparejadas / sum(q$n_quitar), 1), " %)\n", sep = "")
  if (emparejadas < 0.5 * sum(q$n_quitar))
    warning("Muy pocos emparejan: revisa la proyección o el formato del nombre de celda.")

  # ---- presencia estricta ---------------------------------------------------
  nueva <- pres %>%
    left_join(q, by = c("celda", "especie")) %>%
    mutate(n_quitar = coalesce(n_quitar, 0),
           n_reg    = n_reg - n_quitar)

  negativos <- sum(nueva$n_reg < 0)
  if (negativos > 0)
    cat("AVISO:", negativos, "filas con n_reg negativo; se truncan en 0.\n")

  nueva <- nueva %>%
    mutate(n_reg = pmax(n_reg, 0)) %>%
    filter(n_reg > 0) %>%
    select(all_of(names(pres)))

  cat("Combinaciones celda-especie: ", nrow(pres), " -> ", nrow(nueva),
      " (", round(100 * nrow(nueva) / nrow(pres), 1), " % conservadas)\n", sep = "")

  write_csv(nueva, file.path(dir_rej, sprintf("presencia_%dkm_rango_estricto.csv", tam)))

  # ---- base de modelos estricta ---------------------------------------------
  nuevo_S <- nueva %>%
    group_by(celda, genero = genero_analisis) %>%
    summarise(S_e = n_distinct(especie), n_reg_genero_e = sum(n_reg), .groups = "drop")

  col_ses <- intersect(
    c("mpd_obs", "pd_obs", "ses_mpd", "ses_pd", "ses_mpd_de_arboles", "ses_pd_de_arboles",
      "rango_mpd", "rango_pd", "S_ses", "ses_mpd_sin_hermanas", "ses_pd_sin_hermanas",
      "ses_mpd_solo_adn", "ses_pd_solo_adn"),
    names(bm))

  bm_e <- bm %>%
    left_join(nuevo_S, by = c("celda", "genero")) %>%
    mutate(S_e            = coalesce(S_e, 0L),
           n_reg_genero_e = coalesce(n_reg_genero_e, 0),
           S              = S_e,
           n_reg_genero   = n_reg_genero_e,
           log_S          = ifelse(S > 0, log(S), NA_real_),
           esfuerzo_genero_log = ifelse(n_reg_genero > 0, log(n_reg_genero), NA_real_),
           across(all_of(col_ses), ~ NA_real_)) %>%
    select(-S_e, -n_reg_genero_e)

  if ("filo_valida" %in% names(bm_e)) bm_e$filo_valida <- FALSE

  cat("Celdas-género: ", nrow(bm_e), " (sin cambio). S media ",
      round(mean(bm$S), 3), " -> ", round(mean(bm_e$S), 3),
      "; celdas-género que pierden alguna especie: ", sum(bm_e$S < bm$S), "\n", sep = "")

  write_csv(bm_e, file.path(dir_mod, sprintf("base_modelos_%dkm_rango_estricto.csv", tam)))
}

cat("\n\nListo. Ahora corre el script 25 con:\n",
    "  tamanos_km <- c(50, 100, 200)\n",
    "  solo_variantes <- NULL\n",
    "y busca la variante 'rango_estricto' en alfa_robustez_resumen.csv.\n",
    "Solo las filas de riqueza son válidas para esa variante.\n", sep = "")
