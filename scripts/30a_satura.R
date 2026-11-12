# =============================================================================
# 30a_satura.R
# ¿Cuánto muerde la saturación de CHELSA bio12 en la zona de estudio?
#
# Cuenta los píxeles que llegan al tope (6553.5 mm), ve en qué celdas caen y
# qué géneros ocupan esas celdas. Si son pocos y no tocan a las celdas que
# sostienen el resultado, no hace falta cambiar de capa.
#
# Solo lee. No escribe rásteres; deja un CSV con las celdas afectadas.
# =============================================================================

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

suppressPackageStartupMessages({library(terra); library(dplyr); library(readr)})

base    <- RUTA_BASE
dir_amb <- file.path(base, "datos", "ambiente")
dir_rej <- file.path(base, "datos", "rejilla")
dir_mod <- file.path(base, "datos", "modelos")

TOPE <- 6553.5

# ---- encontrar la capa de bio12 ---------------------------------------------
cand <- list.files(dir_amb, pattern = "bio12|bio_12|_12_", recursive = TRUE,
                   full.names = TRUE, ignore.case = TRUE)
cand <- cand[grepl("\\.(tif|tiff)$", cand, ignore.case = TRUE)]
cat("Archivos candidatos de bio12:\n"); print(basename(cand))
if (!length(cand)) stop("No encontré la capa de bio12 en ", dir_amb,
                        ". Pon su ruta a mano en la variable 'ruta'.")
ruta <- cand[1]
cat("\nUso:", ruta, "\n\n")

r <- rast(ruta)
cat("Resolución:", res(r), " | CRS:", crs(r, describe = TRUE)$name, "\n")

# ---- recorte a la zona de estudio -------------------------------------------
celdas <- vect(file.path(dir_rej, "celdas_100km.gpkg"))
zona   <- project(celdas, crs(r))
r <- crop(r, ext(zona))

v <- values(r, mat = FALSE)
v <- v[!is.na(v)]
cat("\nPíxeles con dato en la zona:", length(v),
    "\nMáximo:", round(max(v), 1), "mm\n")

for (u in c(TOPE, TOPE - 0.5, 6000, 5000, 4000)) {
  n <- sum(v >= u)
  cat(sprintf("  píxeles >= %6.1f mm: %8d  (%.4f %%)\n", u, n, 100 * n / length(v)))
}

if (max(v) < TOPE - 0.5) {
  cat("\n=> Ningún píxel llega al tope. La saturación NO afecta a esta zona.\n",
      "   No hace falta cambiar de capa.\n", sep = "")
} else {
  cat("\n=> Hay píxeles en el tope. Veamos en qué celdas caen.\n")

  topados <- ifel(r >= TOPE - 0.5, 1, 0)
  for (tam in c(50, 100, 200)) {
    cv <- vect(file.path(dir_rej, sprintf("celdas_%dkm.gpkg", tam)))
    cz <- project(cv, crs(r))
    n_top <- terra::extract(topados, cz, fun = sum, na.rm = TRUE)[, 2]
    n_tot <- terra::extract(!is.na(r), cz, fun = sum, na.rm = TRUE)[, 2]
    id <- names(cz)[grepl("celda", names(cz), ignore.case = TRUE)][1]
    d <- tibble(celda = as.data.frame(cz)[[id]], tam_km = tam,
                px_topados = n_top, px_total = n_tot) %>%
      mutate(frac = px_topados / pmax(px_total, 1)) %>%
      filter(px_topados > 0)
    cat("\n", tam, " km: ", nrow(d), " celdas con algún píxel topado; ",
        sum(d$frac > 0.1), " con más del 10 % topado.\n", sep = "")
    if (nrow(d)) {
      pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam)),
                       show_col_types = FALSE)
      g <- pres %>% filter(celda %in% d$celda) %>% count(genero_analisis, name = "celdas")
      print(as.data.frame(g), row.names = FALSE)
      write_csv(d, file.path(dir_mod, sprintf("celdas_saturadas_%dkm.csv", tam)))
    }
  }
}
