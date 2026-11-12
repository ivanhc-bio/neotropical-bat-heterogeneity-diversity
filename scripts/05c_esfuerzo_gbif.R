# 05c_esfuerzo_gbif.R
# Esfuerzo de muestreo INDEPENDIENTE de los géneros focales, por celda LAEA (50, 100 y 200 km).
#
# Por qué: el esfuerzo de 05 (ctrl_esfuerzo_log) se calcula con los registros de los propios géneros focales, así que es
# circular con la riqueza (más registros = más especies por pura mecánica). Aquí se usa, como en el manuscrito anterior,
# el número de registros de GBIF de otros grupos que reflejan cuánto se ha muestreado cada celda:
#   esf_mam   : todos los registros de Mammalia (como en el borrador, 1960-2025)
#   esf_quir  : registros de Chiroptera SIN los cinco géneros focales (Anoura, Artibeus, Dermanura, Carollia,
#               Platyrrhinus, Sturnira): refleja mejor el muestreo de murciélagos
#   esf_quir_ej: igual que esf_quir pero solo ejemplares de museo (PRESERVED_SPECIMEN), que es el tipo de dato de los
#               registros focales
# Para cada una se guarda el conteo, su log10(n + 1) y el número de localidades distintas (coordenadas redondeadas a 0.01°).
#
# Pasos (el script detecta qué falta y continúa donde se quedó):
#   1) solicita a GBIF la descarga de Mammalia (requiere cuenta de GBIF; ver 'Credenciales');
#   2) espera, descarga y descomprime;
#   3) cuenta los registros por celda y escribe datos/predictores/esfuerzo_independiente_{tam}km.csv.
#
# Credenciales: en R, usethis::edit_r_environ() y añade
#   GBIF_USER="tu_usuario"  GBIF_PWD="tu_contraseña"  GBIF_EMAIL="tu_correo"
# y reinicia R.
#
# No probado: no tengo R ni acceso a GBIF en mi entorno. La descarga de Mammalia para esta región puede tardar
# (de minutos a un par de horas) y pesar varios GB. Guarda el DOI que imprime para citarlo.
# Requiere: rgbif, data.table, sf, dplyr, readr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(rgbif)
library(data.table)
library(sf)
library(dplyr)
library(readr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif_esfuerzo")
dir_pred <- file.path(dir_base, "datos", "predictores")
dir.create(dir_gbif, showWarnings = FALSE, recursive = TRUE)

crs_laea   <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"   # igual que 04
tamanos_km <- c(50, 100, 200)
anio_min <- 1960; anio_max <- 2025
region_wkt <- "POLYGON((-120 -40, -30 -40, -30 35, -120 35, -120 -40))"   # lon_min, lat_min, lon_max, lat_max
incertidumbre_max_m <- 10000          # se descartan registros con incertidumbre mayor; NA se conserva
generos_focales <- c("Anoura", "Artibeus", "Dermanura", "Carollia", "Platyrrhinus", "Sturnira")
archivo_clave <- file.path(dir_gbif, "descarga_mammalia_clave.txt")

# ---- 1. Solicitar y descargar -------------------------------------------------------------
if (!file.exists(archivo_clave)) {
  cat("Solicitando la descarga de Mammalia a GBIF...\n")
  sol <- occ_download(
    pred("taxonKey", 359),                          # 359 = Mammalia
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_gte("year", anio_min), pred_lte("year", anio_max),
    pred_within(region_wkt),
    pred_not(pred("basisOfRecord", "FOSSIL_SPECIMEN")),
    format = "SIMPLE_CSV")
  writeLines(as.character(sol), archivo_clave)
  cat("Clave de la descarga:", as.character(sol), "\n")
}
clave <- readLines(archivo_clave, n = 1)
zip_local <- file.path(dir_gbif, paste0(clave, ".zip"))
if (!file.exists(zip_local)) {
  occ_download_wait(clave)
  occ_download_get(clave, path = dir_gbif, overwrite = FALSE)
}
meta <- occ_download_meta(clave)
cat("DOI de la descarga (cítalo en métodos):", meta$doi, "| registros:", meta$totalRecords, "\n")

csv_local <- file.path(dir_gbif, paste0(clave, ".csv"))
if (!file.exists(csv_local)) unzip(zip_local, files = paste0(clave, ".csv"), exdir = dir_gbif)

# ---- 2. Leer solo las columnas necesarias ----------------------------------------------------
cols <- c("order", "genus", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "basisOfRecord")
d <- fread(csv_local, select = cols, quote = "", showProgress = TRUE)
cat("Registros leídos:", nrow(d), "\n")
d <- d[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
d <- d[is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= incertidumbre_max_m]
cat("Registros tras el filtro de incertidumbre:", nrow(d), "\n")

d[, quir := order == "Chiroptera" & !(genus %in% generos_focales)]
d[, ejemplar := basisOfRecord == "PRESERVED_SPECIMEN"]
d[, loc := paste(round(decimalLatitude, 2), round(decimalLongitude, 2))]

# ---- 3. Asignar a celdas LAEA y contar --------------------------------------------------
xy <- sf_project(from = "EPSG:4326", to = crs_laea, pts = as.matrix(d[, .(decimalLongitude, decimalLatitude)]))
d[, `:=`(x = xy[, 1], y = xy[, 2])]
d <- d[is.finite(x) & is.finite(y)]

contar <- function(sub, prefijo) {
  r <- sub[, .(n = .N, nloc = uniqueN(loc)), by = celda]
  setnames(r, c("n", "nloc"), c(paste0(prefijo, "_n"), paste0(prefijo, "_loc")))
  r
}
for (tam_km in tamanos_km) {
  tam <- tam_km * 1000
  d[, celda := sprintf("c%d_%d_%d", as.integer(tam_km), as.integer(floor(x / tam)), as.integer(floor(y / tam)))]
  e <- Reduce(function(a, b) merge(a, b, by = "celda", all = TRUE),
              list(contar(d, "esf_mam"), contar(d[quir == TRUE], "esf_quir"), contar(d[quir == TRUE & ejemplar == TRUE], "esf_quir_ej")))
  e <- as_tibble(e) %>% mutate(across(-celda, ~ ifelse(is.na(.x), 0L, .x)))
  for (p in c("esf_mam", "esf_quir", "esf_quir_ej")) e[[paste0(p, "_log")]] <- log10(e[[paste0(p, "_n")]] + 1)
  write_csv(e, file.path(dir_pred, sprintf("esfuerzo_independiente_%dkm.csv", tam_km)))
  cat(sprintf("%d km: %d celdas con registros | mediana de registros de Mammalia: %.0f, Chiroptera sin focales: %.0f, ejemplares: %.0f\n",
              tam_km, nrow(e), median(e$esf_mam_n), median(e$esf_quir_n), median(e$esf_quir_ej_n)))
}

# ---- 4. Comparación con el esfuerzo focal de 05 (solo diagnóstico) ----------------------------
f100 <- file.path(dir_pred, "predictores_100km.csv")
if (file.exists(f100)) {
  p <- read_csv(f100, show_col_types = FALSE) %>% select(celda, ctrl_esfuerzo_log)
  e <- read_csv(file.path(dir_pred, "esfuerzo_independiente_100km.csv"), show_col_types = FALSE) %>%
    right_join(p, by = "celda") %>% mutate(across(starts_with("esf_"), ~ ifelse(is.na(.x), 0, .x)))
  cat("\nCorrelación de Spearman entre el esfuerzo focal (05) y los independientes (100 km, celdas de la rejilla):\n")
  print(round(cor(e[, c("ctrl_esfuerzo_log", "esf_mam_log", "esf_quir_log", "esf_quir_ej_log")], method = "spearman"), 2))
}
cat("\nListo. Siguiente paso: en 09_bases_modelos.R, esfuerzo_independiente <- 'esf_quir_ej_log' (o la que elijas).\n")
