# =============================================================================
# 33_metodos_numeros.R
# Saca los números que faltan en Métodos y reconcilia la descarga de GBIF con
# la Tabla S1. El script SOLO LEE: no modifica ni escribe nada salvo un csv
# de resumen al final.
#
# Resuelve cuatro cosas:
#   (1) Qué contiene realmente la descarga dl.8mtewz (géneros, basisOfRecord,
#       conjuntos de datos) y de dónde pueden salir los 343,553 de la Tabla S1.
#   (2) La n de los modelos de SES-MPD y de riqueza, por tamaño y género.
#   (3) Las celdas que dependen de registros rescatados fuera del polígono.
#   (4) La extensión del área de estudio, para declararla en Métodos.
#
# Requiere: dplyr, readr, tidyr. Para (4) conviene sf, pero no es obligatorio.
# =============================================================================

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

base    <- RUTA_BASE
dir_cur <- file.path(base, "datos", "curacion")
dir_rej <- file.path(base, "datos", "rejilla")
dir_mod <- file.path(base, "datos", "modelos")
dir_out <- file.path(base, "tablas")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

tamanos_km <- c(50, 100, 200)
generos    <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")

# <<< AJUSTA ESTO >>> ruta del zip (o del .csv ya descomprimido) de dl.8mtewz
ruta_descarga <- file.path(base, "datos", "gbif", "0012345-260101000000000.zip")

linea <- function(txt) cat("\n", strrep("=", 70), "\n", txt, "\n",
                           strrep("=", 70), "\n", sep = "")

# =============================================================================
# 1. QUÉ HAY EN LA DESCARGA DE GBIF
# =============================================================================
linea("1. CONTENIDO DE LA DESCARGA (dl.8mtewz)")

leer_descarga <- function(ruta) {
  if (!file.exists(ruta)) {
    message("No encontré ", ruta,
            "\n  -> edita 'ruta_descarga' arriba y vuelve a correr el bloque 1.")
    return(NULL)
  }
  cols <- cols_only(
    gbifID                       = col_character(),
    datasetKey                   = col_character(),
    family                       = col_character(),
    genus                        = col_character(),
    species                      = col_character(),
    basisOfRecord                = col_character(),
    decimalLatitude              = col_double(),
    decimalLongitude             = col_double(),
    coordinateUncertaintyInMeters = col_double(),
    year                         = col_integer()
  )
  if (grepl("\\.zip$", ruta, ignore.case = TRUE)) {
    dentro <- utils::unzip(ruta, list = TRUE)$Name
    cat("Archivos dentro del zip:\n"); print(dentro)
    # formato "Simple": un solo .csv separado por tabuladores
    cand <- dentro[grepl("\\.(csv|txt)$", dentro, ignore.case = TRUE)]
    cand <- cand[!grepl("^dataset/|meta\\.xml|citations|rights", cand)]
    if (!length(cand)) { message("No vi un csv/txt dentro del zip."); return(NULL) }
    cat("Leyendo: ", cand[1], "\n", sep = "")
    con <- unz(ruta, cand[1])
    d <- read_tsv(con, col_types = cols, quote = "", progress = FALSE)
  } else {
    d <- read_tsv(ruta, col_types = cols, quote = "", progress = FALSE)
  }
  d
}

gb <- leer_descarga(ruta_descarga)

if (!is.null(gb)) {
  cat("\nRegistros en el archivo: ", nrow(gb),
      "   (la página dice 148,424)\n", sep = "")
  cat("Conjuntos de datos distintos: ", n_distinct(gb$datasetKey),
      "   (la página dice 194)\n", sep = "")
  cat("Familias: ", paste(sort(unique(na.omit(gb$family))), collapse = ", "), "\n", sep = "")

  cat("\n--- Registros por género ---\n")
  print(as.data.frame(count(gb, genus, name = "registros") |> arrange(desc(registros))),
        row.names = FALSE)

  cat("\n--- Registros por basisOfRecord ---\n")
  print(as.data.frame(count(gb, basisOfRecord, name = "registros")), row.names = FALSE)

  cat("\n--- Géneros x basisOfRecord ---\n")
  print(as.data.frame(count(gb, genus, basisOfRecord, name = "n") |>
                        pivot_wider(names_from = basisOfRecord, values_from = n,
                                    values_fill = 0)), row.names = FALSE)

  cat("\n--- Especies distintas por género ---\n")
  print(as.data.frame(gb |> group_by(genus) |>
                        summarise(especies = n_distinct(species), .groups = "drop")),
        row.names = FALSE)

  cat("\n--- Años: rango observado ---\n")
  cat("min: ", min(gb$year, na.rm = TRUE), "   max: ", max(gb$year, na.rm = TRUE),
      "   sin año: ", sum(is.na(gb$year)), "\n", sep = "")

  # Las dos preguntas que importan para la Tabla S1
  cat("\n>>> ¿Está Anoura?      ", "Anoura"     %in% gb$genus, "\n", sep = "")
  cat(">>> ¿Está Chiroderma?  ", "Chiroderma" %in% gb$genus, "\n", sep = "")
  cat(">>> Registros de los 5 géneros focales: ",
      sum(gb$genus %in% c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira", "Anoura"),
          na.rm = TRUE), "\n", sep = "")
  cat(">>> Solo 'PRESERVED_SPECIMEN' de los 5 focales: ",
      sum(gb$genus %in% c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira", "Anoura") &
          gb$basisOfRecord == "PRESERVED_SPECIMEN", na.rm = TRUE),
      "   (la Tabla S1 dice 142,926)\n", sep = "")
}

# -----------------------------------------------------------------------------
# 1b. Primer paso de la cadena de curación, para comparar
# -----------------------------------------------------------------------------
linea("1b. PRIMER PASO DE LA CURACIÓN (especimenes_con_marcas.rds)")

ruta_marcas <- file.path(dir_cur, "especimenes_con_marcas.rds")
if (file.exists(ruta_marcas)) {
  m <- readRDS(ruta_marcas)
  cat("filas: ", nrow(m), "   columnas: ", ncol(m), "\n", sep = "")
  cat("nombres: ", paste(names(m), collapse = ", "), "\n", sep = "")
  cg <- intersect(c("genero_analisis", "genero", "genus"), names(m))
  if (length(cg)) {
    cat("\n--- Registros por género en este paso ---\n")
    print(as.data.frame(count(m, .data[[cg[1]]], name = "registros")), row.names = FALSE)
  }
  cb <- intersect(c("basisOfRecord", "basis_of_record"), names(m))
  if (length(cb)) {
    cat("\n--- basisOfRecord en este paso ---\n")
    print(as.data.frame(count(m, .data[[cb[1]]], name = "registros")), row.names = FALSE)
  }
} else {
  message("No encontré ", ruta_marcas)
}

# =============================================================================
# 2. n DE LOS MODELOS (riqueza y SES-MPD)
# =============================================================================
linea("2. n DE LOS MODELOS, POR TAMAÑO Y GÉNERO")

hip          <- c("topo_sd_elev", "hist_sd_temp")
ctrl_riq     <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2",
                  "clima_bio01_media", "clima_bio12_media")
ctrl_ses     <- c("ctrl_esfuerzo_log", "clima_bio01_media",
                  "clima_bio12_media", "log_S")
n_min_modelo <- 25

n_modelos <- lapply(tamanos_km, function(tam) {
  ruta <- file.path(dir_mod, paste0("base_modelos_", tam, "km.csv"))
  if (!file.exists(ruta)) { message("No encontré ", ruta); return(NULL) }
  b <- read_csv(ruta, show_col_types = FALSE) |> filter(genero %in% generos)

  vr <- c(hip, ctrl_riq)
  vs <- c(hip, ctrl_ses)

  riq <- b |> filter(if_all(all_of(vr), ~ !is.na(.x))) |>
    count(genero, name = "n_riqueza")

  ses <- b |> filter(filo_valida, !is.na(ses_mpd),
                     if_all(all_of(vs), ~ !is.na(.x))) |>
    count(genero, name = "n_ses_mpd")

  full_join(riq, ses, by = "genero") |>
    mutate(tam_km = tam,
           entra_riqueza = !is.na(n_riqueza) & n_riqueza >= n_min_modelo,
           entra_ses     = !is.na(n_ses_mpd) & n_ses_mpd >= n_min_modelo)
}) |> bind_rows()

if (nrow(n_modelos)) {
  print(as.data.frame(n_modelos), row.names = FALSE)
  cat("\n--- Totales por tamaño (la n que va en Métodos) ---\n")
  print(as.data.frame(
    n_modelos |> group_by(tam_km) |>
      summarise(celdas_riqueza = sum(n_riqueza[entra_riqueza], na.rm = TRUE),
                celdas_ses     = sum(n_ses_mpd[entra_ses],     na.rm = TRUE),
                generos_riqueza = sum(entra_riqueza),
                generos_ses     = sum(entra_ses), .groups = "drop")),
    row.names = FALSE)
  write_csv(n_modelos, file.path(dir_out, "n_modelos.csv"))
}

# =============================================================================
# 3. CELDAS QUE DEPENDEN DE REGISTROS RESCATADOS
# =============================================================================
linea("3. CELDAS RESCATADAS POR EL MARGEN DE 100 km")

esp_ruta <- file.path(dir_rej, "especies_celdas_ocupadas.csv")
esp <- if (file.exists(esp_ruta)) read_csv(esp_ruta, show_col_types = FALSE) else NULL

rescatadas <- lapply(tamanos_km, function(tam) {
  r_std <- file.path(dir_rej, paste0("presencia_", tam, "km.csv"))
  r_est <- file.path(dir_rej, paste0("presencia_", tam, "km_rango_estricto.csv"))
  if (!file.exists(r_std) || !file.exists(r_est)) {
    message("Falta presencia estándar o estricta de ", tam, " km; lo salto.")
    return(NULL)
  }
  a <- read_csv(r_std, show_col_types = FALSE) |> distinct(genero_analisis, celda, especie)
  b <- read_csv(r_est, show_col_types = FALSE) |> distinct(genero_analisis, celda, especie)

  ga <- a |> distinct(genero_analisis, celda)
  gb2 <- b |> distinct(genero_analisis, celda)

  tibble(
    tam_km                   = tam,
    celdas_std               = n_distinct(a$celda),
    celdas_est               = n_distinct(b$celda),
    celdas_perdidas          = n_distinct(a$celda) - n_distinct(b$celda),
    pares_genero_celda_std   = nrow(ga),
    pares_genero_celda_est   = nrow(gb2),
    pares_perdidos           = nrow(anti_join(ga, gb2, by = c("genero_analisis", "celda"))),
    registros_especie_perdidos = nrow(anti_join(a, b,
                                   by = c("genero_analisis", "celda", "especie")))
  )
}) |> bind_rows()

if (nrow(rescatadas)) {
  print(as.data.frame(rescatadas), row.names = FALSE)
  cat("\nCompara 'pares_perdidos' y 'celdas_perdidas' con los 170 / 114 / 76\n",
      "del manuscrito y dime cuál de los dos coincide.\n", sep = "")
  write_csv(rescatadas, file.path(dir_out, "celdas_rescatadas.csv"))
}

# =============================================================================
# 4. EXTENSIÓN DEL ÁREA DE ESTUDIO
# =============================================================================
linea("4. EXTENSIÓN DEL ÁREA DE ESTUDIO")

crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
hay_sf   <- requireNamespace("sf", quietly = TRUE)

extension <- lapply(tamanos_km, function(tam) {
  ruta <- file.path(dir_mod, paste0("base_modelos_", tam, "km.csv"))
  if (!file.exists(ruta)) return(NULL)
  b <- read_csv(ruta, show_col_types = FALSE) |> filter(genero %in% generos)
  if (!all(c("x", "y") %in% names(b))) {
    message("No hay columnas x / y en ", basename(ruta)); return(NULL)
  }
  u <- b |> distinct(celda, x, y)
  area_col <- intersect("ctrl_area_tierra_km2", names(b))

  out <- tibble(
    tam_km      = tam,
    celdas      = nrow(u),
    x_min_km    = round(min(u$x) / 1000), x_max_km = round(max(u$x) / 1000),
    y_min_km    = round(min(u$y) / 1000), y_max_km = round(max(u$y) / 1000),
    area_tierra_km2 = if (length(area_col))
      round(sum(distinct(b, celda, .keep_all = TRUE)[[area_col]], na.rm = TRUE)) else NA_real_
  )

  if (hay_sf) {
    esq <- data.frame(x = c(min(u$x), max(u$x), min(u$x), max(u$x)),
                      y = c(min(u$y), min(u$y), max(u$y), max(u$y)))
    ll  <- sf::st_coordinates(
             sf::st_transform(sf::st_as_sf(esq, coords = c("x", "y"), crs = crs_laea),
                              4326))
    out$lon_min <- round(min(ll[, 1]), 2); out$lon_max <- round(max(ll[, 1]), 2)
    out$lat_min <- round(min(ll[, 2]), 2); out$lat_max <- round(max(ll[, 2]), 2)
  }
  out
}) |> bind_rows()

if (nrow(extension)) {
  print(as.data.frame(extension), row.names = FALSE)
  if (!hay_sf) cat("\n(Instala sf si quieres los límites en grados.)\n")
  write_csv(extension, file.path(dir_out, "extension_estudio.csv"))
}

linea("PENDIENTE QUE NO SALE DE AQUÍ")
cat("La fecha de descarga de la capa de esfuerzo (Mammalia, dl.2u7xpj).\n",
    "Ábrela en GBIF y mándame la fecha que diga 'Download'.\n", sep = "")

cat("\nListo. Salidas en ", dir_out, ":\n",
    "  n_modelos.csv\n  celdas_rescatadas.csv\n  extension_estudio.csv\n", sep = "")
