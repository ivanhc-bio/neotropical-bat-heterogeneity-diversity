# 01_descarga_gbif.R
# Descarga de registros de Phyllostomidae desde GBIF en UNA sola solicitud.
# Sirve para (a) los cinco géneros focales y (b) el esfuerzo de muestreo (todos los filostómidos).
# Sin ejecutar todavía: revisa los parámetros marcados con "AJUSTAR" antes de correrlo.

# ---- 0. Paquetes y credenciales ----------------------------------------------
# install.packages(c("rgbif", "dplyr", "readr"))
# Credenciales de tu cuenta de GBIF: agrégalas a ~/.Renviron y reinicia R.
#   usethis::edit_r_environ()   # abre el archivo
#   GBIF_USER="tu_usuario"
#   GBIF_PWD="tu_contraseña"
#   GBIF_EMAIL="tu_correo"
# rgbif las lee solas; no las escribas dentro del script.
source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(rgbif)
library(dplyr)
library(readr)

stopifnot(
  nzchar(Sys.getenv("GBIF_USER")),
  nzchar(Sys.getenv("GBIF_PWD")),
  nzchar(Sys.getenv("GBIF_EMAIL"))
)

# Carpeta del proyecto (en Windows usa "/" o "\\", no "\")
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir.create(dir_gbif, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Resolver los taxones en el backbone de GBIF ---------------------------
familia <- name_backbone(name = "Phyllostomidae", rank = "family", kingdom = "Animalia")

generos_focales <- c("Sturnira", "Platyrrhinus", "Artibeus", "Dermanura", "Carollia", "Anoura")

campo <- function(x) if (is.null(x)) NA_character_ else as.character(x)

claves <- bind_rows(lapply(generos_focales, function(g) {
  # family = "Phyllostomidae" ayuda a desambiguar homónimos con géneros de otros animales
  r <- name_backbone(name = g, rank = "genus", family = "Phyllostomidae", kingdom = "Animalia")
  tibble(
    genero    = g,
    usageKey  = campo(r$usageKey),
    status    = campo(r$status),
    matchType = campo(r$matchType),
    rango     = campo(r$rank),
    familia   = campo(r$family)
  )
}))
print(claves)  # Revisa: familia == "Phyllostomidae" y matchType == "EXACT" en todos.
# Dermanura puede aparecer como sinónimo de Artibeus según la versión del backbone;
# la descarga es por familia, así que ambos nombres quedan cubiertos.

stopifnot(!is.null(familia$usageKey), identical(familia$matchType, "EXACT"))

problemas <- claves %>%
  filter(is.na(familia) | familia != "Phyllostomidae" | matchType != "EXACT")
if (nrow(problemas) > 0) {
  print(problemas)
  stop(
    "Género(s) sin coincidencia exacta dentro de Phyllostomidae. ",
    "Prueba name_backbone(name = \"<género>\", rank = \"genus\", family = \"Phyllostomidae\", verbose = TRUE) ",
    "para ver las alternativas, o revisa la grafía (p. ej. Platyrrhinus vs Platyrhinus)."
  )
}

# ---- 2. Conteo previo por género (con los mismos filtros de la descarga) -----
conteo_previo <- bind_rows(lapply(seq_len(nrow(claves)), function(i) {
  n <- occ_search(
    taxonKey = claves$usageKey[i],
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE,
    limit = 0
  )$meta$count
  tibble(genero = claves$genero[i], n_registros = n)
}))
print(conteo_previo)

# ---- 3. Solicitud de descarga --------------------------------------------------
# AJUSTAR: polígono amplio del Neotrópico (lon -120 a -30, lat -56 a 33), en sentido
# antihorario como exige GBIF. Se refina por país o elevación en el paso de curación.
wkt_neotropico <- "POLYGON((-120 -56, -30 -56, -30 33, -120 33, -120 -56))"

descarga <- occ_download(
  pred("taxonKey", familia$usageKey),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("occurrenceStatus", "PRESENT"),
  pred_in("basisOfRecord", c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION")),
  pred_within(wkt_neotropico),
  format = "SIMPLE_CSV"  # AJUSTAR: usa "DWCA" si necesitas más campos (p. ej. identifiedBy, typeStatus)
)

occ_download_wait(descarga)  # puede tardar de minutos a una hora

# ---- 4. Descargar, importar y registrar la cita -------------------------------
zip_path <- occ_download_get(descarga, path = dir_gbif, overwrite = TRUE)
datos <- occ_download_import(zip_path)

meta <- occ_download_meta(descarga)
cita <- gbif_citation(meta)$download
writeLines(c(paste("DOI:", meta$doi), cita), file.path(dir_gbif, "descarga_doi.txt"))

# ---- 5. Guardar datos y bitácora inicial --------------------------------------
saveRDS(datos, file.path(dir_gbif, "phyllostomidae_bruto.rds"))

focales <- datos %>% filter(genus %in% generos_focales)
saveRDS(focales, file.path(dir_gbif, "generos_focales_bruto.rds"))

bitacora <- tibble(
  paso = c("Descarga bruta (familia)", "Géneros focales"),
  n_registros = c(nrow(datos), nrow(focales))
)
write_csv(bitacora, file.path(dir_gbif, "bitacora_curacion.csv"))

focales %>%
  count(genus, basisOfRecord) %>%
  print()
