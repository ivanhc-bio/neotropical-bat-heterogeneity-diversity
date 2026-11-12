# 03_curacion_coordenadas.R
# Curación de coordenadas de los ESPECÍMENES ya armonizados con MDD v2.5
# (salida de 02_taxonomia_referencia.R).
#
# Qué hace:
#   1. Filtra por incertidumbre de la coordenada.
#   2. Aplica pruebas de CoordinateCleaner (ceros, lon = lat, centroides, capitales, instituciones).
#   3. Compara la elevación del registro con un DEM y con rangos altitudinales por especie
#      (los marca para revisión; no los elimina salvo que lo pidas).
#   4. Guarda los registros con sus marcas, los curados y la bitácora.
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.
# Instalar una vez: install.packages(c("CoordinateCleaner", "terra"))

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(CoordinateCleaner)
library(terra)

# terra también exporta select() y, al cargarse después, tapa a la de dplyr
# ("el argumento no fue usado"). Se fija la de dplyr para todo el script.
select <- dplyr::select

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir_tax  <- file.path(dir_base, "datos", "taxonomia")
dir_amb  <- file.path(dir_base, "datos", "ambiente")
dir_cur  <- file.path(dir_base, "datos", "curacion")
dir.create(dir_amb, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_cur, recursive = TRUE, showWarnings = FALSE)

# PARÁMETROS (valores propuestos por mí, no de la literatura: confírmalos o ajústalos)
umbral_incertidumbre_m <- 10000  # se eliminan registros con incertidumbre mayor; los NA se conservan
umbral_dif_elev_m      <- 500    # diferencia (m) entre elevación del registro y DEM a partir de la cual se MARCA
umbral_dif_sens_m      <- 1000   # diferencia (m) a partir de la cual se marca para la sensibilidad (no elimina)
umbral_dif_elim_m      <- 2000   # diferencia (m) a partir de la cual se ELIMINA (coordenada probablemente errónea)
tolerancia_alt_m       <- 200    # tolerancia sobre los rangos altitudinales por especie
eliminar_revision      <- FALSE  # TRUE: elimina también los marcados por elevación (hazlo tras revisarlos)

# DEM: cualquier ráster de elevación global o neotropical en metros (idealmente el mismo que
# usarás para la SD de elevación). Si no existe el archivo, se omiten las pruebas de elevación.
dem_ruta <- file.path(dir_amb, "dem.tif")
hay_dem  <- file.exists(dem_ruta)

# ---- 1. Entrada ---------------------------------------------------------------
esp <- readRDS(file.path(dir_tax, "especimenes_taxonomia_armonizada.rds"))
n_entrada <- nrow(esp)

requeridas <- c("decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters",
                "elevation", "nombre_final", "genero_analisis")
faltan <- setdiff(requeridas, names(esp))
if (length(faltan) > 0) {
  stop("Faltan columnas: ", paste(faltan, collapse = ", "),
       ". Vuelve a correr 02_taxonomia_referencia.R (versión con 'genero_analisis').")
}

esp <- esp %>%
  mutate(decimalLongitude = as.numeric(decimalLongitude),
         decimalLatitude  = as.numeric(decimalLatitude)) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude),
         between(decimalLongitude, -180, 180), between(decimalLatitude, -90, 90))

cat("\nEspecímenes de entrada:", n_entrada, "| con coordenadas válidas:", nrow(esp), "\n")
n_validas <- nrow(esp)

# ---- 2. Incertidumbre de la coordenada ---------------------------------------
cat("\nIncertidumbre (m): registros sin dato =", sum(is.na(esp$coordinateUncertaintyInMeters)),
    "| mayores que el umbral =",
    sum(esp$coordinateUncertaintyInMeters > umbral_incertidumbre_m, na.rm = TRUE), "\n")
print(summary(esp$coordinateUncertaintyInMeters))

esp <- esp %>%
  mutate(flag_incertidumbre = !is.na(coordinateUncertaintyInMeters) &
           coordinateUncertaintyInMeters > umbral_incertidumbre_m)

# ---- 3. Pruebas de CoordinateCleaner -----------------------------------------
# Con value = "flagged" cada función devuelve TRUE para los registros que PASAN la prueba;
# por eso se niega para guardar la marca (TRUE = marcado).
x <- as.data.frame(esp)
prueba <- function(f, ...) {
  !f(x, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged", verbose = FALSE, ...)
}

esp$flag_cc_zero <- prueba(cc_zero)                    # en (0, 0) o a menos de 0,5 grados
esp$flag_cc_equ  <- prueba(cc_equ, test = "absolute")  # |lon| == |lat|
esp$flag_cc_cen  <- prueba(cc_cen)                     # centroides de países y provincias (1 km)
esp$flag_cc_cap  <- prueba(cc_cap)                     # capitales (10 km)
esp$flag_cc_inst <- prueba(cc_inst)                    # instituciones de biodiversidad (100 m)
# DECISIÓN (4-nov-2026): flag_cc_cap SOLO MARCA, no elimina. De 2 619 registros marcados,
# 764 son de Panamá y el resto se concentra en países pequeños e islas (TT, GD, VC, CR), donde
# 10 km alrededor de la capital cubren buena parte del territorio colectable; la prueba
# descartaría localidades legítimas. La marca se conserva para poder repetir el análisis sin ellos.

# ---- 4. Elevación: DEM y rangos altitudinales ---------------------------------
esp$elev_dem <- NA_real_
if (hay_dem) {
  dem <- rast(dem_ruta)
  pts <- vect(as.data.frame(esp[, c("decimalLongitude", "decimalLatitude")]),
              geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  pts <- project(pts, crs(dem))
  esp$elev_dem <- terra::extract(dem, pts)[, 2]
} else {
  message("No hay DEM en ", dem_ruta, ": se omiten las pruebas de elevación. ",
          "Ponlo ahí (GeoTIFF en metros) y vuelve a correr desde esta sección.")
}

# elevation == 0 se trata como dato faltante (en estos especímenes casi siempre es un
# valor por omisión, no el nivel del mar): no entra en la comparación con el DEM.
esp <- esp %>%
  mutate(
    elev_cmp = if_else(!is.na(elevation) & elevation == 0, NA_real_, as.numeric(elevation)),
    dif_elev_abs = abs(elev_cmp - elev_dem),
    flag_dem_sin_dato = hay_dem & is.na(elev_dem),
    flag_dif_elev        = !is.na(dif_elev_abs) & dif_elev_abs > umbral_dif_elev_m,
    flag_dif_elev_sens   = !is.na(dif_elev_abs) & dif_elev_abs > umbral_dif_sens_m &
                           dif_elev_abs <= umbral_dif_elim_m,
    flag_dif_elev_extrema = !is.na(dif_elev_abs) & dif_elev_abs > umbral_dif_elim_m
  )
cat("\nElevación = 0 tratada como faltante en", sum(!is.na(esp$elevation) & esp$elevation == 0),
    "registros\n")

# Rangos altitudinales por especie (de la literatura): plantilla para completar.
archivo_alt <- file.path(dir_cur, "rangos_altitudinales.csv")
if (!file.exists(archivo_alt)) {
  plantilla_alt <- esp %>%
    distinct(genero_analisis, nombre_final) %>%
    arrange(genero_analisis, nombre_final) %>%
    transmute(genero_analisis, especie = nombre_final,
              elev_min = NA_real_, elev_max = NA_real_, fuente = NA_character_)
  write_csv(plantilla_alt, archivo_alt)
  message("Se creó la plantilla ", archivo_alt, ". Complétala (m s. n. m., con la fuente) ",
          "y vuelve a correr desde esta sección.")
  esp$flag_fuera_rango_alt <- FALSE
} else {
  alt <- read_csv(archivo_alt, show_col_types = FALSE) %>%
    distinct(especie, .keep_all = TRUE) %>%
    select(especie, elev_min, elev_max)
  # Se quitan las columnas de una corrida anterior de esta sección; si no, la unión las
  # duplica (elev_min.x / elev_min.y) y mutate() no encuentra elev_min.
  esp <- esp %>%
    select(-any_of(c("elev_min", "elev_max", "elev_ref"))) %>%
    left_join(alt, by = c("nombre_final" = "especie")) %>%
    mutate(
      elev_ref = coalesce(elev_dem, elevation),
      flag_fuera_rango_alt = !is.na(elev_min) & !is.na(elev_max) & !is.na(elev_ref) &
        (elev_ref < elev_min - tolerancia_alt_m | elev_ref > elev_max + tolerancia_alt_m)
    )
  cat("\nEspecies con rango altitudinal completo:", sum(!is.na(alt$elev_min) & !is.na(alt$elev_max)),
      "de", nrow(alt), "\n")
}

# ---- 5. Consolidación ---------------------------------------------------------
flags_auto <- c("flag_incertidumbre", "flag_cc_zero", "flag_cc_equ", "flag_cc_cen",
                "flag_cc_inst", "flag_dif_elev_extrema")  # esta última: |dif| > umbral_dif_elim_m
flags_rev  <- c("flag_dif_elev", "flag_fuera_rango_alt")
flags_info <- "flag_cc_cap"   # solo informativa: se conserva el registro
flags_sens <- "flag_dif_elev_sens"  # |dif| entre 1 000 y 2 000 m: se conserva; para repetir modelos sin ellos

esp <- esp %>%
  mutate(
    quitar_auto = if_any(all_of(flags_auto)),
    # "revisar" excluye lo ya eliminado automáticamente
    revisar     = if_any(all_of(flags_rev)) & !quitar_auto,
    conservar   = !quitar_auto & !(eliminar_revision & revisar)
  )

pruebas_resumen <- c(flags_auto, flags_info, flags_rev, flags_sens, "flag_dem_sin_dato")
resumen <- tibble(
  prueba = pruebas_resumen,
  n_marcados = vapply(pruebas_resumen, function(f) sum(esp[[f]]), numeric(1)),
  accion = c(rep("elimina", length(flags_auto)), rep("solo marca", length(flags_info)),
             rep("revisar", length(flags_rev)), rep("sensibilidad", length(flags_sens)), "revisar")
)
cat("\nRegistros marcados por prueba (una misma fila puede estar en varias):\n")
print(resumen)

# Pérdida acumulada en el orden de las pruebas (para la bitácora)
acum <- rep(TRUE, nrow(esp))
acumulado <- tibble(paso = character(), n_registros = numeric())
for (f in flags_auto) {
  acum <- acum & !esp[[f]]
  acumulado <- bind_rows(acumulado, tibble(paso = paste("Tras", f), n_registros = sum(acum)))
}
cat("\nRegistros que quedan tras cada prueba automática:\n")
print(acumulado)

cat("\nPor género de análisis:\n")
print(esp %>%
        group_by(genero_analisis) %>%
        summarise(entrada = n(), conservados = sum(conservar),
                  porcentaje = round(100 * conservados / entrada, 1),
                  para_revisar = sum(revisar), .groups = "drop"))

# ---- 6. Salidas ---------------------------------------------------------------
cols_rev <- intersect(c("gbifID", "nombre_final", "genero_analisis", "countryCode", "stateProvince",
                        "locality", "decimalLatitude", "decimalLongitude",
                        "coordinateUncertaintyInMeters", "elevation", "elev_dem",
                        "elev_min", "elev_max", "institutionCode", "catalogNumber",
                        "dif_elev_abs", flags_rev, flags_sens, "flag_dem_sin_dato"), names(esp))
write_csv(esp %>% filter(revisar | flag_dem_sin_dato) %>% select(all_of(cols_rev)),
          file.path(dir_cur, "registros_para_revision_altitud.csv"))
write_csv(resumen, file.path(dir_cur, "resumen_marcas_coordenadas.csv"))

saveRDS(esp, file.path(dir_cur, "especimenes_con_marcas.rds"))
saveRDS(esp %>% filter(conservar), file.path(dir_cur, "especimenes_curados_coord.rds"))

bitacora_archivo <- file.path(dir_gbif, "bitacora_curacion.csv")
nueva <- tibble(
  paso = c("Coordenadas válidas", "Curación de coordenadas (filtros automáticos)",
           "Marcados para revisión de elevación (no eliminados)"),
  n_registros = c(n_validas, sum(esp$conservar), sum(esp$revisar))
)
bit <- if (file.exists(bitacora_archivo)) {
  bind_rows(read_csv(bitacora_archivo, show_col_types = FALSE) %>% filter(!paso %in% nueva$paso), nueva)
} else nueva
write_csv(bit, bitacora_archivo)
print(bit)
