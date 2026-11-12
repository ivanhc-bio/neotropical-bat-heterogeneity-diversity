# 04b_filtro_rango.R
# Filtro de rango: distancia de cada registro al polígono de rango de su especie (nombre_auditado).
#
# Lógica (decidida el 5-nov-2026): se da un margen (buffer) alrededor de los polígonos porque los
# mapas de expertos tienen bordes imprecisos; lo que cae más allá se marca. Los polígonos también
# tienen huecos reales (p. ej. Anoura geoffroyi, Carollia castanea, Sturnira lilium en la Amazonia
# y los Andes), así que NO se excluye automáticamente todo lo que queda fuera: solo los casos
# aislados y muy lejanos. El resto queda en "revisar" y se decide por especie y país.
#
# Salidas (datos/rejilla y datos/curacion):
#   registros_fuera_rango.csv              registros más allá del buffer, con distancia
#   fuera_rango_por_especie_pais.csv       resumen para decidir (patrones taxonómicos vs. errores)
#   decisiones_rango_manual.csv            plantilla que editas: especie, país, decision, motivo
#   especimenes_curados_rango.rds          principal, sin los excluidos
#   especimenes_sensibilidad_rango.rds     sensibilidad taxonómica, sin los mismos excluidos
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(sf)

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir_cur  <- file.path(dir_base, "datos", "curacion")
dir_rang <- file.path(dir_base, "datos", "rangos")
dir_out  <- file.path(dir_base, "datos", "rejilla")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# PARÁMETROS (valores propuestos por mí: ajústalos tras ver el resumen)
buffer_km         <- 100   # margen alrededor del polígono: dentro de él, el registro se conserva
umbral_extremo_km <- 300   # más allá de esta distancia el registro se considera muy lejano
n_aislado_max     <- 3     # grupo especie-país con <= este número de registros fuera = caso aislado
crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# ---- 1. Entradas -------------------------------------------------------------
rec <- readRDS(file.path(dir_cur, "especimenes_curados_pais.rds"))
req <- c("gbifID", "nombre_auditado", "countryCode", "decimalLatitude", "decimalLongitude")
if (length(setdiff(req, names(rec))) > 0) stop("Faltan columnas: ", paste(setdiff(req, names(rec)), collapse = ", "))
rec <- rec %>% filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  mutate(especie = nombre_auditado)
n_entrada <- nrow(rec)
cat("\nRegistros de entrada:", n_entrada, "\n")

rangos <- st_read(file.path(dir_rang, "rangos_focales_v25.gpkg"), quiet = TRUE) %>%
  st_make_valid() %>% st_transform(crs_laea)

sin_poligono <- setdiff(unique(rec$especie), rangos$sciname)
if (length(sin_poligono) > 0) {
  cat("Especies SIN polígono (no se filtran; dibuja o asigna uno):", paste(sin_poligono, collapse = "; "), "\n")
}

# ---- 2. Distancia al polígono (por coordenada única y especie) -----------------
coords <- rec %>% distinct(especie, decimalLatitude, decimalLongitude)
pts <- st_as_sf(coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs_laea)
coords$dist_km <- NA_real_
for (sp in setdiff(unique(coords$especie), sin_poligono)) {
  i <- which(coords$especie == sp)
  poli <- st_geometry(rangos[rangos$sciname == sp, ])
  coords$dist_km[i] <- as.numeric(st_distance(pts[i, ], poli[1])) / 1000   # 0 si está dentro
}
rec <- rec %>% left_join(coords, by = c("especie", "decimalLatitude", "decimalLongitude"))

rec <- rec %>% mutate(
  estado_rango = case_when(
    is.na(dist_km)             ~ "sin_poligono",
    dist_km <= 0               ~ "dentro",
    dist_km <= buffer_km       ~ "buffer",
    dist_km <= umbral_extremo_km ~ "fuera_moderado",
    TRUE                       ~ "fuera_extremo"),
  fuera = estado_rango %in% c("fuera_moderado", "fuera_extremo"))
print(rec %>% count(estado_rango) %>% mutate(porcentaje = round(100 * n / sum(n), 2)))

# ---- 3. Resumen por especie y país --------------------------------------------
tot_esp <- rec %>% count(especie, name = "n_total_especie")
grupo <- rec %>% filter(fuera) %>%
  group_by(especie, countryCode) %>%
  summarise(n_fuera = n(), n_extremo = sum(estado_rango == "fuera_extremo"),
            dist_mediana_km = round(median(dist_km)), dist_max_km = round(max(dist_km)),
            .groups = "drop") %>%
  left_join(tot_esp, by = "especie") %>%
  mutate(pct_de_la_especie = round(100 * n_fuera / n_total_especie, 1),
         aislado = n_fuera <= n_aislado_max,
         decision_auto = if_else(aislado & n_extremo == n_fuera, "excluir_auto", "revisar")) %>%
  arrange(desc(n_fuera))
write_csv(grupo, file.path(dir_out, "fuera_rango_por_especie_pais.csv"))
cat("\nGrupos especie-país con registros fuera del buffer:", nrow(grupo),
    "| excluidos en automático:", sum(grupo$decision_auto == "excluir_auto"),
    "| a revisar:", sum(grupo$decision_auto == "revisar"), "\n")
print(head(grupo, 25))

# ---- 4. Decisiones manuales (especie + país) -----------------------------------
# Edita este archivo: decision = "excluir" o "conservar". Tiene prioridad sobre la regla automática.
arch_man <- file.path(dir_cur, "decisiones_rango_manual.csv")
if (!file.exists(arch_man)) {
  write_csv(tibble(nombre_auditado = character(), countryCode = character(),
                   decision = character(), motivo = character()), arch_man)
  cat("\nCreé la plantilla vacía:", arch_man, "\n")
}
man <- read_csv(arch_man, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  distinct(nombre_auditado, countryCode, .keep_all = TRUE) %>%
  select(especie = nombre_auditado, countryCode, decision_manual = decision)

rec <- rec %>%
  left_join(grupo %>% select(especie, countryCode, decision_auto), by = c("especie", "countryCode")) %>%
  left_join(man, by = c("especie", "countryCode")) %>%
  mutate(decision_rango = case_when(
    !fuera                                  ~ "conservar",
    !is.na(decision_manual)                 ~ decision_manual,
    decision_auto == "excluir_auto"         ~ "excluir",
    TRUE                                    ~ "revisar"))   # se conserva hasta decidir
print(rec %>% count(decision_rango))

write_csv(rec %>% filter(fuera) %>%
            select(any_of(c("gbifID", "especie", "countryCode", "stateProvince", "locality",
                            "decimalLatitude", "decimalLongitude", "institutionCode", "catalogNumber",
                            "dist_km", "estado_rango", "decision_rango"))) %>% arrange(desc(dist_km)),
          file.path(dir_out, "registros_fuera_rango.csv"))

# ---- 5. Guardar -----------------------------------------------------------------
excluidos <- rec %>% filter(decision_rango == "excluir") %>% pull(gbifID)
curados <- rec %>% filter(decision_rango != "excluir")
saveRDS(curados, file.path(dir_cur, "especimenes_curados_rango.rds"))
# Variante estricta para la sensibilidad: además quita todo lo que quedó más allá del buffer.
saveRDS(curados %>% filter(!fuera), file.path(dir_cur, "especimenes_curados_rango_estricto.rds"))
sens <- readRDS(file.path(dir_cur, "especimenes_sensibilidad_pais.rds")) %>% filter(!gbifID %in% excluidos)
saveRDS(sens, file.path(dir_cur, "especimenes_sensibilidad_rango.rds"))
cat("\nConservados tras el filtro de rango:", nrow(curados), "de", n_entrada,
    "| a revisar (conservados):", sum(curados$decision_rango == "revisar"), "\n")

bitacora_archivo <- file.path(dir_gbif, "bitacora_curacion.csv")
nueva <- tibble(
  paso = c(sprintf("Filtro de rango (excluye aislados a >%d km del polígono y decisiones manuales)", umbral_extremo_km),
           "Registros fuera del buffer pendientes de revisión (conservados)"),
  n_registros = c(nrow(curados), sum(curados$decision_rango == "revisar")))
bit <- if (file.exists(bitacora_archivo)) {
  bind_rows(read_csv(bitacora_archivo, show_col_types = FALSE) %>% filter(!paso %in% nueva$paso), nueva)
} else nueva
write_csv(bit, bitacora_archivo)
print(bit)

# A partir de aquí, el script 04 debe leer especimenes_curados_rango.rds (y _sensibilidad_rango.rds).
