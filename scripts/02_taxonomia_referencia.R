# 02_taxonomia_referencia.R
# Lista taxonómica de referencia (Mammal Diversity Database, MDD) para los cinco géneros focales
# y armonización de los nombres de los ESPECÍMENES de GBIF con esa lista.
# No probado todavía: corre por secciones y revisa lo que imprime cada una.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir_tax  <- file.path(dir_base, "datos", "taxonomia")
dir.create(dir_tax, recursive = TRUE, showWarnings = FALSE)

# Versión de MDD fijada para todo el proyecto (cámbiala solo de forma deliberada).
mdd_version <- "v2.5"
mdd_doi     <- "10.5281/zenodo.21654811"
mdd_archivo <- file.path(dir_tax, "MDD_v2.5_6904species.csv")
mdd_url     <- "https://zenodo.org/records/21654811/files/MDD_v2.5_6904species.csv?download=1"

generos_focales <- c("Sturnira", "Platyrrhinus", "Artibeus", "Dermanura", "Carollia", "Anoura")

# ---- 1. Descargar y leer MDD --------------------------------------------------
if (!file.exists(mdd_archivo)) {
  options(timeout = 600)
  ok <- tryCatch(
    download.file(mdd_url, mdd_archivo, mode = "wb"),
    error = function(e) e
  )
  if (inherits(ok, "error") || !file.exists(mdd_archivo)) {
    stop("No se pudo descargar MDD. Bájalo a mano de https://doi.org/", mdd_doi,
         " (archivo MDD_v2.5_6904species.csv) y guárdalo en: ", mdd_archivo)
  }
}
writeLines(c(paste("MDD", mdd_version), paste("DOI:", mdd_doi), paste("Archivo:", basename(mdd_archivo))),
           file.path(dir_tax, "MDD_version.txt"))

mdd <- read_csv(mdd_archivo, show_col_types = FALSE, guess_max = 10000)

requeridas <- c("sciName", "genus", "specificEpithet", "family")
faltan <- setdiff(requeridas, names(mdd))
if (length(faltan) > 0) {
  stop("Faltan columnas en MDD: ", paste(faltan, collapse = ", "),
       ". Columnas disponibles: ", paste(names(mdd), collapse = ", "))
}

# ---- 2. Lista de referencia de los géneros focales ---------------------------
opcionales <- intersect(c("subfamily", "tribe", "subgenus", "iucnStatus", "extinct", "domestic",
                          "countryDistribution", "subregionDistribution", "biogeographicRealm", "nominalNames"),
                        names(mdd))

ref <- mdd %>%
  filter(genus %in% generos_focales) %>%
  select(all_of(c("sciName", "genus", "specificEpithet", "family", opcionales))) %>%
  mutate(
    nombre_mdd  = str_replace_all(str_squish(sciName), "_", " "),
    # Epíteto sin la terminación de género gramatical (glaucus/glauca), para casar
    # nombres de GBIF y MDD que difieren solo en eso.
    epiteto_norm = str_remove(tolower(specificEpithet), "(us|a|um)$"),
    # Artibeus y Dermanura se casan juntos por epíteto (GBIF los trata como un género)
    grupo_genero = if_else(genus %in% c("Artibeus", "Dermanura"), "Artibeus_s.l.", genus),
    # DECISIÓN (30-oct-2026): Artibeus y Dermanura se analizan como UN solo género ("Artibeus").
    # MDD v2.5 los mantiene separados (13 + 11 especies); 'genus' conserva el género de MDD
    # y 'genero_analisis' es el que se usa en todos los análisis.
    genero_analisis = if_else(genus %in% c("Artibeus", "Dermanura"), "Artibeus", genus)
  )

cat("\nEspecies por género en MDD", mdd_version, ":\n")
print(ref %>% count(genus, name = "n_especies"))
cat("\nEspecies marcadas como extintas o domésticas (revisar):\n")
if (all(c("extinct", "domestic") %in% names(ref))) {
  print(ref %>% filter(extinct == 1 | domestic == 1) %>% select(nombre_mdd, extinct, domestic))
}

# Epítetos repetidos dentro de un mismo grupo: el casamiento por epíteto sería ambiguo
ambiguos <- ref %>% count(grupo_genero, epiteto_norm) %>% filter(n > 1)
if (nrow(ambiguos) > 0) {
  cat("\nAVISO: epítetos normalizados repetidos (el casamiento aproximado no se aplicará a estos):\n")
  print(ambiguos)
}

write_csv(ref, file.path(dir_tax, "lista_referencia_mdd.csv"))

# ---- 3. Especímenes de GBIF ---------------------------------------------------
focales <- readRDS(file.path(dir_gbif, "generos_focales_bruto.rds"))

esp <- focales %>%
  filter(basisOfRecord == "PRESERVED_SPECIMEN")

n_inicial <- nrow(esp)

# Rango taxonómico de los especímenes (para ver cuántos no llegan a especie)
cat("\nRango taxonómico de los especímenes (GBIF):\n")
print(esp %>% count(taxonRank, sort = TRUE))

# Solo registros con especie asignada. Se incluyen los de rango SUBSPECIES/VARIETY/FORM porque
# la columna 'species' trae la especie. Ojo: en los registros de rango GENUS/SUBGENUS la columna
# 'species' llega como cadena vacía (""), no como NA; por eso se exige que no esté vacía.
esp <- esp %>%
  filter(taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY", "FORM"),
         !is.na(species), nzchar(str_squish(species))) %>%
  mutate(
    nombre_gbif   = str_squish(str_replace_all(species, "_", " ")),
    genero_gbif   = word(nombre_gbif, 1),
    epiteto_gbif  = word(nombre_gbif, 2),
    epiteto_norm  = str_remove(tolower(epiteto_gbif), "(us|a|um)$"),
    grupo_genero  = if_else(genero_gbif %in% c("Artibeus", "Dermanura"), "Artibeus_s.l.", genero_gbif)
  )

cat("\nEspecímenes:", n_inicial, "-> identificados a especie:", nrow(esp), "\n")

# ---- 4. Casamiento en tres niveles -------------------------------------------
# 1) exacto por nombre completo; 2) mismo grupo de género + epíteto sin terminación (marcado);
# 3) sin coincidencia.
ref_exacto <- ref %>% distinct(nombre_mdd) %>% mutate(exacto = TRUE)

ref_aprox <- ref %>%
  anti_join(ambiguos, by = c("grupo_genero", "epiteto_norm")) %>%
  select(grupo_genero, epiteto_norm, nombre_aprox = nombre_mdd)

esp <- esp %>%
  left_join(ref_exacto, by = c("nombre_gbif" = "nombre_mdd")) %>%
  left_join(ref_aprox, by = c("grupo_genero", "epiteto_norm")) %>%
  mutate(
    nombre_final = case_when(
      !is.na(exacto)       ~ nombre_gbif,
      !is.na(nombre_aprox) ~ nombre_aprox,
      TRUE                 ~ NA_character_
    ),
    tipo_coincidencia = case_when(
      !is.na(exacto)       ~ "exacta",
      !is.na(nombre_aprox) ~ "aproximada (revisar)",
      TRUE                 ~ "sin coincidencia"
    ),
    excluir = FALSE
  )

# ---- 4b. Registros de rango subespecífico: respaldo en los nombres nominales de MDD ----
# GBIF asigna 'species' con su propio backbone, que puede seguir un concepto de especie anterior
# (p. ej. Sturnira lilium parvidens -> S. lilium; Anoura geoffroyi lasiopyga -> A. geoffroyi),
# mientras que MDD ya los reconoce como especies distintas. Se usa 'nominalNames' de MDD
# (epítetos de los nombres asociados a cada especie): el registro se reasigna solo si el epíteto
# subespecífico NO pertenece a la especie de GBIF y pertenece a UNA sola especie del mismo género.
if (!"nominalNames" %in% names(ref)) stop("MDD no trae la columna nominalNames.")

nominales <- ref %>%
  select(nombre_mdd, grupo_genero, nominalNames) %>%
  filter(!is.na(nominalNames), nominalNames != "NA") %>%
  mutate(nominal = str_split(nominalNames, fixed("|"))) %>%
  unnest(nominal) %>%
  mutate(
    nominal = tolower(word(str_squish(nominal), 1)),
    epiteto_sub_norm = str_remove(nominal, "(us|a|um)$")
  ) %>%
  distinct(nombre_mdd, grupo_genero, epiteto_sub_norm)

esp <- esp %>%
  mutate(
    infra = if_else(!is.na(infraspecificEpithet) & nzchar(str_squish(infraspecificEpithet)),
                    tolower(str_squish(infraspecificEpithet)), NA_character_),
    epiteto_sub_norm = str_remove(infra, "(us|a|um)$")
  )

# ¿El epíteto subespecífico ya pertenece a la especie asignada?
propios <- nominales %>%
  transmute(nombre_final = nombre_mdd, epiteto_sub_norm, soportado = TRUE)

# Especies candidatas dentro del mismo grupo de género
candidatas <- esp %>%
  filter(!is.na(infra)) %>%
  distinct(grupo_genero, epiteto_sub_norm) %>%
  inner_join(nominales, by = c("grupo_genero", "epiteto_sub_norm")) %>%
  group_by(grupo_genero, epiteto_sub_norm) %>%
  summarise(n_cand = n_distinct(nombre_mdd), cand = first(nombre_mdd), .groups = "drop")

esp <- esp %>%
  left_join(propios, by = c("nombre_final", "epiteto_sub_norm")) %>%
  left_join(candidatas, by = c("grupo_genero", "epiteto_sub_norm")) %>%
  mutate(
    nombre_antes_sub = nombre_final,
    nombre_final = case_when(
      is.na(infra) ~ nombre_final,
      !is.na(soportado) ~ nombre_final,
      n_cand %in% 1 ~ cand,
      TRUE ~ nombre_final
    ),
    tipo_coincidencia = case_when(
      is.na(infra) | !is.na(soportado) ~ tipo_coincidencia,
      n_cand %in% 1 ~ "reasignada por subespecie (revisar)",
      !is.na(n_cand) & n_cand > 1 ~ "subespecie ambigua en MDD (revisar)",
      TRUE ~ "subespecie sin respaldo en MDD (revisar)"
    )
  )

cat("\nTipo de coincidencia (número de especímenes):\n")
print(esp %>% count(tipo_coincidencia))

cat("\nRegistros reasignados por su epíteto subespecífico (nombre de GBIF -> nombre de MDD):\n")
print(esp %>%
        filter(tipo_coincidencia == "reasignada por subespecie (revisar)") %>%
        count(nombre_antes_sub, infra, nombre_final, sort = TRUE), n = 50)

# ---- 5. Tabla de revisión manual (cruce con decisiones previas) --------------
revision <- esp %>%
  filter(tipo_coincidencia != "exacta") %>%
  count(nombre_gbif, infraspecificEpithet, tipo_coincidencia, nombre_final, name = "n_especimenes") %>%
  arrange(desc(n_especimenes))

archivo_cruce <- file.path(dir_tax, "crosswalk_manual.csv")
if (!file.exists(archivo_cruce)) {
  # Plantilla: completa 'nombre_mdd' y 'accion' (aceptar / agrupar / excluir) y vuelve a correr.
  plantilla <- revision %>%
    distinct(nombre_gbif, .keep_all = TRUE) %>%
    transmute(nombre_gbif, nombre_mdd = nombre_final, accion = NA_character_, nota = NA_character_)
  write_csv(plantilla, archivo_cruce)
  cat("\nSe creó la plantilla de decisiones:", archivo_cruce,
      "\nCompleta 'accion' (aceptar / agrupar / excluir) y 'nombre_mdd', y vuelve a correr desde la sección 3.\n")
} else {
  cruce <- read_csv(archivo_cruce, show_col_types = FALSE) %>%
    distinct(nombre_gbif, .keep_all = TRUE)
  esp <- esp %>%
    left_join(cruce %>% select(nombre_gbif, nombre_mdd_manual = nombre_mdd, accion), by = "nombre_gbif") %>%
    mutate(
      nombre_final = if_else(!is.na(nombre_mdd_manual), nombre_mdd_manual, nombre_final),
      excluir = accion %in% "excluir"
    )
}

# ---- 6. Salidas y bitácora ----------------------------------------------------
esp_ok <- esp %>%
  filter(!is.na(nombre_final), !(excluir %in% TRUE)) %>%
  left_join(ref %>% select(nombre_mdd, genus_mdd = genus, genero_analisis, family_mdd = family),
            by = c("nombre_final" = "nombre_mdd"))

saveRDS(esp_ok, file.path(dir_tax, "especimenes_taxonomia_armonizada.rds"))
write_csv(revision, file.path(dir_tax, "nombres_para_revisar.csv"))

cat("\nEspecies de MDD con al menos un espécimen, por género:\n")
print(esp_ok %>%
        group_by(genero_analisis) %>%
        summarise(n_especies_con_especimenes = n_distinct(nombre_final),
                  n_especimenes = n(), .groups = "drop"))

sin_registros <- ref %>% filter(!(nombre_mdd %in% esp_ok$nombre_final)) %>% select(genus, nombre_mdd)
cat("\nEspecies de MDD SIN especímenes en GBIF (", nrow(sin_registros), "):\n", sep = "")
print(sin_registros, n = Inf)

bitacora_archivo <- file.path(dir_gbif, "bitacora_curacion.csv")
nueva <- tibble(
  paso = c("Solo PRESERVED_SPECIMEN", "Identificados a especie", "Con nombre de MDD (tras casamiento y decisiones)"),
  n_registros = c(n_inicial, nrow(esp), nrow(esp_ok))
)
bit <- if (file.exists(bitacora_archivo)) {
  bind_rows(read_csv(bitacora_archivo, show_col_types = FALSE) %>% filter(!paso %in% nueva$paso), nueva)
} else nueva
write_csv(bit, bitacora_archivo)
print(bit)
