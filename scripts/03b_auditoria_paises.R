# 03b_auditoria_paises.R
# Auditoría de países: compara el país de cada espécimen curado (script 03) con la lista de
# países de cada especie en MDD v2.5 (columna countryDistribution).
#
# Por qué va ANTES de la rejilla y del filtro de rango: los registros con nombre de taxonomía
# anterior (p. ej. "Sturnira lilium" de México o Colombia) extenderían artificialmente los rangos
# cuando se agregue la envolvente de los registros curados.
#
# Qué hace:
#   1. Aplica reasignaciones explícitas de nombre por país (reasignaciones_pais.csv).
#   2. Clasifica cada registro: listado / dudoso (MDD con "?") / no_listado / sin_pais / sin_equivalencia.
#   3. Guarda el resumen por especie y país (con candidatos a rescate) y los registros curados.
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir_tax  <- file.path(dir_base, "datos", "taxonomia")
dir_cur  <- file.path(dir_base, "datos", "curacion")

# PARÁMETROS (valores propuestos por mí: confírmalos o ajústalos)
excluir_no_listado <- TRUE   # registros de países que MDD v2.5 no lista para la especie
excluir_dudoso     <- TRUE   # registros de países que MDD lista con "?" (distribución dudosa)
# Regla (decidida el 5-nov-2026): se excluyen los registros cuestionables SOLO si la especie conserva
# al menos este número de registros tras excluirlos. Si le quedarían menos, se conservan y la
# especie pasa a "revisar_caso" para decidir a mano.
n_min_restantes <- 30

generos <- c("Sturnira", "Platyrrhinus", "Artibeus", "Dermanura", "Carollia", "Anoura")

# ---- 1. Entradas ---------------------------------------------------------------
esp <- readRDS(file.path(dir_cur, "especimenes_curados_coord.rds"))
requeridas <- c("nombre_final", "genero_analisis", "countryCode", "locality", "institutionCode",
                "elevation", "elev_dem")
faltan <- setdiff(requeridas, names(esp))
if (length(faltan) > 0) stop("Faltan columnas: ", paste(faltan, collapse = ", "))
n_entrada <- nrow(esp)
cat("\nEspecímenes de entrada:", n_entrada, "\n")

arch_mdd <- list.files(dir_tax, pattern = "MDD.*\\.csv$", ignore.case = TRUE, full.names = TRUE)
arch_mdd <- arch_mdd[!str_detect(basename(arch_mdd), "lista_referencia")]
if (length(arch_mdd) == 0) stop("No encuentro el CSV de MDD v2.5 en ", dir_tax)
if (length(arch_mdd) > 1) {
  message("Hay varios CSV de MDD; uso el primero: ", basename(arch_mdd[1]))
  arch_mdd <- arch_mdd[1]
}
mdd <- read_csv(arch_mdd, show_col_types = FALSE, guess_max = 20000)
if (!all(c("sciName", "genus", "countryDistribution") %in% names(mdd))) {
  stop("El CSV de MDD no tiene sciName, genus y countryDistribution.")
}

# Lista de países por especie. Un "?" al final indica distribución dudosa en MDD.
mdd_paises <- mdd %>%
  filter(genus %in% generos) %>%
  transmute(especie = str_replace_all(sciName, "_", " "), pais = countryDistribution) %>%
  separate_rows(pais, sep = "\\|") %>%
  mutate(pais = str_squish(pais),
         dudoso = str_detect(pais, "\\?$"),
         pais_mdd = str_squish(str_remove(pais, "\\?$"))) %>%
  filter(!is.na(pais_mdd), pais_mdd != "", pais_mdd != "NA") %>%
  distinct(especie, pais_mdd, .keep_all = TRUE) %>%
  transmute(especie, pais_mdd, dudoso, en_lista = TRUE)

# ---- 2. Equivalencias código de país (GBIF, ISO-2) -> nombre en MDD ------------
equiv <- tribble(
  ~countryCode, ~pais_mdd,
  "MX", "Mexico", "GT", "Guatemala", "BZ", "Belize", "HN", "Honduras",
  "SV", "El Salvador", "NI", "Nicaragua", "CR", "Costa Rica", "PA", "Panama",
  "CO", "Colombia", "VE", "Venezuela", "GY", "Guyana", "SR", "Suriname",
  "GF", "French Guiana", "EC", "Ecuador", "PE", "Peru", "BO", "Bolivia",
  "BR", "Brazil", "PY", "Paraguay", "AR", "Argentina", "UY", "Uruguay", "CL", "Chile",
  "TT", "Trinidad and Tobago", "GD", "Grenada", "VC", "Saint Vincent and the Grenadines",
  "LC", "Saint Lucia", "DM", "Dominica", "MQ", "Martinique", "GP", "Guadeloupe",
  "BB", "Barbados", "JM", "Jamaica", "CU", "Cuba", "HT", "Haiti", "DO", "Dominican Republic",
  "PR", "Puerto Rico", "BS", "Bahamas", "AG", "Antigua and Barbuda",
  "KN", "Saint Kitts and Nevis", "US", "United States",
  "KY", "Cayman Islands", "MS", "Montserrat"
)
# Sin equivalencia en MDD para estas especies (quedan como "sin_equivalencia" y no se excluyen):
# VI (Islas Vírgenes de EE. UU.), TC (Islas Turks y Caicos), CW (Curazao).
# Nombres de MDD (para estas especies) que no están en la tabla: revisa que ninguno sea relevante.
sin_equiv <- setdiff(unique(mdd_paises$pais_mdd), equiv$pais_mdd)
if (length(sin_equiv) > 0) {
  message("Países de MDD sin equivalencia en la tabla: ", paste(sin_equiv, collapse = ", "))
}

# ---- 3. Reasignaciones explícitas de nombre por país --------------------------
# Solo donde hay respaldo en la literatura. Se aplican ANTES de comparar con MDD.
arch_reas <- file.path(dir_cur, "reasignaciones_pais.csv")
reas_defecto <- bind_rows(
  tibble(especie_original = "Anoura geoffroyi",
         countryCode = c("MX", "GT", "BZ", "HN", "SV", "NI", "CR", "PA"),
         especie_nueva = "Anoura peruana",
         fuente = "Molinari et al. 2023: las poblaciones de México y Centroamérica son A. peruana lasiopyga"),
  # Panamá queda fuera: MDD lista C. brevicauda y C. sowelli (el límite está en el oeste de Panamá).
  tibble(especie_original = "Carollia brevicauda",
         countryCode = c("MX", "BZ", "GT", "HN", "NI", "CR"),
         especie_nueva = "Carollia sowelli",
         fuente = "Baker, Solari y Hoffmann 2002 (Occas. Pap. Mus. Texas Tech Univ.): C. sowelli, de Panamá occidental a Veracruz, San Luis Potosí y Oaxaca"),
  tibble(especie_original = "Sturnira mistratensis",
         countryCode = "PE",
         especie_nueva = "Sturnira giannae",
         fuente = "de Oliveira et al. 2024, Zoologia 41:e23061 (apéndice): FMNH 203416, Moyobamba (San Martín), es S. giannae"),
  tibble(especie_original = "Sturnira lilium",
         countryCode = "MX",
         especie_nueva = "Sturnira parvidens",
         fuente = "Velazco y Patterson 2013, Mol. Phylogenet. Evol. 68:683-698: ejemplares secuenciados de México (Sonora a Chiapas, Campeche, Quintana Roo) son S. parvidens, resucitada de la sinonimia de S. lilium")
)
# Se combinan con las reglas que ya estén en el archivo (las tuyas tienen prioridad).
reas_prev <- if (file.exists(arch_reas)) read_csv(arch_reas, show_col_types = FALSE) else NULL
reas_todas <- bind_rows(reas_prev, reas_defecto) %>%
  distinct(especie_original, countryCode, .keep_all = TRUE)
write_csv(reas_todas, arch_reas)
reas <- reas_todas %>% select(especie_original, countryCode, especie_nueva)

esp <- esp %>%
  select(-any_of(c("especie_nueva", "nombre_auditado", "reasignado", "pais_mdd",
                   "dudoso", "en_lista", "estatus_pais"))) %>%
  left_join(reas, by = c("nombre_final" = "especie_original", "countryCode")) %>%
  mutate(nombre_auditado = coalesce(especie_nueva, nombre_final),
         reasignado = !is.na(especie_nueva)) %>%
  select(-especie_nueva)

# Fusiones taxonómicas del análisis principal (decidido el 5-nov-2026).
# Artibeus intermedius se trata como parte de A. lituratus: MDD v2.5 la reconoce, pero su propia nota
# dice que casi siempre se considera sinónimo de A. lituratus y que las diferencias moleculares y
# morfológicas son limitadas. En la sensibilidad se mantiene separada (nombre_sensibilidad).
fusiones <- tibble(nombre_origen = "Artibeus intermedius", nombre_fusion = "Artibeus lituratus")
esp <- esp %>%
  mutate(nombre_pre_fusion = nombre_auditado) %>%
  left_join(fusiones, by = c("nombre_auditado" = "nombre_origen")) %>%
  mutate(nombre_auditado = coalesce(nombre_fusion, nombre_auditado)) %>%
  select(-nombre_fusion)
cat("Registros fusionados:", sum(esp$nombre_pre_fusion != esp$nombre_auditado), "\n")

# Anoura geoffroyi en los países andinos: regla por elevación (decidida el 5-nov-2026).
# Molinari et al. 2023 (Zootaxa 5297): las poblaciones andinas sudamericanas son A. peruana peruana,
# mientras que A. geoffroyi es de tierras bajas (típicamente < 1,000-1,500 m); en la pendiente
# amazónica peruana hay ejemplares de ambas a altitudes parecidas. Los umbrales son una aproximación
# mía, no un criterio de los autores:
#   > 1500 m -> A. peruana | 1000-1500 m o sin elevación -> se excluye (ambiguo) | < 1000 m -> A. geoffroyi
andinos <- c("CO", "EC", "PE", "BO")
elev_alta <- 1500; elev_baja <- 1000
esp <- esp %>%
  mutate(elev_eval = coalesce(elev_dem, suppressWarnings(as.numeric(elevation))),
         geo_andino = nombre_auditado == "Anoura geoffroyi" & countryCode %in% andinos,
         clase_geo = case_when(
           !geo_andino ~ NA_character_,
           is.na(elev_eval) ~ "ambiguo_sin_elevacion",
           elev_eval > elev_alta ~ "peruana",
           elev_eval >= elev_baja ~ "ambiguo_1000_1500",
           TRUE ~ "geoffroyi"),
         nombre_auditado = if_else(!is.na(clase_geo) & clase_geo == "peruana", "Anoura peruana", nombre_auditado),
         reasignado = reasignado | (!is.na(clase_geo) & clase_geo == "peruana"),
         excl_elev_geoffroyi = !is.na(clase_geo) & clase_geo %in% c("ambiguo_sin_elevacion", "ambiguo_1000_1500"))
cat("\nA. geoffroyi andino (CO, EC, PE, BO) por elevación:\n")
print(esp %>% filter(geo_andino) %>% count(countryCode, clase_geo) %>%
        pivot_wider(names_from = clase_geo, values_from = n, values_fill = 0))
esp <- esp %>% select(-geo_andino)
cat("Registros reasignados:", sum(esp$reasignado), "\n")
print(esp %>% filter(reasignado) %>% count(nombre_final, nombre_auditado, countryCode))

# ---- 4. Clasificación de cada registro ----------------------------------------
esp <- esp %>%
  left_join(equiv, by = "countryCode") %>%
  left_join(mdd_paises, by = c("nombre_auditado" = "especie", "pais_mdd")) %>%
  mutate(
    en_lista = coalesce(en_lista, FALSE),
    dudoso   = coalesce(dudoso, FALSE),
    estatus_pais = case_when(
      is.na(countryCode) | countryCode == "" ~ "sin_pais",
      is.na(pais_mdd)                        ~ "sin_equivalencia",
      en_lista & dudoso                      ~ "dudoso",
      en_lista                               ~ "listado",
      TRUE                                   ~ "no_listado"
    )
  )

cat("\nRegistros por estatus de país:\n")
print(esp %>% count(estatus_pais) %>% mutate(porcentaje = round(100 * n / sum(n), 1)))
cat("\nCódigos de país sin equivalencia (revisa si hay que agregarlos):\n")
print(esp %>% filter(estatus_pais == "sin_equivalencia") %>% count(countryCode, sort = TRUE))

# ---- 5. Decisión por registro: excluir solo si la especie conserva suficientes ------
esp <- esp %>%
  mutate(cuestionable = (excluir_no_listado & estatus_pais == "no_listado") |
                        (excluir_dudoso & estatus_pais == "dudoso"))

restantes <- esp %>%
  group_by(nombre_auditado) %>%
  summarise(n_total = n(), n_cuestionables = sum(cuestionable),
            n_restantes = n_total - n_cuestionables, .groups = "drop")

esp <- esp %>%
  left_join(restantes %>% select(nombre_auditado, n_restantes), by = "nombre_auditado") %>%
  mutate(decision_pais = case_when(
    !cuestionable                  ~ "conservar",
    n_restantes >= n_min_restantes ~ "excluir",
    TRUE                           ~ "revisar_caso"   # se conserva hasta decidir a mano
  ))

# Decisiones manuales (por especie y país) sobre los registros cuestionables. Tienen prioridad
# sobre la regla del mínimo de registros. Edita este archivo para cambiarlas.
arch_man <- file.path(dir_cur, "decisiones_manuales_pais.csv")
if (!file.exists(arch_man)) {
  write_csv(tibble(
    nombre_auditado = c("Platyrrhinus guianensis", "Sturnira sorianoi", "Platyrrhinus aquilus"),
    countryCode     = c("BR", "CO", "CO"),
    decision        = c("excluir", "excluir", "conservar"),
    motivo = c("17 registros sin localidad ni institución, lejos del Escudo Guayanés; MDD no lista Brasil",
               "MDD no lista Colombia",
               "Cerro Tacarcuna (frontera con Panamá), altitud coherente; MDD pone Colombia con '?'")
  ), arch_man)
}
man <- read_csv(arch_man, show_col_types = FALSE) %>%
  distinct(nombre_auditado, countryCode, .keep_all = TRUE) %>%
  select(nombre_auditado, countryCode, decision_manual = decision)
esp <- esp %>%
  left_join(man, by = c("nombre_auditado", "countryCode")) %>%
  mutate(decision_pais = if_else(cuestionable & !is.na(decision_manual),
                                 decision_manual, decision_pais)) %>%
  select(-decision_manual)

# Registros ambiguos de A. geoffroyi andino (regla por elevación): se excluyen en principal y sensibilidad
esp <- esp %>% mutate(decision_pais = if_else(excl_elev_geoffroyi, "excluir", decision_pais))

cat("\nDecisión por registro:\n")
print(esp %>% count(decision_pais) %>% mutate(porcentaje = round(100 * n / sum(n), 1)))

casos <- esp %>% filter(decision_pais == "revisar_caso") %>%
  count(nombre_auditado, countryCode, estatus_pais, name = "n_registros")
cat("\nCasos pendientes (a la especie le quedarían menos de", n_min_restantes,
    "registros y no hay decisión manual):\n")
print(casos)

# Reasignaciones SOLO para la sensibilidad (decidido el 5-nov-2026): en el análisis principal esos
# registros siguen excluidos; en la sensibilidad entran con el nombre reasignado.
arch_sens <- file.path(dir_cur, "reasignaciones_sensibilidad.csv")
if (!file.exists(arch_sens)) {
  write_csv(tibble(
    especie_original = "Sturnira lilium",
    countryCode = c("VE", "GY", "TT", "SR", "GF"),
    especie_nueva = "Sturnira giannae",
    fuente = "Velazco y Patterson 2014, 2019: la especie confundida durante años con S. lilium en estos países es S. giannae; asignación incierta por registro (también ocurre S. tildae)"
  ), arch_sens)
}
sens <- read_csv(arch_sens, show_col_types = FALSE) %>%
  distinct(especie_original, countryCode, .keep_all = TRUE) %>%
  select(especie_original, countryCode, especie_sens = especie_nueva)
esp <- esp %>%
  select(-any_of("especie_sens")) %>%
  left_join(sens, by = c("nombre_auditado" = "especie_original", "countryCode")) %>%
  mutate(nombre_sensibilidad = coalesce(especie_sens, nombre_auditado),
         # en la sensibilidad las especies fusionadas vuelven a separarse
         nombre_sensibilidad = if_else(is.na(especie_sens) & nombre_pre_fusion != nombre_auditado,
                                       nombre_pre_fusion, nombre_sensibilidad)) %>%
  select(-especie_sens)
cat("\nRegistros con reasignación solo en la sensibilidad:", sum(esp$nombre_sensibilidad != esp$nombre_auditado), "\n")

# ---- 6. Resumen y salidas ------------------------------------------------------
resumen <- esp %>%
  group_by(nombre_auditado, countryCode, pais_mdd, estatus_pais, decision_pais) %>%
  summarise(n = n(),
            n_localidades = n_distinct(locality),
            n_instituciones = n_distinct(institutionCode),
            .groups = "drop") %>%
  arrange(desc(decision_pais != "conservar"), desc(n))

por_especie <- esp %>%
  group_by(nombre_auditado) %>%
  summarise(n_total = n(),
            n_no_listado = sum(estatus_pais == "no_listado"),
            n_dudoso = sum(estatus_pais == "dudoso"),
            n_excluidos = sum(decision_pais == "excluir"),
            n_restantes = sum(decision_pais != "excluir"),
            pct_no_listado = round(100 * n_no_listado / n_total, 1),
            .groups = "drop") %>%
  arrange(desc(n_excluidos))

cat("\nEspecies con más registros excluidos:\n")
print(head(por_especie, 15))

write_csv(resumen, file.path(dir_cur, "auditoria_paises_resumen.csv"))
write_csv(por_especie, file.path(dir_cur, "auditoria_paises_por_especie.csv"))
write_csv(esp %>% filter(decision_pais == "revisar_caso") %>%
            select(any_of(c("gbifID", "nombre_final", "nombre_auditado", "countryCode", "stateProvince",
                            "locality", "decimalLatitude", "decimalLongitude", "institutionCode",
                            "catalogNumber", "estatus_pais"))),
          file.path(dir_cur, "auditoria_paises_casos_a_revisar.csv"))
saveRDS(esp, file.path(dir_cur, "especimenes_con_pais.rds"))

curados <- esp %>% filter(decision_pais != "excluir")
saveRDS(curados, file.path(dir_cur, "especimenes_curados_pais.rds"))
# Base para la sensibilidad: lo conservado más los registros excluidos que se reasignan en ella.
# Usa la columna "nombre_sensibilidad".
saveRDS(esp %>% filter(decision_pais != "excluir" | nombre_sensibilidad != nombre_auditado),
        file.path(dir_cur, "especimenes_sensibilidad_pais.rds"))
cat("\nConservados tras la auditoría de países:", nrow(curados), "de", n_entrada, "\n")

bitacora_archivo <- file.path(dir_gbif, "bitacora_curacion.csv")
nueva <- tibble(
  paso = c("Auditoría de países (excluye países no listados o dudosos en MDD v2.5)",
           "Registros de especies a revisar por caso (conservados)"),
  n_registros = c(nrow(curados), sum(esp$decision_pais == "revisar_caso"))
)
bit <- if (file.exists(bitacora_archivo)) {
  # también se quitan las filas de la versión anterior de este script
  viejas <- c("Auditoría de países (excluye países no listados en MDD v2.5)",
              "Registros de países dudosos en MDD (conservados)")
  bind_rows(read_csv(bitacora_archivo, show_col_types = FALSE) %>%
              filter(!paso %in% c(nueva$paso, viejas)), nueva)
} else nueva
write_csv(bit, bitacora_archivo)
print(bit)

# A partir de aquí, los análisis usan "nombre_auditado" (no "nombre_final").
