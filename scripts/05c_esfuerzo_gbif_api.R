# 05c_esfuerzo_gbif_api.R
# Esfuerzo de muestreo INDEPENDIENTE de los géneros focales, por celda LAEA (50, 100 y 200 km), contando directamente en
# los servidores de GBIF (petición a la API con limit = 0: no se descarga ningún registro). Es la adaptación a la rejilla LAEA
# de 03_sampling_effort.R del manuscrito anterior (que usaba la rejilla Behrmann).
#
# Variables (conteo n y log10(n + 1)); mismo periodo que el borrador (1960-2025), con coordenadas y sin problemas
# geoespaciales:
#   esf_mam    : todos los registros de Mammalia (idéntico a la definición del borrador)
#   esf_quir   : Chiroptera menos los cinco géneros focales (Anoura, Artibeus, Dermanura, Carollia, Platyrrhinus,
#                Sturnira): muestreo específico de murciélagos
#   esf_quir_ej: igual, solo ejemplares de museo (PRESERVED_SPECIMEN), el mismo tipo de dato que los registros focales
# Para Chiroptera se pide un solo conteo con facetas por género (facet = "genusKey") y se restan los géneros focales.
#
# Solo se consultan las celdas incluidas en algún criterio (script 06; columnas incl_*), no toda la rejilla.
# Reanuda desde un checkpoint si se interrumpe. No necesita cuenta de GBIF.
#
# No probado: no tengo R ni acceso a GBIF en mi entorno. Prueba con n_prueba <- 5 antes de lanzar todo, y revisa que
# esf_quir salga distinto de esf_mam y sin NA (si las facetas no vienen en el resultado, el script lo avisa).
# Requiere: rgbif (solo para name_backbone), crul, jsonlite, sf, dplyr, readr (crul y jsonlite vienen con rgbif)

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(rgbif)
library(sf)
library(dplyr)
library(readr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_pred <- file.path(dir_base, "datos", "predictores")
dir_chk  <- file.path(dir_base, "datos", "gbif_esfuerzo")
dir.create(dir_chk, showWarnings = FALSE, recursive = TRUE)

crs_laea   <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"   # igual que 04
tamanos_km <- c(100, 50, 200)
anios      <- "1960,2025"
generos_focales <- c("Anoura", "Artibeus", "Dermanura", "Carollia", "Platyrrhinus", "Sturnira")
n_prueba   <- Inf        # pon 5 para probar con las primeras 5 celdas de cada tamaño
pausa_s    <- 0.25       # pausa entre consultas, por cortesía con la API

# Consulta directa a la API de GBIF (sin rgbif::occ_search, que en versiones recientes limita las consultas masivas y
# devolvió ceros). Se usa crul, que rgbif ya instala. Reintenta con espera si GBIF responde 429 o error de servidor.
cliente <- crul::HttpClient$new("https://api.gbif.org", headers = list(`User-Agent` = "esfuerzo-muestreo-tesis"))
pedir <- function(query, intentos = 6) {
  for (k in seq_len(intentos)) {
    r <- tryCatch(cliente$get("v1/occurrence/search", query = query), error = function(e) NULL)
    if (!is.null(r) && r$status_code == 200) return(jsonlite::fromJSON(r$parse("UTF-8"), simplifyVector = FALSE))
    if (!is.null(r) && r$status_code == 400) stop("GBIF rechazó la consulta (400): ", substr(r$parse("UTF-8"), 1, 200))
    Sys.sleep(2^k / 4)                       # espera creciente: 0.5, 1, 2, 4, 8, 16 s
  }
  stop("GBIF no respondió tras ", intentos, " intentos")
}

# ---- 1. Claves taxonómicas ------------------------------------------------------------
# Se piden directamente a la API (species/match); name_backbone() cambió el nombre de sus columnas entre versiones de rgbif
# y devolvió NA. Si una clave no se resuelve, el script se detiene.
cliente_sp <- crul::HttpClient$new("https://api.gbif.org")
# Los géneros se buscan dentro del orden Chiroptera: Anoura y Platyrrhinus también son nombres de géneros de insectos y, sin
# restringir, GBIF devolvía el reino (clave 1). Se exige que el resultado sea exactamente del rango pedido.
clave <- function(nombre, rango, orden = NULL) {
  q <- list(name = nombre, rank = rango, kingdom = "Animalia", strict = "true")
  if (!is.null(orden)) q$order <- orden
  r <- cliente_sp$get("v1/species/match", query = q)
  x <- jsonlite::fromJSON(r$parse("UTF-8"))
  k <- suppressWarnings(as.integer(x$usageKey))
  ok <- r$status_code == 200 && length(k) == 1 && !is.na(k) && identical(toupper(x$rank), toupper(rango)) &&
    identical(x$matchType, "EXACT")
  if (!ok) stop("No se pudo resolver la clave de ", nombre, " (", rango, "): matchType = ", x$matchType, ", rank = ", x$rank)
  k
}
key_mam   <- clave("Mammalia", "class")
key_chiro <- clave("Chiroptera", "order")

# Claves de los géneros focales sin buscar por nombre: se piden a GBIF las facetas por género de Chiroptera (las claves que usa
# el propio índice de registros) y se lee el nombre de cada una. Así no dependen del emparejamiento de nombres, que falla con
# homónimos (Anoura, Platyrrhinus también son géneros de insectos) y con los identificadores nuevos del backbone.
fac0 <- pedir(list(taxonKey = key_chiro, hasCoordinate = "true", hasGeospatialIssue = "false", year = anios,
                   limit = 0, facet = "genusKey", facetLimit = 1000))
cl_gen <- vapply(fac0$facets[[1]]$counts, function(z) as.character(z$name), character(1))
n_gen  <- vapply(fac0$facets[[1]]$counts, function(z) as.numeric(z$count), numeric(1))
nom_gen <- vapply(cl_gen, function(k) {
  r <- cliente_sp$get(paste0("v1/species/", k))
  x <- if (r$status_code == 200) jsonlite::fromJSON(r$parse("UTF-8")) else list()
  nm <- if (!is.null(x$canonicalName)) x$canonicalName else x$scientificName
  if (is.null(nm)) NA_character_ else sub(" .*$", "", as.character(nm))   # solo el nombre del género
}, character(1))
tabla_gen <- data.frame(clave = cl_gen, genero = nom_gen, registros = n_gen)
keys_focales <- tabla_gen$clave[tabla_gen$genero %in% generos_focales]
cat("Géneros focales encontrados en las facetas de Chiroptera:\n")
print(tabla_gen[tabla_gen$genero %in% generos_focales, ], row.names = FALSE)
obligatorios <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")   # Dermanura puede estar dentro de Artibeus
faltan <- setdiff(obligatorios, tabla_gen$genero)
if (length(faltan)) stop("No se encontraron en las facetas: ", paste(faltan, collapse = ", "))
cat("Claves: Mammalia", key_mam, "(esperado 359) | Chiroptera", key_chiro, "\n")
if (key_mam != 359L) stop("La clave de Mammalia no es 359; revisa la resolución de nombres")

# ---- 2. Funciones -------------------------------------------------------------------
# polígono de la celda en grados (WKT); se densifican los bordes porque la celda LAEA no es un rectángulo en lon/lat;
# GBIF exige el anillo exterior en sentido antihorario
wkt_celda <- function(geom, tam_km) {
  g <- st_transform(st_segmentize(st_sfc(geom, crs = crs_laea), tam_km * 1000 / 5), 4326)[[1]]
  m <- st_coordinates(g)[, 1:2]
  area <- sum(m[-nrow(m), 1] * m[-1, 2] - m[-1, 1] * m[-nrow(m), 2]) / 2
  if (area < 0) m <- m[nrow(m):1, ]
  paste0("POLYGON((", paste(sprintf("%.5f %.5f", m[, 1], m[, 2]), collapse = ", "), "))")
}

conteo <- function(wkt, taxon, basis = NULL, faceta = FALSE) {
  q <- list(taxonKey = taxon, geometry = wkt, hasCoordinate = "true", hasGeospatialIssue = "false",
            year = anios, limit = 0)
  if (!is.null(basis)) q$basisOfRecord <- basis
  if (faceta) { q$facet <- "genusKey"; q$facetLimit <- 1000 }
  res <- pedir(q)
  total <- res$count
  if (is.null(total)) stop("La respuesta de GBIF no trae 'count'")
  if (!faceta) return(total)
  fac <- if (length(res$facets)) res$facets[[1]]$counts else list()
  focal <- sum(vapply(fac, function(z) if (as.character(z$name) %in% as.character(keys_focales)) z$count else 0, numeric(1)))
  c(total = total, focal = focal)
}

consultar_celda <- function(wkt) {
  tryCatch({
    mam <- conteo(wkt, key_mam)
    q <- conteo(wkt, key_chiro, faceta = TRUE)
    qe <- conteo(wkt, key_chiro, basis = "PRESERVED_SPECIMEN", faceta = TRUE)
    tibble(esf_mam_n = mam, esf_quir_n = q[["total"]] - q[["focal"]], esf_quir_ej_n = qe[["total"]] - qe[["focal"]])
  }, error = function(e) {
    warning("Falló una celda: ", conditionMessage(e))
    tibble(esf_mam_n = NA_real_, esf_quir_n = NA_real_, esf_quir_ej_n = NA_real_)
  })
}

# ---- 2b. Prueba de humo: un rectángulo sobre México central debe tener miles de registros de Mammalia ----
prueba <- conteo("POLYGON((-100 15, -90 15, -90 20, -100 20, -100 15))", key_mam)
cat("Prueba de humo (rectángulo sobre el sur de México): Mammalia =", prueba, "\n")
if (!is.finite(prueba) || prueba < 1000) stop("La API devolvió un conteo implausible; no sigas hasta revisar esto.")

# ---- 3. Un tamaño de celda a la vez -----------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  inc <- read_csv(file.path(dir_mod, sprintf("celdas_incluidas_%dkm.csv", tam_km)), show_col_types = FALSE)
  cols_inc <- grep("^incl_", names(inc), value = TRUE)
  usadas <- inc$celda[rowSums(inc[, cols_inc], na.rm = TRUE) > 0]
  rej <- st_read(file.path(dir_rej, sprintf("celdas_%dkm.gpkg", tam_km)), quiet = TRUE) %>% filter(celda %in% usadas)
  cat("Celdas a consultar:", nrow(rej), "\n")
  if (is.finite(n_prueba)) rej <- rej[seq_len(min(n_prueba, nrow(rej))), ]

  archivo_chk <- file.path(dir_chk, sprintf("esfuerzo_parcial_v2_%dkm.rds", tam_km))
  res <- if (file.exists(archivo_chk)) readRDS(archivo_chk) else tibble(celda = character(), esf_mam_n = numeric(),
                                                                         esf_quir_n = numeric(), esf_quir_ej_n = numeric())
  hechas <- res$celda[!is.na(res$esf_mam_n) & !is.na(res$esf_quir_n) & !is.na(res$esf_quir_ej_n)]
  res <- res %>% filter(celda %in% hechas)
  pendientes <- setdiff(rej$celda, hechas)
  cat("Ya hechas:", length(hechas), "| pendientes:", length(pendientes), "\n")

  for (i in seq_along(pendientes)) {
    cid <- pendientes[i]
    wkt <- wkt_celda(st_geometry(rej)[[which(rej$celda == cid)]], tam_km)
    res <- bind_rows(res, consultar_celda(wkt) %>% mutate(celda = cid, .before = 1))
    Sys.sleep(pausa_s)
    if (i %% 20 == 0 || i == length(pendientes)) {
      cat(sprintf("  [%d/%d] %s -> Mammalia %s | Chiroptera sin focales %s | ejemplares %s\n", i, length(pendientes), cid,
                  tail(res$esf_mam_n, 1), tail(res$esf_quir_n, 1), tail(res$esf_quir_ej_n, 1)))
      saveRDS(res, archivo_chk)
    }
  }
  saveRDS(res, archivo_chk)

  fallos <- sum(is.na(res$esf_mam_n) | is.na(res$esf_quir_n) | is.na(res$esf_quir_ej_n))
  if (fallos > 0) cat("AVISO:", fallos, "celdas con NA (fallos de la API o facetas ausentes); vuelve a correr para reintentar\n")
  salida <- res %>%
    mutate(esf_mam_log = log10(esf_mam_n + 1), esf_quir_log = log10(esf_quir_n + 1), esf_quir_ej_log = log10(esf_quir_ej_n + 1))
  write_csv(salida, file.path(dir_pred, sprintf("esfuerzo_independiente_%dkm.csv", tam_km)))
  if (median(res$esf_mam_n, na.rm = TRUE) == 0) stop("La mediana de registros de Mammalia por celda es 0: algo falla en la consulta")
  cat(sprintf("Medianas de registros por celda: Mammalia %.0f | Chiroptera sin focales %.0f | ejemplares %.0f\n",
              median(res$esf_mam_n, na.rm = TRUE), median(res$esf_quir_n, na.rm = TRUE), median(res$esf_quir_ej_n, na.rm = TRUE)))
}

# ---- 4. Comparación con el esfuerzo focal de 05 (solo diagnóstico, 100 km) ----------------------------
f100 <- file.path(dir_pred, "predictores_100km.csv")
e100 <- file.path(dir_pred, "esfuerzo_independiente_100km.csv")
if (file.exists(f100) && file.exists(e100)) {
  e <- read_csv(e100, show_col_types = FALSE) %>%
    inner_join(read_csv(f100, show_col_types = FALSE) %>% select(celda, ctrl_esfuerzo_log), by = "celda")
  cat("\nCorrelación de Spearman entre el esfuerzo focal (05) y los independientes (100 km):\n")
  print(round(cor(e[, c("ctrl_esfuerzo_log", "esf_mam_log", "esf_quir_log", "esf_quir_ej_log")],
                  method = "spearman", use = "complete.obs"), 2))
}
cat("\nListo. Después corre 09 y 10 (esfuerzo_independiente en 09 indica qué variable usar).\n")
