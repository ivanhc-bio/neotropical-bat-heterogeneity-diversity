# 00_recorte_trace_remoto.R
# Descarga SOLO la ventana de interés de CHELSA-TraCE21k (bio01 y bio12) leyendo los GeoTIFF remotos con GDAL
# (/vsicurl/), sin bajar los archivos globales completos, y la agrega a ~5 km al vuelo.
# Resultado: unas 90 capas pequeñas (unos MB cada una) en <RUTA_TRACE>/bio01 y /bio12, con los
# mismos nombres que los originales, listas para 05b_historia_climatica.R.
#
# Fuente: EnviDat / CHELSA (Karger et al. 2023), https://www.doi.org/10.16904/envidat.211
# URLs tomadas de tu archivo envidatS3paths.txt (ruta S3 oficial de EnviDat):
#   https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim/bio01/CHELSA_TraCE21k_bio01_-200_V.1.0.tif
# El timeID va con 3 dígitos si es negativo (-001, -099, -200) y con 4 si es 0 o positivo (0000, 0020).
# No probado: desde mi entorno no hay acceso a ese servidor. Corre primero con modo_prueba <- TRUE.
#
# Es reanudable: salta lo que ya existe, y se puede cortar y volver a correr.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(terra)

# ---- Parámetros ----------------------------------------------------------------
modo_prueba <- FALSE     # TRUE: solo timeID 20 y 0 de bio01, para medir tiempo y tamaño. Luego pon FALSE.
dir_salida  <- RUTA_TRACE
dir_tmp     <- RUTA_TMP_TERRA
ext_ventana <- ext(-118, -34, -38, 33)     # xmin, xmax, ymin, ymax (grados); cubre todas las celdas + margen
fact_agregacion <- 5                       # 30 arcsec (~1 km) -> ~5 km; 05b debe usar fact_agregacion <- 1
time_ids <- seq(-200, 20, by = 5)          # 45 capas, cada 500 años (1 unidad = 100 años; 20 = 1900-1990)
vars <- list(
  bio01 = list(relleno_mayor = NULL),
  bio12 = list(relleno_mayor = 6500)       # igual que en 05b; si el máximo real pasa de 6500, súbelo
)
base_url <- "https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim"
fmt_id <- function(id) ifelse(id < 0, sprintf("-%03d", abs(id)), sprintf("%04d", id))
intentos <- 3

dir.create(dir_tmp, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = dir_tmp, memfrac = 0.5, todisk = FALSE)
# Lectura remota robusta: reintentos, sin listar el directorio remoto, caché de bloques
setGDALconfig("GDAL_HTTP_MAX_RETRY", "6")
setGDALconfig("GDAL_HTTP_RETRY_DELAY", "5")
setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
setGDALconfig("GDAL_HTTP_TIMEOUT", "120")
setGDALconfig("VSI_CACHE", "TRUE")

if (modo_prueba) { time_ids <- c(20, 0); vars <- vars["bio01"] }

# ---- Función: una capa -----------------------------------------------------------
procesar_capa <- function(v, id, cfg) {
  nom <- sprintf("CHELSA_TraCE21k_%s_%s_V.1.0.tif", v, fmt_id(id))
  url <- file.path(base_url, v, nom)
  dir_v <- file.path(dir_salida, v); dir.create(dir_v, showWarnings = FALSE, recursive = TRUE)
  salida <- file.path(dir_v, nom)
  if (file.exists(salida)) return(invisible("ya existe"))
  tmp_crop <- file.path(dir_tmp, paste0("crop_", nom))
  t0 <- Sys.time()
  r <- rast(paste0("/vsicurl/", url))                       # no descarga nada todavía
  x <- crop(r, ext_ventana, filename = tmp_crop, overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))      # lee solo la ventana, por bloques
  if (!is.null(cfg$relleno_mayor)) {
    x2 <- file.path(dir_tmp, paste0("relleno_", nom))
    x <- ifel(x > cfg$relleno_mayor, NA, x, filename = x2, overwrite = TRUE)
  }
  if (fact_agregacion > 1) {
    x <- aggregate(x, fact = fact_agregacion, fun = "mean", na.rm = TRUE,
                   filename = paste0(salida, ".part.tif"), overwrite = TRUE,
                   gdal = c("COMPRESS=DEFLATE"))
  } else {
    x <- writeRaster(x, paste0(salida, ".part.tif"), overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  }
  rm(x, r); gc(verbose = FALSE)
  file.rename(paste0(salida, ".part.tif"), salida)
  unlink(list.files(dir_tmp, pattern = paste0("(crop|relleno)_", gsub("\\.", "\\\\.", nom)), full.names = TRUE))
  cat(sprintf("  %s listo en %.1f min (%.1f MB)\n", nom, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              file.size(salida) / 1e6))
  invisible("ok")
}

# ---- Bucle con reintentos ---------------------------------------------------------
fallos <- list()
for (v in names(vars)) {
  cat("\n==", v, "==\n")
  for (id in time_ids) {
    ok <- FALSE
    for (k in seq_len(intentos)) {
      res <- tryCatch(procesar_capa(v, id, vars[[v]]), error = function(e) {
        message("  intento ", k, " falló para ", v, " ", id, ": ", conditionMessage(e)); NULL })
      if (!is.null(res)) { ok <- TRUE; break }
      Sys.sleep(10)
    }
    if (!ok) fallos[[length(fallos) + 1]] <- paste(v, id)
  }
}
cat("\nTerminado. Capas fallidas:", if (length(fallos)) paste(unlist(fallos), collapse = "; ") else "ninguna", "\n")

# ---- Revisión rápida ----------------------------------------------------------------
fs <- list.files(dir_salida, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
fs <- fs[!grepl("\\.part\\.tif$", fs)]
cat("Archivos en", dir_salida, ":", length(fs), "\n")
if (length(fs)) {
  r1 <- rast(fs[1]); print(ext(r1)); print(res(r1))
  print(global(r1, c("min", "mean", "max"), na.rm = TRUE))
}
# Cuando termines con modo_prueba <- TRUE y todo se vea bien (extensión, valores, tiempo), cambia a FALSE y vuelve a correr.
