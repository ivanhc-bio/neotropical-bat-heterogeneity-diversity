# 05b_historia_climatica.R
# Anomalía y velocidad del cambio climático desde el LGM a partir de las series de CHELSA-TraCE21k
# (capas cada 500 años), resumidas por celda de 50, 100 y 200 km.
#
# Métricas por píxel (temperatura bio01 y precipitación bio12):
#   hist_anom_temp        |T presente - T LGM| (°C)
#   hist_anom_prec        |log((P presente + 1) / (P LGM + 1))| (cambio relativo)
#   hist_vel_neta_log10   log10 de la velocidad neta LGM -> presente (m/año): (|dT| / años) / gradiente espacial
#   hist_vel_media_log10  media, sobre los intervalos de 500 años, del log10 de la velocidad de cada
#                         intervalo; gradiente espacial de la media de las dos capas del intervalo
#   hist_sd_temp          SD de la temperatura a lo largo de toda la serie (variabilidad total)
#   hist_sd_prec_mm       SD de la precipitación (mm) a lo largo de toda la serie. Es la definición del borrador 6
#                         (bio12_sd_hist); crece con la precipitación media, así que se confunde con ella
#   hist_sd_prec_log      SD de log(P + 1) a lo largo de la serie: variabilidad relativa, sin depender de la media.
#                         Es la candidata principal para la hipótesis de inestabilidad de precipitación
# Velocidad = tasa temporal (°C/año) / gradiente espacial (°C/km), como en Loarie et al. (2009).
# El paquete VoCC es la implementación de referencia; aquí se hace con terra para controlar cada paso.
#
# Se trabaja en la proyección de las celdas (LAEA) con píxeles de ~5 km: las capas de 1 km de una
# serie de ~44 pasos no caben en memoria y para celdas de 50-200 km no hace falta tanto detalle.
# Limitación: a 5 km el gradiente espacial se suaviza en montañas (subestima velocidades locales altas).
# El piso del gradiente (gradiente_min) es una decisión arbitraria: conviene una sensibilidad con otros valores.
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.
# Antes de correr: pon todas las capas de la serie en una carpeta (una variable puede ir en subcarpeta),
# revisa que `tiempo_desde_nombre` lea bien el tiempo de tus nombres y que el rango impreso sea sensato.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(sf)
library(terra)
library(exactextractr)

# Memoria y disco: terra decide solo cuándo escribir a disco (memfrac = fracción de RAM que puede usar) y los
# archivos temporales van a D: porque en C: se llenó el disco ("No space left on device").
dir_tmp_terra <- RUTA_TMP_TERRA; dir.create(dir_tmp_terra, showWarnings = FALSE)
terraOptions(tempdir = dir_tmp_terra, memfrac = 0.4)

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base  <- RUTA_BASE
dir_trace <- RUTA_TRACE   # salida de 00_recorte_trace_remoto.R (bio01/, bio12/), ya a ~5 km y con extensión completa
archivo_dem <- file.path(dir_base, "datos", "ambiente", "dem.tif")   # para enmascarar el océano (el DEM es NA en el mar)
dir_rej   <- file.path(dir_base, "datos", "rejilla")
dir_out   <- file.path(dir_base, "datos", "predictores")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
tamanos_km <- c(50, 100, 200)

# PARÁMETROS (propuestos por mí: confírmalos o ajústalos)
vars <- list(
  # OJO: tu bio01 de clima actual salió en décimas de °C (mediana 2.4 en lugar de 24); revisa si pasa lo mismo aquí
  temp = list(patron = "TraCE21k_bio0?1_", escala = 1, desplazamiento = 0),
  # valores mayores a relleno_mayor se tratan como dato faltante; 6500 era demasiado bajo: la precipitación real llega a ~9600 mm en el Chocó
  prec = list(patron = "TraCE21k_bio12_",  escala = 1, desplazamiento = 0, relleno_mayor = 20000)
)
# Tiempo de cada capa a partir del nombre. Con nombres tipo "CHELSA_TraCE21k_bio01_-200_V1.0.tif"
# toma -200. CONFIRMA en la documentación de CHELSA-TraCE21k la unidad de ese número
# (años por unidad) y el signo; solo importan las diferencias entre capas.
regex_tiempo <- "_(-?[0-9]+)_"
# Según la documentación técnica de CHELSA-TraCE21k (Tabla 4) cada unidad de timeID son 100 años:
# timeID 20 = 1900-1990 (0 k BP), 0 = -100 a -1 (2 k BP), -200 = -20100 a -20001 (22 k BP).
# Tus capas van de -200 a 20 cada 5 unidades (500 años). Solo importan las diferencias entre capas.
anios_por_unidad <- 100
fact_agregacion <- 1          # 00_recorte_trace_remoto.R ya agregó a ~5 km; usa 5 solo si trabajas con capas de 1 km
res_m <- 5000                 # resolución de trabajo en la proyección LAEA (m)
gradiente_min <- 0.005        # °C/km: piso del gradiente espacial para no dividir entre ~0 en zonas planas
crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# ---- 1. Plantilla de trabajo (extensión de las celdas) -----------------------------
celdas_100 <- st_read(file.path(dir_rej, "celdas_100km.gpkg"), quiet = TRUE)
bb <- st_bbox(celdas_100)
bb_m <- c(floor((unname(bb["xmin"]) - 5e4) / res_m) * res_m, ceiling((unname(bb["xmax"]) + 5e4) / res_m) * res_m,
          floor((unname(bb["ymin"]) - 5e4) / res_m) * res_m, ceiling((unname(bb["ymax"]) + 5e4) / res_m) * res_m)
plantilla <- rast(terra::ext(bb_m[1], bb_m[2], bb_m[3], bb_m[4]), resolution = res_m, crs = crs_laea)
bb_wgs <- st_bbox(st_transform(st_as_sfc(st_bbox(celdas_100)), 4326))
ext_wgs <- terra::ext(unname(bb_wgs["xmin"]) - 2, unname(bb_wgs["xmax"]) + 2,
                      unname(bb_wgs["ymin"]) - 2, unname(bb_wgs["ymax"]) + 2)
cat("Plantilla de trabajo:", ncol(plantilla), "x", nrow(plantilla), "píxeles\n")

# ---- 2. Funciones ------------------------------------------------------------
tiempo_desde_nombre <- function(f) as.numeric(str_match(basename(f), regex_tiempo)[, 2])

# Lee, recorta, agrega y proyecta toda la serie de una variable; devuelve capas ordenadas de la más antigua a la más reciente
leer_serie <- function(v, nombre) {
  cache <- file.path(dir_out, sprintf("serie_%s_laea.tif", nombre))
  fs <- list.files(dir_trace, pattern = paste0(v$patron, ".*\\.tif$"), recursive = TRUE,
                   full.names = TRUE, ignore.case = TRUE)
  if (length(fs) < 3) stop("Encontré ", length(fs), " capas para '", nombre, "' con el patrón ", v$patron, " en ", dir_trace)
  t_id <- tiempo_desde_nombre(fs)
  if (anyNA(t_id)) stop("No pude leer el tiempo de: ", paste(basename(fs[is.na(t_id)])[1:3], collapse = ", "),
                         ". Ajusta regex_tiempo.")
  o <- order(t_id); fs <- fs[o]; t_id <- t_id[o]
  anios <- t_id * anios_por_unidad
  cat(sprintf("%s: %d capas, de %s a %s; pasos entre capas: %s años\n", nombre, length(fs),
              min(anios), max(anios), paste(unique(diff(anios)), collapse = ", ")))
  # Firma de los archivos de origen (nombre, tamaño, fecha, extensión del primero). Si cambió (p. ej. re-recortaste las
  # capas TraCE a otra extensión), se borran la caché y las capas intermedias para que no se reutilice nada viejo.
  e1 <- ext(rast(fs[1]))
  firma <- c(paste(basename(fs), file.size(fs), format(file.mtime(fs)), collapse = "|"),
             paste(round(c(e1$xmin, e1$xmax, e1$ymin, e1$ymax), 3), collapse = ","),
             paste(round(c(ext_wgs$xmin, ext_wgs$xmax, ext_wgs$ymin, ext_wgs$ymax), 3), collapse = ","),
             as.character(fact_agregacion))
  f_firma <- file.path(dir_out, sprintf("firma_%s.txt", nombre))
  if (!file.exists(f_firma) || !identical(readLines(f_firma), firma)) {
    cat("  los archivos de origen o la extensión cambiaron: se borran la caché y las capas intermedias de", nombre, "\n")
    unlink(cache); unlink(file.path(dir_out, sprintf("capas_%s_laea", nombre)), recursive = TRUE)
    writeLines(firma, f_firma)
  }
  cat(sprintf("  extensión de las capas de origen: %s | extensión pedida: %s\n",
              paste(round(c(e1$xmin, e1$xmax, e1$ymin, e1$ymax), 1), collapse = ", "),
              paste(round(c(ext_wgs$xmin, ext_wgs$xmax, ext_wgs$ymin, ext_wgs$ymax), 1), collapse = ", ")))
  if (file.exists(cache)) {
    cat("  usando la caché", basename(cache), "\n")
    r <- rast(cache)
  } else {
    # Una capa a la vez y guardando cada resultado en disco: las capas originales (~1 km, ~90 millones de píxeles)
    # no caben en memoria si se acumulan. Si R se cierra, al volver a correr retoma desde las capas ya hechas.
    dir_tmp <- file.path(dir_out, sprintf("capas_%s_laea", nombre)); dir.create(dir_tmp, showWarnings = FALSE)
    capas <- character(length(fs))
    for (k in seq_along(fs)) {
      capas[k] <- file.path(dir_tmp, sprintf("t%s.tif", t_id[k]))
      if (file.exists(capas[k])) next
      x <- crop(rast(fs[k])[[1]], ext_wgs)
      if (!is.null(v$relleno_mayor)) x <- ifel(x > v$relleno_mayor, NA, x)
      if (fact_agregacion > 1) x <- aggregate(x, fact = fact_agregacion, fun = "mean", na.rm = TRUE)
      x <- project(x, plantilla, method = "bilinear") * v$escala + v$desplazamiento
      writeRaster(x, capas[k], overwrite = TRUE)
      rm(x); gc(verbose = FALSE)
      cat("  capa", k, "de", length(fs), "lista\n")
    }
    r <- rast(capas)
    writeRaster(r, cache, overwrite = TRUE)
    r <- rast(cache)
  }
  names(r) <- paste0("t", anios)
  # Control de calidad por capa: min, media, máx y % de NA (revisa las capas vacías o constantes)
  qc <- as_tibble(global(r, c("min", "mean", "max"), na.rm = TRUE)) %>%
    mutate(anios = anios, archivo = basename(fs), pct_na = round(100 * as.numeric(unlist(global(is.na(r), "mean"))), 2))
  write_csv(qc, file.path(dir_out, sprintf("qc_serie_%s.csv", nombre)))
  malas <- which(is.na(qc$mean) | qc$max == qc$min | qc$pct_na > 50)
  if (length(malas) > 0) {
    warning("Capas sospechosas en '", nombre, "' (vacías, constantes o con >50 % de NA): ",
            paste(qc$archivo[malas], collapse = ", "))
  }
  list(r = r, anios = anios)
}

# Gradiente espacial (°C/km) con el filtro de Horn sobre un ráster proyectado en metros
gradiente_espacial <- function(r) {
  kx <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), 3, 3, byrow = TRUE) / (8 * res_m)
  ky <- matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), 3, 3, byrow = TRUE) / (8 * res_m)
  gx <- focal(r, w = kx, fun = "sum"); gy <- focal(r, w = ky, fun = "sum")
  sqrt(gx^2 + gy^2) * 1000      # °C por m -> °C por km
}

# ---- 3. Temperatura: anomalía y velocidad ----------------------------------------------
st <- leer_serie(vars$temp, "temp")
tr <- st$r; ta <- st$anios; n <- nlyr(tr)
cat("Temperatura, capa más antigua y más reciente:\n"); print(global(tr[[c(1, n)]], c("min", "mean", "max"), na.rm = TRUE))
cat("  (revisa que los valores sean °C razonables para el Neotrópico; si no, ajusta escala/desplazamiento)\n")
if (max(abs(unlist(global(tr[[n]], c("min", "max"), na.rm = TRUE)))) < 6) {
  warning("La temperatura de la capa más reciente no llega a 6: parece estar en décimas de °C. Usa escala = 10.")
}

anom_temp <- abs(tr[[n]] - tr[[1]])
# Nota sobre el tiempo: el id de archivo de CHELSA-TraCE21k va de -200 a 20 en pasos de 100 años. Según el
# paquete pastclim, id 20 = presente (1950) e id -200 = 22 000 años antes; la columna "anios" de este script
# es id*100 (de -20000 a 2000) y está desplazada 2000 años respecto a esa escala. Como solo se usan
# diferencias de tiempo (22 000 años entre la primera y la última capa), la velocidad no cambia.
names(anom_temp) <- "hist_anom_temp"
sd_temp <- app(tr, "sd", na.rm = TRUE); names(sd_temp) <- "hist_sd_temp"

# velocidad neta: tasa media / gradiente espacial de la capa presente
anios_tot <- ta[n] - ta[1]
g_pres <- gradiente_espacial(tr[[n]])
vel_neta <- (anom_temp / anios_tot) / ifel(g_pres < gradiente_min, gradiente_min, g_pres) * 1000   # km/año -> m/año
names(vel_neta) <- "hist_vel_neta_log10"
vel_neta <- log10(vel_neta + 1e-12)

# velocidad por intervalos y su media
vint <- vector("list", n - 1)
for (i in seq_len(n - 1)) {
  dT <- abs(tr[[i + 1]] - tr[[i]]) / (ta[i + 1] - ta[i])
  g <- gradiente_espacial((tr[[i]] + tr[[i + 1]]) / 2)
  vint[[i]] <- log10(dT / ifel(g < gradiente_min, gradiente_min, g) * 1000 + 1e-12)
}
vel_media <- app(rast(vint), "mean", na.rm = TRUE); names(vel_media) <- "hist_vel_media_log10"

# ---- 4. Precipitación: anomalía --------------------------------------------------------
sp <- leer_serie(vars$prec, "prec")
pr <- sp$r; npr <- nlyr(pr)
cat("Precipitación, capa más antigua y más reciente:\n"); print(global(pr[[c(1, npr)]], c("min", "mean", "max"), na.rm = TRUE))
anom_prec <- abs(log((pr[[npr]] + 1) / (pr[[1]] + 1))); names(anom_prec) <- "hist_anom_prec"

# variabilidad temporal de la precipitación: absoluta (mm, como el borrador 6) y relativa (log)
pr_pos <- clamp(pr, lower = 0, values = TRUE)
sd_prec_mm  <- app(pr_pos, "sd", na.rm = TRUE);          names(sd_prec_mm)  <- "hist_sd_prec_mm"
sd_prec_log <- app(log(pr_pos + 1), "sd", na.rm = TRUE); names(sd_prec_log) <- "hist_sd_prec_log"

metricas <- c(anom_temp, anom_prec, vel_neta, vel_media, sd_temp, sd_prec_mm, sd_prec_log)

# Solo tierra actual: las capas de CHELSA traen valores también sobre el océano, y en las celdas
# costeras esos píxeles contaminarían el promedio. El DEM es NA en el mar.
if (file.exists(archivo_dem)) {
  dem_c <- crop(rast(archivo_dem)[[1]], ext_wgs)
  # tierra = DEM con valor y > 0 (si el DEM trae 0 o batimetría negativa en el mar, esos píxeles cuentan como mar)
  es_mar <- is.na(dem_c) | dem_c <= 0
  tierra <- project(ifel(es_mar, 0, 1), plantilla, method = "average")
  metricas <- mask(metricas, ifel(tierra >= 0.5, 1, NA))
  cat("Enmascarado a tierra con el DEM (fracción de píxeles de mar según el DEM:",
      round(as.numeric(unlist(global(es_mar, "mean"))), 2), ")\n")
} else message("Sin DEM: no se enmascaró el océano.")
writeRaster(metricas, file.path(dir_out, "historia_metricas_laea.tif"), overwrite = TRUE)

# ---- 5. Resumen por celda -------------------------------------------------------------
for (tam_km in tamanos_km) {
  celdas <- st_read(file.path(dir_rej, sprintf("celdas_%dkm.gpkg", tam_km)), quiet = TRUE) %>% filter(mantener)
  z <- exact_extract(metricas, celdas, fun = "mean", progress = FALSE)
  names(z) <- str_remove(names(z), "^mean\\.")
  out <- bind_cols(tibble(celda = celdas$celda), z)
  write_csv(out, file.path(dir_out, sprintf("historia_%dkm.csv", tam_km)))
  cat(sprintf("\nCeldas de %d km: %d; NA por variable:\n", tam_km, nrow(out))); print(colSums(is.na(out[, -1])))
}

# Revisión recomendada: mapa de cada métrica para comprobar que el patrón tiene sentido
# (p. ej. anomalía mayor en latitudes altas y zonas montañosas; velocidad alta en zonas planas).
# plot(metricas)
# Después corre 05_predictores_celdas.R: si encuentra historia_{tam}km.csv, lo usa en lugar de la anomalía simple.
