# 05c_topo_rango.R
# Rango de elevación de cada celda = (máximo de la capa de elevación máxima) - (mínimo de la capa de elevación mínima) dentro de la celda,
# más el rango entre los percentiles 5 y 95 de la capa mediana (menos sensible a un pixel extremo). Medida alterna de heterogeneidad
# topográfica para la robustez (script 25, variante topo_rango).
# Capas EarthEnv (GMTED2010, 30 arc-seg, ~1 km, mundiales). Usa las MISMAS que el script 05 para que todo salga de la misma fuente:
#   ruta_max : elevación máxima  (tu dem.tif: "Maximum Statistic")
#   ruta_min : elevación mínima  (EarthEnv "Minimum Statistic"; si no la tienes, pon NULL y se usa la mediana como piso, lo que SUBESTIMA algo el rango)
#   ruta_med : elevación mediana (la que usaste para topo_elev_media; sirve para los percentiles)
# Como los rasters son mundiales (43200 x 16800), primero se recortan a la zona de las celdas y solo entonces se extrae (evita agotar memoria).
# Las celdas se reconstruyen desde la base (centro y tamaño) en la azimutal equiárea de Lambert del análisis (lon_0 = -75, lat_0 = 0).
# Salida: topo_rango_{tam}km.csv (celda, tam_km, topo_rango_elev, topo_rango_p5_p95, n_pixeles). El rango crece con el tamaño de celda:
# se compara dentro de cada tamaño. No probado (no hay R en mi entorno). Requiere: terra, dplyr, readr.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(terra); library(dplyr); library(readr)

dir_base <- RUTA_BASE
dir_mod  <- file.path(dir_base, "datos", "modelos")
ruta_max <- file.path(dir_amb, "elevation_1KMma_GMTEDma.tif")             # elevación MÁXIMA
ruta_min <- file.path(dir_amb, "elevation_1KMmi_GMTEDmi.tif")             # elevación MÍNIMA
ruta_med <- file.path(dir_amb, "elevation_1KMmd_GMTEDmd.tif")             # elevación MEDIANA
tamanos_km <- c(100, 50, 200)
crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
terraOptions(memfrac = 0.5)

poligonos <- function(tam_km, crs_salida) {
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE) %>% distinct(celda, x_centro_km, y_centro_km)
  h <- tam_km * 1000 / 2
  xy <- do.call(rbind, lapply(seq_len(nrow(base)), function(i) {
    cx <- base$x_centro_km[i] * 1000; cy <- base$y_centro_km[i] * 1000
    cbind(object = i, part = 1, x = c(cx - h, cx + h, cx + h, cx - h, cx - h), y = c(cy - h, cy - h, cy + h, cy + h, cy - h), hole = 0)
  }))
  v <- vect(xy, type = "polygons", crs = crs_laea)
  list(base = base, v = project(densify(v, interval = 10000), crs_salida))
}

r_max <- rast(ruta_max)[[1]]
cat(sprintf("Capa máxima: %s | %d x %d pixeles | resolución ~%.0f m | CRS: %s\n", basename(ruta_max), ncol(r_max), nrow(r_max), res(r_max)[1] * 111320, crs(r_max, describe = TRUE)$name))
# recorte a la zona de estudio (celdas de 200 km, que cubren todo) con 1 grado de margen
e <- ext(poligonos(200, crs(r_max))$v) + 1
r_max <- crop(r_max, e)
r_min <- if (!is.null(ruta_min)) crop(rast(ruta_min)[[1]], e) else NULL
r_med <- if (!is.null(ruta_med)) crop(rast(ruta_med)[[1]], e) else NULL
piso <- if (!is.null(r_min)) r_min else r_med
if (is.null(piso)) stop("Necesito al menos la capa mínima o la mediana para calcular el rango.")
if (is.null(r_min)) cat("AVISO: sin capa mínima; el rango se calcula contra la mediana (subestima el rango real).\n")
mu <- spatSample(r_max, 5e4, method = "regular", na.rm = TRUE)[[1]]
cat(sprintf("Muestra de la capa máxima -> min: %.0f | max: %.0f | media: %.0f (elevación: min cerca de 0, max de miles)\n", min(mu), max(mu), mean(mu)))

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  p <- poligonos(tam_km, crs(r_max)); base <- p$base; v <- p$v
  emax <- extract(r_max, v, fun = "max", na.rm = TRUE, ID = FALSE)[[1]]
  emin <- extract(piso,  v, fun = "min", na.rm = TRUE, ID = FALSE)[[1]]
  out <- tibble(ID = seq_len(nrow(base)), celda = base$celda, tam_km = tam_km,
                topo_rango_elev = ifelse(is.finite(emax) & is.finite(emin), emax - emin, NA_real_))
  if (!is.null(r_med)) {
    px <- extract(r_med, v, ID = TRUE); names(px) <- c("ID", "z")
    res <- px %>% filter(!is.na(z)) %>% group_by(ID) %>%
      summarise(topo_rango_p5_p95 = diff(quantile(z, c(0.05, 0.95), names = FALSE)), n_pixeles = n(), .groups = "drop")
    out <- out %>% left_join(res, by = "ID") %>% mutate(n_pixeles = ifelse(is.na(n_pixeles), 0L, n_pixeles))
  } else out <- out %>% mutate(topo_rango_p5_p95 = NA_real_, n_pixeles = NA_integer_)
  out <- select(out, celda, tam_km, topo_rango_elev, topo_rango_p5_p95, n_pixeles)
  cat(sprintf("  celdas: %d | sin dato: %d | rango medio: %.0f m | mediana del rango: %.0f m | p5-p95 medio: %.0f m\n", nrow(out),
              sum(is.na(out$topo_rango_elev)), mean(out$topo_rango_elev, na.rm = TRUE), median(out$topo_rango_elev, na.rm = TRUE),
              mean(out$topo_rango_p5_p95, na.rm = TRUE)))
  write_csv(out, file.path(dir_mod, sprintf("topo_rango_%dkm.csv", tam_km)))
}
