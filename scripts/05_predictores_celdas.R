# 05_predictores_celdas.R
# Predictores ambientales por celda (50, 100 y 200 km) a partir de las celdas del script 04.
#
# Bloques (según el plan):
#   H1 Topografía : SD y rango de elevación (DEM), rugosidad (TRI) opcional
#   H2 Clima actual: SD de bio01 y bio12 (CHELSA), más medias como controles
#   H3 Historia   : anomalía de temperatura y precipitación desde el LGM (CHELSA-TraCE21k)
#   Hábitat       : SD y media del EVI, número y diversidad de clases de cobertura, cobertura forestal
#   Controles     : esfuerzo de muestreo (log de especímenes de los géneros focales) y área de tierra
#
# Cada bloque se calcula solo si encuentra su capa; si falta, avisa y sigue con los demás.
# Las estadísticas por celda se ponderan por la fracción de cada píxel dentro de la celda
# (paquete exactextractr). Las celdas se transforman al CRS de cada capa.
#
# Capas que debes colocar en datos/ambiente (nombres configurables abajo). Todas en lon-lat (EPSG:4326):
#   dem.tif                                  ya la tienes (EarthEnv GMTED2010, mediana a 1 km)
#   chelsa/bio01_actual.tif, bio12_actual.tif   CHELSA v2.1, 1981-2010
#   trace/bio01_lgm.tif, bio12_lgm.tif          CHELSA-TraCE21k en el LGM. Confirma en su documentación
#                                               la versión y el identificador de tiempo del LGM
#   modis/evi_anual.tif                         EVI medio anual (tú decides producto y años)
#   cobertura/cobertura.tif                     clases de cobertura (p. ej. ESA WorldCover o Copernicus, ~100 m)
#   bosque/cobertura_forestal.tif               % de cobertura arbórea (p. ej. Hansen, agregada a <= 1 km)
# Registra versión y fecha de descarga de cada capa en el documento de planeación.
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.
# Paquetes nuevos: install.packages("exactextractr")

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(sf)
library(terra)
library(exactextractr)
select <- dplyr::select

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir_amb  <- file.path(dir_base, "datos", "ambiente")
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_out  <- file.path(dir_base, "datos", "predictores")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

tamanos_km <- c(50, 100, 200)

dir_clima_actual <- RUTA_CLIMA_ACTUAL   # CHELSA v2.1 recortada al Neotrópico
archivos <- list(
  dem        = file.path(dir_amb, "dem.tif"),
  bio01      = file.path(dir_clima_actual, "CHELSA_bio01_neotropics.tif"),
  bio12      = file.path(dir_clima_actual, "CHELSA_bio12_neotropics.tif"),
  bio01_lgm  = file.path(dir_amb, "trace", "bio01_lgm.tif"),
  bio12_lgm  = file.path(dir_amb, "trace", "bio12_lgm.tif"),
  evi        = file.path(dir_amb, "modis", "evi_anual.tif"),
  cobertura  = file.path(dir_amb, "cobertura", "cobertura.tif"),
  bosque     = file.path(dir_amb, "bosque", "cobertura_forestal.tif")
)

# PARÁMETROS (propuestos por mí: confírmalos o ajústalos)
escala_bio01 <- 10                # tu bio01 actual viene en décimas de °C (mediana 2.4): x10 la lleva a °C
escala_bio12 <- 1                 # factor para llevar bio12 a mm
escala_evi   <- 1                 # factor para llevar el EVI a su escala natural (-1 a 1)
umbral_clase_cobertura <- 0.01    # una clase cuenta si cubre al menos esta fracción de la celda
clases_excluir_cobertura <- integer(0)   # códigos a excluir (p. ej. agua), según tu capa
calcular_tri <- TRUE              # rugosidad (TRI) a partir del DEM
# Hábitat: por defecto solo la media del EVI como control de productividad. La heterogeneidad de hábitat
# (SD del EVI, clases de cobertura, cobertura forestal) refleja el paisaje actual y es más reciente que los
# patrones de diversidad; se calcula solo como exploración o sensibilidad si lo activas.
calcular_habitat_exploratorio <- FALSE
n_min_registros_exploracion <- 10 # solo para la correlación y el VIF exploratorios
crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

cargar <- function(ruta, nombre) {
  if (!file.exists(ruta)) { message("Falta la capa '", nombre, "': ", ruta); return(NULL) }
  rast(ruta)
}
capas <- Map(cargar, archivos, names(archivos))

# ---- 1. Funciones ------------------------------------------------------------
recortar <- function(r, pol_r) {
  bb <- st_bbox(pol_r)
  e <- terra::ext(unname(bb["xmin"]) - 0.5, unname(bb["xmax"]) + 0.5,
                  unname(bb["ymin"]) - 0.5, unname(bb["ymax"]) + 0.5)
  terra::crop(r, e)
}

# Estadísticas por celda ponderadas por la fracción de píxel; devuelve un data.frame con columnas prefijo_fun
zonal <- function(r, pol, funs, prefijo) {
  res <- exact_extract(r, pol, fun = funs, progress = FALSE)
  res <- as.data.frame(res)
  if (length(funs) == 1) names(res) <- funs
  names(res) <- paste0(prefijo, "_", names(res))
  res
}

# Deja solo la tierra actual: las capas de CHELSA traen valores también sobre el océano (el DEM es NA en el mar)
a_tierra <- function(r, dem_c) {
  d <- terra::resample(dem_c, r, method = "near")
  terra::mask(r, terra::ifel(is.na(d), NA, 1))
}

resumen_capa <- function(r, nombre) {
  g <- global(r, c("min", "max"), na.rm = TRUE)
  cat(sprintf("  %s: min %.2f, max %.2f\n", nombre, g$min, g$max))
}

id_celda <- function(ix, iy, tam_km) sprintf("c%d_%d_%d", as.integer(tam_km), as.integer(ix), as.integer(iy))

vif <- function(df) {
  sapply(names(df), function(v) {
    r2 <- summary(lm(reformulate(setdiff(names(df), v), response = v), data = df))$r.squared
    1 / (1 - r2)
  })
}

# ---- 2. Esfuerzo de muestreo (especímenes de los géneros focales) -----------------
# Se usan todos los especímenes con coordenadas de los géneros focales, incluidos los que no
# llegaron a especie, como medida de esfuerzo de colecta. Cambia el archivo si prefieres la familia completa.
esf <- readRDS(file.path(dir_gbif, "generos_focales_bruto.rds"))
if ("basisOfRecord" %in% names(esf)) esf <- esf %>% filter(basisOfRecord == "PRESERVED_SPECIMEN")
esf <- esf %>% filter(!is.na(decimalLatitude), !is.na(decimalLongitude))
cat("\nEspecímenes para el esfuerzo:", nrow(esf), "\n")
if ("year" %in% names(esf)) {
  cat("Años de colecta (cuantiles) y proporción de especímenes con año:\n")
  print(quantile(esf$year, c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE))
  cat("Sin año:", round(100 * mean(is.na(esf$year)), 1), "% | desde 2000:",
      round(100 * mean(esf$year >= 2000, na.rm = TRUE), 1), "% | antes de 1980:",
      round(100 * mean(esf$year < 1980, na.rm = TRUE), 1), "%\n")
}
esf_xy <- st_as_sf(esf, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs_laea) %>% st_coordinates()

# ---- 3. Bucle por tamaño de celda ----------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  celdas <- st_read(file.path(dir_rej, sprintf("celdas_%dkm.gpkg", tam_km)), quiet = TRUE) %>%
    filter(mantener)
  cat("Celdas conservadas:", nrow(celdas), "\n")
  base <- st_drop_geometry(celdas) %>%
    select(any_of(c("celda", "tam_km", "celda_100", "celda_200", "lon_centro", "lat_centro",
                    "area_tierra_km2", "frac_tierra", "motivo", "n_reg", "S_obs", "cob_good")))

  pred <- base %>% select(celda)

  # Celdas en el CRS de cada capa (todas en lon-lat; se toma el de la primera disponible)
  cap_ref <- Filter(Negate(is.null), capas)[[1]]
  pol <- st_transform(celdas[, "celda"], st_crs(terra::crs(cap_ref)))
  bbp <- st_bbox(pol)
  for (nm in names(capas)) {
    if (is.null(capas[[nm]])) next
    r_ <- capas[[nm]]
    if (unname(bbp["xmin"]) < terra::xmin(r_) || unname(bbp["xmax"]) > terra::xmax(r_) ||
        unname(bbp["ymin"]) < terra::ymin(r_) || unname(bbp["ymax"]) > terra::ymax(r_)) {
      message("AVISO: hay celdas fuera de la extensión de la capa '", nm, "'; esas celdas saldrán con NA parcial o total")
    }
  }

  # H1 Topografía
  if (!is.null(capas$dem)) {
    cat("H1 Topografía\n")
    dem <- recortar(capas$dem, pol)
    z <- zonal(dem, pol, c("mean", "stdev", "min", "max"), "topo_elev")
    pred <- bind_cols(pred, tibble(topo_elev_media = z$topo_elev_mean,
                                   topo_sd_elev = z$topo_elev_stdev,
                                   topo_rango_elev = z$topo_elev_max - z$topo_elev_min))
    if (calcular_tri) {
      tri <- terra::terrain(dem, v = "TRI")
      pred$topo_tri_media <- zonal(tri, pol, "mean", "tri")$tri_mean
    }
  }

  # H2 Clima actual
  b1 <- b12 <- NULL
  if (!is.null(capas$bio01)) {
    cat("H2 Clima actual\n")
    b1 <- recortar(capas$bio01, pol) * escala_bio01
    if (!is.null(capas$dem)) b1 <- a_tierra(b1, recortar(capas$dem, pol))
    resumen_capa(b1, "bio01 (°C esperados)")
    z <- zonal(b1, pol, c("mean", "stdev"), "bio01")
    pred$clima_bio01_media <- z$bio01_mean; pred$clima_bio01_sd <- z$bio01_stdev
  }
  if (!is.null(capas$bio12)) {
    b12 <- recortar(capas$bio12, pol) * escala_bio12
    # CHELSA v2.1 guarda bio12 como entero de 16 bits en décimas de mm: el máximo posible es 6553.5, y ese valor
    # significa "6553.5 o más" (saturado), NO dato faltante. Se comprobó: 9 157 píxeles reales entre 6500 y 6553 y
    # 252 236 píxeles saturados en 6553.5 (Chocó y Pacífico oriental). Se conservan; en el Chocó y el Caribe de
    # Costa Rica y Nicaragua la lluvia real es mayor, así que la media y la SD de esas celdas están subestimadas.
    if (!is.null(capas$dem)) b12 <- a_tierra(b12, recortar(capas$dem, pol))
    resumen_capa(b12, "bio12 (mm esperados)")
    z <- zonal(b12, pol, c("mean", "stdev"), "bio12")
    pred$clima_bio12_media <- z$bio12_mean; pred$clima_bio12_sd <- z$bio12_stdev
  }
  # SD de bio01 libre de la SD de elevación (residual), calculada con las celdas disponibles
  if (all(c("clima_bio01_sd", "topo_sd_elev") %in% names(pred))) {
    ok <- complete.cases(pred$clima_bio01_sd, pred$topo_sd_elev)
    ajuste <- lm(clima_bio01_sd ~ topo_sd_elev, data = pred[ok, ])
    pred$clima_bio01_sd_res <- NA_real_
    pred$clima_bio01_sd_res[ok] <- resid(ajuste)
    cat("  R2 de SD bio01 ~ SD elevación:", round(summary(ajuste)$r.squared, 2), "\n")
  }

  # H3 Historia: si existe historia_{tam}km.csv (script 05b) se usa; si no, anomalía simple con las capas del LGM
  arch_hist <- file.path(dir_out, sprintf("historia_%dkm.csv", tam_km))
  usar_hist <- file.exists(arch_hist)
  if (usar_hist) {
    cat("H3 Historia: usando", basename(arch_hist), "\n")
    pred <- left_join(pred, read_csv(arch_hist, show_col_types = FALSE), by = "celda")
  }
  if (!usar_hist && !is.null(b1) && !is.null(capas$bio01_lgm)) {
    cat("H3 Historia (temperatura)\n")
    l1 <- recortar(capas$bio01_lgm, pol) * escala_bio01
    l1 <- terra::resample(l1, b1)
    resumen_capa(l1, "bio01 LGM")
    pred$hist_anom_temp <- zonal(abs(b1 - l1), pol, "mean", "anom")$anom_mean
  }
  if (!usar_hist && !is.null(b12) && !is.null(capas$bio12_lgm)) {
    cat("H3 Historia (precipitación)\n")
    l12 <- recortar(capas$bio12_lgm, pol) * escala_bio12
    l12 <- terra::resample(l12, b12)
    resumen_capa(l12, "bio12 LGM")
    # diferencia absoluta del logaritmo de la precipitación (cambio relativo)
    dif <- abs(log((b12 + 1) / (l12 + 1)))
    pred$hist_anom_prec <- zonal(dif, pol, "mean", "anomp")$anomp_mean
  }

  # Hábitat
  if (!is.null(capas$evi)) {
    cat("Hábitat: EVI\n")
    evi <- recortar(capas$evi, pol) * escala_evi
    resumen_capa(evi, "EVI")
    z <- zonal(evi, pol, c("mean", "stdev"), "evi")
    pred$hab_evi_media <- z$evi_mean
    if (calcular_habitat_exploratorio) pred$hab_evi_sd <- z$evi_stdev
  }
  if (calcular_habitat_exploratorio && !is.null(capas$cobertura)) {
    cat("Hábitat: cobertura\n")
    cob <- recortar(capas$cobertura, pol)
    fr <- zonal(cob, pol, "frac", "cob")
    fr <- fr[, !names(fr) %in% paste0("cob_frac_", clases_excluir_cobertura), drop = FALSE]
    p <- as.matrix(fr); p[is.na(p)] <- 0
    tot <- rowSums(p); tot[tot == 0] <- NA
    p <- p / tot
    pred$hab_n_clases <- rowSums(p >= umbral_clase_cobertura)
    pred$hab_shannon <- -rowSums(ifelse(p > 0, p * log(p), 0))
  }
  if (calcular_habitat_exploratorio && !is.null(capas$bosque)) {
    cat("Hábitat: cobertura forestal\n")
    pred$hab_bosque_pct <- zonal(recortar(capas$bosque, pol), pol, "mean", "bosque")$bosque_mean
  }

  # Controles
  ind <- list(ix = floor(esf_xy[, 1] / (tam_km * 1000)), iy = floor(esf_xy[, 2] / (tam_km * 1000)))
  n_esf <- table(id_celda(ind$ix, ind$iy, tam_km))
  pred$ctrl_n_esfuerzo <- as.integer(n_esf[pred$celda]); pred$ctrl_n_esfuerzo[is.na(pred$ctrl_n_esfuerzo)] <- 0L
  pred$ctrl_esfuerzo_log <- log10(pred$ctrl_n_esfuerzo + 1)

  salida <- base %>% left_join(pred, by = "celda") %>%
    mutate(ctrl_area_tierra_km2 = area_tierra_km2)
  write_csv(salida, file.path(dir_out, sprintf("predictores_%dkm.csv", tam_km)))

  # Resumen y exploración de colinealidad (solo en celdas con registros suficientes)
  vars <- names(salida)[str_detect(names(salida), "^(topo|clima|hist|hab|ctrl)_")]
  res <- salida %>% summarise(across(all_of(vars), list(
    n_NA = ~sum(is.na(.x)), media = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE),
    min = ~min(.x, na.rm = TRUE), max = ~max(.x, na.rm = TRUE)), .names = "{.col}__{.fn}")) %>%
    pivot_longer(everything(), names_to = c("variable", "estadistico"), names_sep = "__") %>%
    pivot_wider(names_from = estadistico, values_from = value)
  write_csv(res, file.path(dir_out, sprintf("resumen_predictores_%dkm.csv", tam_km)))
  print(as.data.frame(res), digits = 3)

  expl <- salida %>% filter(n_reg >= n_min_registros_exploracion) %>%
    select(all_of(setdiff(vars, c("ctrl_n_esfuerzo")))) %>% select(where(~ !all(is.na(.x)))) %>% drop_na()
  cat("Celdas para la exploración (n_reg >=", n_min_registros_exploracion, "y sin NA):", nrow(expl), "\n")
  if (nrow(expl) > 30 && ncol(expl) > 2) {
    write_csv(as_tibble(cor(expl), rownames = "variable"), file.path(dir_out, sprintf("correlacion_%dkm.csv", tam_km)))
    v <- vif(expl)
    write_csv(tibble(variable = names(v), vif = round(v, 2)) %>% arrange(desc(vif)),
              file.path(dir_out, sprintf("vif_%dkm.csv", tam_km)))
    cat("VIF más altos:\n"); print(sort(round(v, 1), decreasing = TRUE)[seq_len(min(6, length(v)))])
  }
}

# Los predictores se estandarizan (media 0, SD 1) y se transforman (log) ya en el script de modelos,
# sobre el conjunto final de celdas, una vez fijado el criterio de inclusión de celdas.
