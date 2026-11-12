# 20_topo_rugosidad_fina.R
# Paso 1 de la idea "desviación estándar de la pendiente": ¿aporta información distinta de topo_sd_elev y se separa de la
# temperatura? Solo calcula y compara; NO modifica los modelos.
#
# Capas de Amatulli et al. (2018, Sci. Data 5:180040), GMTED2010, a 1 km. Nombre = variable_agregación_capa fuente
# (Tabla 1 del artículo); cada píxel de 1 km agrega 4 x 4 píxeles de 250 m, sobre variables derivadas con ventana de 3 x 3:
#   slope_1KMsd_GMTEDmd       SD de la pendiente (derivada de la mediana de elevación) entre los 16 píxeles de 250 m de cada
#                             píxel de 1 km. Es la "desviación estándar de la pendiente" de la idea original.
#   tri_1KMmd_GMTEDmd         mediana del TRI (diferencia absoluta media de elevación con los 8 vecinos de 250 m).
#   vrm_1KMmd_GMTEDmd         (opcional) mediana de la Vector Ruggedness Measure (0 = plano, 1 = muy rugoso): variación local de la
#                             pendiente y la orientación, más independiente de la pendiente media que el TRI (Amatulli et al., pág. 6).
#   elevation_1KMsd_GMTEDsd   OJO: agrega con "sd" la capa fuente GMTEDsd (que ya es la SD de la elevación nativa dentro de cada
#                             píxel de 250 m), así que es la variación entre píxeles de 250 m de esas SD, NO la SD de elevación
#                             dentro del kilómetro. La SD "combinada" correcta sería elevation_1KMpsd_GMTEDsd (o 1KMsd_GMTEDmd,
#                             la SD de las medianas de 250 m). Se incluye solo como comparación; interprétala con cautela.
# Por celda se toma la MEDIA ponderada de esos valores de 1 km (exactextractr): rugosidad "fina" (escala sub-kilométrica),
# distinta de topo_sd_elev, que es la SD de la elevación ENTRE píxeles de 1 km dentro de la celda (relieve a escala de celda).
#
# Salidas (datos/modelos): topo_fina_{tam}km.csv (valores por celda), topo_fina_correlaciones_{tam}km.csv (Pearson y
#   Spearman entre predictores topográficos y climáticos), topo_fina_resumen_{tam}km.csv (por variable nueva: r con
#   topo_sd_elev, r con la SD de temperatura, R2, VIF dentro del conjunto de 10 predictores).
# No probado: no tengo R en mi entorno. Requiere: sf, terra, exactextractr, dplyr, readr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(sf)
library(terra)
library(exactextractr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_pred <- file.path(dir_base, "datos", "predictores")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_topo <- RUTA_TOPO_AMATULLI      # capas de Amatulli

todas_las_escalas <- FALSE   # TRUE = corre 100, 50 y 200 km (en ese orden); FALSE = solo 100 km
tamanos_km <- if (todas_las_escalas) c(100, 50, 200) else c(100)

# nombre en la columna = inicio del nombre del archivo (la extensión se busca sola)
capas_topo <- c(topo_fina_sd_elev  = "elevation_1KMsd_GMTEDsd",
                topo_fina_sd_pend  = "slope_1KMsd_GMTEDmd",
                topo_fina_tri_med  = "tri_1KMmd_GMTEDmd",
                topo_fina_vrm_med  = "vrm_1KMmd_GMTEDmd")   # opcional: se omite si no la has descargado
referencia <- c("topo_sd_elev", "topo_tri_media", "topo_elev_media", "topo_rango_elev",
                "clima_bio01_sd", "clima_bio01_sd_res", "clima_bio12_sd", "hist_sd_temp", "hist_anom_prec")
conjunto_10 <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec",
                 "clima_bio01_media", "clima_bio12_media", "hab_evi_media", "ctrl_esfuerzo_log", "ctrl_area_tierra_km2")

# ---- 1. Lectura de capas --------------------------------------------------------------------
buscar <- function(prefijo) {
  f <- list.files(dir_topo, pattern = paste0("^", prefijo, "\\.(tif|tiff|vrt|asc|img)$"), full.names = TRUE, ignore.case = TRUE)
  if (!length(f)) { message("No encuentro '", prefijo, "' en ", dir_topo, "; se omite esa capa"); return(NULL) }
  f[1]
}
rutas <- lapply(capas_topo, buscar)
capas_topo <- capas_topo[!vapply(rutas, is.null, logical(1))]
if (!length(capas_topo)) stop("No se encontró ninguna capa en ", dir_topo, ". Archivos presentes:\n", paste(list.files(dir_topo), collapse = "\n"))
capas <- lapply(rutas[names(capas_topo)], rast)
cat("Capas leídas:\n"); for (n in names(capas)) { g <- global(capas[[n]], c("min", "max"), na.rm = TRUE)
  cat(sprintf("  %s: CRS %s; rango %.2f a %.2f\n", n, crs(capas[[n]], describe = TRUE)$code, g$min, g$max)) }

vif <- function(df) sapply(names(df), function(v) 1 / (1 - summary(lm(reformulate(setdiff(names(df), v), response = v), data = df))$r.squared))

# ---- 2. Por tamaño de celda --------------------------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  celdas <- st_read(file.path(dir_rej, sprintf("celdas_%dkm.gpkg", tam_km)), quiet = TRUE) %>% filter(mantener)
  pred   <- read_csv(file.path(dir_pred, sprintf("predictores_%dkm.csv", tam_km)), show_col_types = FALSE)
  base   <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  usadas <- unique(base$celda)

  nuevos <- tibble(celda = celdas$celda)
  for (n in names(capas)) {
    r <- capas[[n]]
    pol <- st_transform(celdas[, "celda"], st_crs(crs(r)))
    bb <- st_bbox(pol)
    r <- crop(r, ext(unname(bb["xmin"]) - 0.5, unname(bb["xmax"]) + 0.5, unname(bb["ymin"]) - 0.5, unname(bb["ymax"]) + 0.5))
    nuevos[[n]] <- exact_extract(r, pol, fun = "mean", progress = FALSE)
  }
  d <- pred %>% inner_join(nuevos, by = "celda") %>% filter(celda %in% usadas)
  cat("Celdas usadas en modelos:", nrow(d), "; con NA en alguna capa nueva:", sum(!complete.cases(d[names(capas_topo)])), "\n")
  write_csv(nuevos, file.path(dir_mod, sprintf("topo_fina_%dkm.csv", tam_km)))

  # correlaciones
  vars <- c(intersect(referencia, names(d)), names(capas_topo))
  dd <- d[vars] %>% filter(complete.cases(.))
  cor_long <- function(metodo) {
    m <- cor(dd, method = metodo); as.data.frame(as.table(m)) %>% setNames(c("var1", "var2", "r")) %>% mutate(metodo = metodo)
  }
  cors <- bind_rows(cor_long("pearson"), cor_long("spearman")) %>% mutate(tam_km = tam_km)
  write_csv(cors, file.path(dir_mod, sprintf("topo_fina_correlaciones_%dkm.csv", tam_km)))

  # resumen por variable nueva (y por topo_sd_elev como referencia)
  base_10 <- intersect(conjunto_10, names(d))
  res <- bind_rows(lapply(c("topo_sd_elev", names(capas_topo)), function(v) {
    x <- d[[v]]
    out <- tibble(variable = v, n = sum(!is.na(x)),
                  r_topo_sd_elev = if (v == "topo_sd_elev") NA_real_ else cor(x, d$topo_sd_elev, use = "complete.obs"),
                  r_topo_tri_media = if ("topo_tri_media" %in% names(d)) cor(x, d$topo_tri_media, use = "complete.obs") else NA_real_,
                  r_temp_sd = if ("clima_bio01_sd" %in% names(d)) cor(x, d$clima_bio01_sd, use = "complete.obs") else NA_real_,
                  r_temp_sd_res = if ("clima_bio01_sd_res" %in% names(d)) cor(x, d$clima_bio01_sd_res, use = "complete.obs") else NA_real_,
                  r_prec_sd = if ("clima_bio12_sd" %in% names(d)) cor(x, d$clima_bio12_sd, use = "complete.obs") else NA_real_)
    # VIF de la variable dentro del conjunto de 10 predictores (la variable nueva se suma; topo_sd_elev ya está en el conjunto)
    cols <- unique(c(base_10, v)); dv <- d[cols] %>% filter(complete.cases(.))
    out$VIF_en_conjunto <- if (nrow(dv) > length(cols) + 5) unname(vif(dv)[v]) else NA_real_
    out$VIF_topo_sd_elev_con_ella <- if (v != "topo_sd_elev" && nrow(dv) > length(cols) + 5) unname(vif(dv)["topo_sd_elev"]) else NA_real_
    out
  })) %>% mutate(tam_km = tam_km)
  write_csv(res, file.path(dir_mod, sprintf("topo_fina_resumen_%dkm.csv", tam_km)))

  r2 <- function(x) round(x, 2)
  cat("\nCorrelación (Pearson) de las variables topográficas finas con los predictores existentes:\n")
  print(as.data.frame(res %>% mutate(across(where(is.numeric), r2))), row.names = FALSE)
  cat("\nLectura: |r| > 0.8 con topo_sd_elev = no aporta información nueva; r_temp_sd_res cerca de 0 = no arrastra temperatura\n",
      "independiente de la elevación; VIF > 5 = colinealidad alta si se añade al conjunto.\n", sep = "")
}
