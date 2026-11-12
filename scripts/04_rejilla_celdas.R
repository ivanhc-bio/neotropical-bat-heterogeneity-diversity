# 04_rejilla_celdas.R
# Rejilla equivalente (equal-area) anidada y asignación de los especímenes curados a celdas.
#
# Decisiones (5-nov-2026):
#   - Celdas CUADRADAS en proyección Lambert azimutal de áreas iguales (LAEA).
#   - Tamaños anidados 50 / 100 / 200 km, alineados a múltiplos de 200 km. Principal: 100 km.
#   - Se conservan celdas con >= 50 % de tierra; las demás se rescatan solo si tienen registros
#     (islas, costa) y se usa el área de tierra como covariable.
#   - Rangos (Marsh et al. 2022 + mapas propios, MDD v2.5) comparados con registros por celda.
#   - Los registros vienen de especimenes_curados_pais.rds (usar "nombre_auditado").
#     La versión de sensibilidad taxonómica usa "nombre_sensibilidad".
#
# Qué produce (carpeta datos/rejilla):
#   celdas_{50,100,200}km.gpkg / .csv      celdas con esfuerzo, cobertura y covariables de tierra
#   presencia_{tam}km.csv                  celda x especie (n de registros), versión principal
#   presencia_sens_{tam}km.csv             igual con nombre_sensibilidad
#   rango_vs_registros_{tam}km.csv         celda x especie: en polígono / con registro
#   rango_vs_registros_resumen_{tam}km.csv por especie
#   especies_celdas_ocupadas.csv           insumo para el criterio de inclusión de especies
#   registros_sin_celda_{tam}km.csv        registros fuera de la rejilla (mar, fuera de región)
#   resumen_rejilla.csv, mapa_registros_100km.png
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.
# Paquetes nuevos: install.packages(c("sf", "rnaturalearth", "rnaturalearthdata"))

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(sf)
library(rnaturalearth)
library(ggplot2)

# ---- 0. Rutas y parámetros ---------------------------------------------------
dir_base <- RUTA_BASE
dir_gbif <- file.path(dir_base, "datos", "gbif")
dir_cur  <- file.path(dir_base, "datos", "curacion")
dir_rang <- file.path(dir_base, "datos", "rangos")
dir_out  <- file.path(dir_base, "datos", "rejilla")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# PARÁMETROS (propuestos por mí: confírmalos o ajústalos)
tamanos_km        <- c(50, 100, 200)  # anidados; deben dividir a 200
tam_principal_km  <- 100
frac_tierra_min   <- 0.5              # fracción de tierra para conservar la celda
reg_min_rescate   <- 5                # registros mínimos para rescatar una celda con < 50 % de tierra
crs_laea <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

archivo_rangos <- file.path(dir_rang, "rangos_focales_v25.gpkg")
if (!file.exists(archivo_rangos)) stop("Falta ", archivo_rangos)
stopifnot(all(200 %% tamanos_km == 0), tam_principal_km %in% tamanos_km)

# ---- 1. Entradas -------------------------------------------------------------
leer_registros <- function(archivo, col_nombre) {
  x <- readRDS(file.path(dir_cur, archivo))
  req <- c("decimalLatitude", "decimalLongitude", col_nombre, "genero_analisis", "institutionCode")
  faltan <- setdiff(req, names(x))
  if (length(faltan) > 0) stop(archivo, ": faltan columnas ", paste(faltan, collapse = ", "))
  x <- x %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(.data[[col_nombre]])) %>%
    mutate(especie = .data[[col_nombre]])
  pts <- st_as_sf(x, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE) %>%
    st_transform(crs_laea)
  xy <- st_coordinates(pts)
  x$x <- xy[, 1]; x$y <- xy[, 2]
  x$loc_key <- paste(round(x$decimalLatitude, 2), round(x$decimalLongitude, 2))
  x
}

# Si corriste 04b_filtro_rango.R, se usan sus salidas; si no, las de la auditoría de países.
sufijo <- if (file.exists(file.path(dir_cur, "especimenes_curados_rango.rds"))) "rango" else "pais"
cat("Insumo de registros:", sufijo, "\n")
rec  <- leer_registros(sprintf("especimenes_curados_%s.rds", sufijo), "nombre_auditado")
rec_s <- leer_registros(sprintf("especimenes_sensibilidad_%s.rds", sufijo), "nombre_sensibilidad")
cat("\nRegistros principal:", nrow(rec), "| sensibilidad:", nrow(rec_s), "\n")

rangos <- st_read(archivo_rangos, quiet = TRUE) %>% st_make_valid() %>% st_transform(crs_laea)
cat("Polígonos de rango:", nrow(rangos), "\n")
sin_rango <- setdiff(unique(rec$especie), rangos$sciname)
if (length(sin_rango) > 0) cat("Especies con registros pero sin polígono:", paste(sin_rango, collapse = "; "), "\n")

# Tierra: Natural Earth 1:50m (incluye islas medianas; las muy pequeñas se rescatan por registros)
sf_use_s2(FALSE)
tierra <- ne_countries(scale = 50, returnclass = "sf") %>%
  st_make_valid() %>% st_union() %>% st_transform(crs_laea) %>% st_make_valid()
sf_use_s2(TRUE)

# ---- 2. Extensión de la región (alineada a 200 km) -----------------------------
bb_r <- st_bbox(rangos)
ext <- c(xmin = unname(min(bb_r["xmin"], min(rec$x, rec_s$x))),
         ymin = unname(min(bb_r["ymin"], min(rec$y, rec_s$y))),
         xmax = unname(max(bb_r["xmax"], max(rec$x, rec_s$x))),
         ymax = unname(max(bb_r["ymax"], max(rec$y, rec_s$y))))
ext <- c(xmin = floor(ext[["xmin"]] / 2e5) * 2e5, ymin = floor(ext[["ymin"]] / 2e5) * 2e5,
         xmax = ceiling(ext[["xmax"]] / 2e5) * 2e5, ymax = ceiling(ext[["ymax"]] / 2e5) * 2e5)
cat("Extensión (km):", round(ext / 1000), "\n")

# ---- 3. Funciones ------------------------------------------------------------
# ID de celda a partir de índices enteros (sirve para registros y para el anidamiento)
id_celda <- function(ix, iy, tam_km) sprintf("c%d_%d_%d", as.integer(tam_km), as.integer(ix), as.integer(iy))

indices <- function(x, y, tam_km) list(ix = floor(x / (tam_km * 1000)), iy = floor(y / (tam_km * 1000)))

crear_rejilla <- function(tam_km) {
  tam <- tam_km * 1000
  nx <- as.integer(round((ext[["xmax"]] - ext[["xmin"]]) / tam))
  ny <- as.integer(round((ext[["ymax"]] - ext[["ymin"]]) / tam))
  g <- st_make_grid(st_sfc(st_point(c(0, 0)), crs = crs_laea), cellsize = tam,
                    offset = c(ext[["xmin"]], ext[["ymin"]]), n = c(nx, ny), what = "polygons")
  g <- st_sf(geometry = g)
  cen <- st_coordinates(st_centroid(st_geometry(g)))
  ind <- indices(cen[, 1], cen[, 2], tam_km)
  g$ix <- ind$ix; g$iy <- ind$iy
  g$celda <- id_celda(g$ix, g$iy, tam_km)
  g$tam_km <- tam_km
  g$x_centro_km <- cen[, 1] / 1000
  g$y_centro_km <- cen[, 2] / 1000
  # anidamiento
  if (tam_km < 100) g$celda_100 <- id_celda(floor(g$ix / (100 / tam_km)), floor(g$iy / (100 / tam_km)), 100)
  if (tam_km < 200) g$celda_200 <- id_celda(floor(g$ix / (200 / tam_km)), floor(g$iy / (200 / tam_km)), 200)
  # solo celdas que tocan tierra
  g <- g[lengths(st_intersects(g, tierra)) > 0, ]
  # área de tierra por celda
  inter <- suppressWarnings(st_intersection(g[, "celda"], st_sf(geometry = tierra)))
  at <- tapply(as.numeric(st_area(inter)) / 1e6, inter$celda, sum)
  g$area_tierra_km2 <- as.numeric(at[g$celda]); g$area_tierra_km2[is.na(g$area_tierra_km2)] <- 0
  g$area_celda_km2 <- tam_km^2
  g$frac_tierra <- pmin(g$area_tierra_km2 / g$area_celda_km2, 1)
  # coordenadas geográficas del centro (para mapas y covariables)
  cg <- st_coordinates(st_transform(st_centroid(st_geometry(g)), 4326))
  g$lon_centro <- cg[, 1]; g$lat_centro <- cg[, 2]
  g
}

estadisticos_celda <- function(r) {
  # r: registros con columnas celda, especie, genero_analisis, loc_key, institutionCode
  por_esp <- r %>% count(celda, especie, name = "k")
  base <- por_esp %>% group_by(celda) %>%
    summarise(n_reg = sum(k), S_obs = n(), f1 = sum(k == 1), f2 = sum(k == 2), .groups = "drop") %>%
    mutate(chao1 = S_obs + f1 * (f1 - 1) / (2 * (f2 + 1)),   # Chao1 con corrección de sesgo
           cob_good = 1 - f1 / n_reg)                         # cobertura de Good (por registros)
  esf <- r %>% group_by(celda) %>%
    summarise(n_localidades = n_distinct(loc_key), n_instituciones = n_distinct(institutionCode), .groups = "drop")
  gen <- r %>% count(celda, genero_analisis, especie) %>%
    group_by(celda, genero_analisis) %>% summarise(S = n(), .groups = "drop") %>%
    pivot_wider(names_from = genero_analisis, values_from = S, names_prefix = "S_", values_fill = 0)
  gen_n <- r %>% count(celda, genero_analisis, name = "k") %>%
    pivot_wider(names_from = genero_analisis, values_from = k, names_prefix = "n_reg_", values_fill = 0)
  base %>% left_join(esf, by = "celda") %>% left_join(gen, by = "celda") %>% left_join(gen_n, by = "celda")
}

# ---- 4. Bucle por tamaño de celda ----------------------------------------------
resumen <- list(); ocupadas <- list()

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  g <- crear_rejilla(tam_km)

  # asignar registros
  asignar <- function(r) {
    ind <- indices(r$x, r$y, tam_km)
    r$celda <- id_celda(ind$ix, ind$iy, tam_km)
    r
  }
  r  <- asignar(rec)
  rs <- asignar(rec_s)

  est <- estadisticos_celda(r)
  g <- g %>% left_join(est, by = "celda")
  cols_num <- setdiff(names(est), "celda")
  for (cc in cols_num) if (cc %in% names(g) && !cc %in% c("chao1", "cob_good")) g[[cc]][is.na(g[[cc]])] <- 0

  g <- g %>% mutate(
    motivo = case_when(
      frac_tierra >= frac_tierra_min ~ "tierra_>=50%",
      area_tierra_km2 > 0 & n_reg >= reg_min_rescate ~ "rescate_registros",
      TRUE ~ "descartada_mar"),
    mantener = motivo != "descartada_mar",
    sin_registros = n_reg == 0)

  kept <- g %>% filter(mantener)
  r_ok <- r %>% filter(celda %in% kept$celda)
  r_fuera <- r %>% filter(!celda %in% kept$celda)
  cat("Celdas con tierra:", nrow(g), "| conservadas:", nrow(kept),
      "| rescatadas:", sum(kept$motivo == "rescate_registros"),
      "| con registros:", sum(kept$n_reg > 0), "\n")
  cat("Registros asignados:", nrow(r_ok), "| sin celda conservada:", nrow(r_fuera), "\n")

  # guardar celdas
  st_write(g, file.path(dir_out, sprintf("celdas_%dkm.gpkg", tam_km)), delete_dsn = TRUE, quiet = TRUE)
  write_csv(st_drop_geometry(g), file.path(dir_out, sprintf("celdas_%dkm.csv", tam_km)))
  write_csv(r_fuera %>% select(any_of(c("gbifID", "especie", "countryCode", "locality",
                                         "decimalLatitude", "decimalLongitude", "institutionCode", "catalogNumber", "celda"))),
            file.path(dir_out, sprintf("registros_sin_celda_%dkm.csv", tam_km)))

  # presencia (principal y sensibilidad)
  pres <- r_ok %>% count(celda, genero_analisis, especie, name = "n_reg")
  write_csv(pres, file.path(dir_out, sprintf("presencia_%dkm.csv", tam_km)))
  pres_s <- rs %>% filter(celda %in% kept$celda) %>% count(celda, genero_analisis, especie, name = "n_reg")
  write_csv(pres_s, file.path(dir_out, sprintf("presencia_sens_%dkm.csv", tam_km)))

  # rango (polígono) vs registros
  hits <- st_intersects(rangos, kept)
  en_pol <- tibble(especie = rep(rangos$sciname, lengths(hits)), celda = kept$celda[unlist(hits)]) %>%
    mutate(en_poligono = TRUE)
  comp <- full_join(en_pol, pres %>% select(celda, especie, n_reg), by = c("celda", "especie")) %>%
    mutate(en_poligono = coalesce(en_poligono, FALSE), con_registro = !is.na(n_reg), n_reg = coalesce(n_reg, 0L))
  write_csv(comp, file.path(dir_out, sprintf("rango_vs_registros_%dkm.csv", tam_km)))
  comp_res <- comp %>% group_by(especie) %>%
    summarise(n_celdas_poligono = sum(en_poligono), n_celdas_registro = sum(con_registro),
              n_ambas = sum(en_poligono & con_registro),
              n_registro_fuera_poligono = sum(con_registro & !en_poligono),
              pct_poligono_con_registro = round(100 * n_ambas / pmax(n_celdas_poligono, 1), 1),
              .groups = "drop") %>% arrange(desc(n_registro_fuera_poligono))
  write_csv(comp_res, file.path(dir_out, sprintf("rango_vs_registros_resumen_%dkm.csv", tam_km)))

  ocupadas[[as.character(tam_km)]] <- pres %>% group_by(especie, genero_analisis) %>%
    summarise(n_celdas = n(), n_registros = sum(n_reg), .groups = "drop") %>% mutate(tam_km = tam_km)

  resumen[[as.character(tam_km)]] <- tibble(
    tam_km = tam_km, celdas_con_tierra = nrow(g), conservadas = nrow(kept),
    rescatadas_por_registros = sum(kept$motivo == "rescate_registros"),
    descartadas_mar = sum(!g$mantener),
    conservadas_con_registros = sum(kept$n_reg > 0),
    con_ge10_registros = sum(kept$n_reg >= 10), con_ge30_registros = sum(kept$n_reg >= 30),
    con_ge5_especies = sum(kept$S_obs >= 5),
    registros_asignados = nrow(r_ok), registros_sin_celda = nrow(r_fuera))

  if (tam_km == tam_principal_km) {
    n_principal <- nrow(r_ok)
    g_principal <- kept
  }
}

# ---- 5. Resúmenes ---------------------------------------------------------------
resumen <- bind_rows(resumen)
write_csv(resumen, file.path(dir_out, "resumen_rejilla.csv"))
cat("\n"); print(as.data.frame(resumen))

ocupadas <- bind_rows(ocupadas) %>%
  pivot_wider(names_from = tam_km, values_from = c(n_celdas, n_registros), values_fill = 0) %>%
  arrange(n_celdas_100)
write_csv(ocupadas, file.path(dir_out, "especies_celdas_ocupadas.csv"))
cat("\nEspecies con menos de 5 celdas ocupadas a 100 km:",
    sum(ocupadas$n_celdas_100 < 5), "de", nrow(ocupadas), "\n")

# ---- 6. Mapa de verificación -----------------------------------------------------
tryCatch({
  p <- ggplot(g_principal) +
    geom_sf(aes(fill = log10(n_reg + 1)), colour = NA) +
    scale_fill_viridis_c(name = "log10(registros + 1)", option = "magma", direction = -1) +
    labs(title = sprintf("Registros por celda de %d km (celdas conservadas)", tam_principal_km)) +
    theme_minimal()
  ggsave(file.path(dir_out, "mapa_registros_100km.png"), p, width = 7, height = 8, dpi = 200)
}, error = function(e) message("No se pudo hacer el mapa: ", conditionMessage(e)))

# ---- 7. Bitácora -------------------------------------------------------------------
bitacora_archivo <- file.path(dir_gbif, "bitacora_curacion.csv")
nueva <- tibble(
  paso = sprintf("Registros asignados a celdas conservadas de %d km (principal)", tam_principal_km),
  n_registros = n_principal)
bit <- if (file.exists(bitacora_archivo)) {
  bind_rows(read_csv(bitacora_archivo, show_col_types = FALSE) %>% filter(!paso %in% nueva$paso), nueva)
} else nueva
write_csv(bit, bitacora_archivo)
print(bit)
