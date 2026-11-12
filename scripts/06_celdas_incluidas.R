# 06_celdas_incluidas.R
# Decide qué celdas entran en los modelos, a partir de la rejilla (04), los predictores (05) y los criterios
# acordados el 6 de noviembre de 2026.
#
# Entradas (en dir_rej y dir_pred):
#   celdas_{tam}km.gpkg            rejilla con n_reg, S_obs, cob_good y, por género, n_reg_<género> y S_<género>
#   predictores_{tam}km.csv        tabla de predictores de 05
# Salidas (en dir_out):
#   celdas_incluidas_{tam}km.csv   una fila por celda con la bandera de cada criterio, y por género
#   resumen_inclusion.csv          número de celdas por tamaño, criterio y género
#   vif_final_{tam}km.csv          VIF del conjunto final de predictores entre las celdas incluidas
#
# No probado: no tengo R en mi entorno. Corre por secciones y revisa lo que imprime cada una.
# Limitación: la cobertura de Good por género necesita el número de registros por especie y celda, que no está en
# la rejilla; aquí la cobertura (cob_good) es la de todos los géneros juntos, y el criterio por género usa solo el
# número de registros del género (n_reg_<género>). Si quieres cobertura por género, hay que calcularla en 04.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)
library(sf)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_pred <- file.path(dir_base, "datos", "predictores")
dir_out  <- file.path(dir_base, "datos", "modelos")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

tamanos_km <- c(50, 100, 200)

# Criterios (propuestos; confírmalos o ajústalos). Se evalúan los tres para comparar el tamaño de muestra.
criterios <- tribble(
  ~criterio,  ~n_reg_min, ~cob_min, ~n_reg_genero_min,
  "principal", 10,        0.8,      5,
  "estricto",  20,        0.9,      10,
  "laxo",      5,         0.7,      3)
frac_tierra_min <- 0.5     # se excluyen las celdas rescatadas por registros

# Conjunto final de predictores del modelo principal
predictores <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res",
                 "hist_sd_temp", "hist_anom_prec",
                 "clima_bio01_media", "clima_bio12_media", "ctrl_esfuerzo_log", "ctrl_area_tierra_km2")
# media del EVI (control de productividad); se llama así en la tabla de 05
extra_ctrl <- "hab_evi_media"

vif_desde_datos <- function(d) {
  m <- cor(d, use = "complete.obs")
  v <- diag(solve(m))
  tibble(variable = colnames(m), vif = round(v, 2)) %>% arrange(desc(vif))
}

# ---- 1. Bucle por tamaño de celda ---------------------------------------------------------
resumenes <- list()
resumenes_S <- list()

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  rej <- st_read(file.path(dir_rej, sprintf("celdas_%dkm.gpkg", tam_km)), quiet = TRUE) %>% st_drop_geometry()
  pred <- read_csv(file.path(dir_pred, sprintf("predictores_%dkm.csv", tam_km)), show_col_types = FALSE)

  # columnas por género, tal como las escribe 04
  cols_n_gen <- grep("^n_reg_", names(rej), value = TRUE)
  cols_n_gen <- setdiff(cols_n_gen, "n_reg_")
  generos <- sub("^n_reg_", "", cols_n_gen)
  cat("Géneros encontrados en la rejilla:", paste(generos, collapse = ", "), "\n")

  vars_ok <- intersect(c(predictores, if (extra_ctrl %in% names(pred)) extra_ctrl), names(pred))
  if (!extra_ctrl %in% names(pred)) warning("No encontré '", extra_ctrl, "' en predictores_", tam_km, "km.csv; el EVI queda fuera")
  falta <- setdiff(predictores, names(pred))
  if (length(falta)) stop("Faltan predictores en la tabla de 05: ", paste(falta, collapse = ", "))

  cols_S_gen <- paste0("S_", generos)
  if (!all(cols_S_gen %in% names(rej))) stop("Faltan columnas S_<género> en la rejilla: ",
                                             paste(setdiff(cols_S_gen, names(rej)), collapse = ", "))

  d <- pred %>%
    select(celda, frac_tierra, n_reg, S_obs, cob_good, all_of(vars_ok)) %>%
    left_join(rej %>% select(celda, all_of(cols_n_gen), all_of(cols_S_gen)), by = "celda") %>%
    mutate(pred_completos = if_all(all_of(vars_ok), ~ !is.na(.x)),
           tierra_ok = frac_tierra >= frac_tierra_min)

  out <- d %>% select(celda, frac_tierra, n_reg, S_obs, cob_good, all_of(cols_n_gen), all_of(cols_S_gen), pred_completos)

  for (i in seq_len(nrow(criterios))) {
    cr <- criterios[i, ]
    base <- d$tierra_ok & d$pred_completos & !is.na(d$cob_good) &
      d$n_reg >= cr$n_reg_min & d$cob_good >= cr$cob_min
    out[[paste0("incl_", cr$criterio)]] <- base
    for (g in generos) {
      out[[paste0("incl_", cr$criterio, "_", g)]] <- base & d[[paste0("n_reg_", g)]] >= cr$n_reg_genero_min
    }
    n_base <- sum(base)
    por_gen <- vapply(generos, function(g) sum(out[[paste0("incl_", cr$criterio, "_", g)]]), numeric(1))
    cat(sprintf("%-10s todos los géneros: %d celdas | por género: %s\n", cr$criterio, n_base,
                paste(sprintf("%s=%d", generos, por_gen), collapse = ", ")))
    resumenes[[length(resumenes) + 1]] <- tibble(
      tam_km = tam_km, criterio = cr$criterio, grupo = c("todos", generos), n_celdas = c(n_base, por_gen))

    # De las celdas incluidas por género, cuántas tienen al menos 2 y 3 especies observadas del género
    # (mínimo para diversidad filogenética: MPD/PD necesitan S >= 2; SES-MPD suele pedir S >= 3)
    for (g in generos) {
      inc <- out[[paste0("incl_", cr$criterio, "_", g)]]
      Sg <- d[[paste0("S_", g)]]
      resumenes_S[[length(resumenes_S) + 1]] <- tibble(
        tam_km = tam_km, criterio = cr$criterio, genero = g,
        n_incluidas = sum(inc), n_S_ge2 = sum(inc & Sg >= 2, na.rm = TRUE), n_S_ge3 = sum(inc & Sg >= 3, na.rm = TRUE),
        S_max = suppressWarnings(max(Sg[inc], na.rm = TRUE)))
    }
  }

  write_csv(out, file.path(dir_out, sprintf("celdas_incluidas_%dkm.csv", tam_km)))

  # VIF del conjunto final en las celdas del criterio principal
  sel <- d[out$incl_principal, vars_ok]
  v <- vif_desde_datos(sel)
  write_csv(v, file.path(dir_out, sprintf("vif_final_%dkm.csv", tam_km)))
  cat("VIF del conjunto final (criterio principal, n =", nrow(sel), "):\n"); print(as.data.frame(v), row.names = FALSE)

  # Comprobación del sesgo de muestreo hacia relieve (mediana de topo_sd_elev: incluidas frente a toda la tierra)
  cat(sprintf("Mediana de topo_sd_elev: incluidas %.0f | todas las celdas de tierra %.0f\n",
              median(d$topo_sd_elev[out$incl_principal], na.rm = TRUE),
              median(d$topo_sd_elev[d$tierra_ok], na.rm = TRUE)))
}

resumen <- bind_rows(resumenes)
write_csv(resumen, file.path(dir_out, "resumen_inclusion.csv"))
cat("\nResumen guardado en", file.path(dir_out, "resumen_inclusion.csv"), "\n")
print(as.data.frame(resumen %>% pivot_wider(names_from = criterio, values_from = n_celdas)), row.names = FALSE)

# Celdas con riqueza suficiente para diversidad filogenética, por género
resumen_S <- bind_rows(resumenes_S) %>% mutate(pct_S_ge2 = round(100 * n_S_ge2 / n_incluidas), pct_S_ge3 = round(100 * n_S_ge3 / n_incluidas))
write_csv(resumen_S, file.path(dir_out, "resumen_S_minimo.csv"))
cat("\nCeldas incluidas con S >= 2 y S >= 3 por género (criterio principal):\n")
print(as.data.frame(resumen_S %>% filter(criterio == "principal")), row.names = FALSE)
