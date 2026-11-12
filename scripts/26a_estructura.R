# =============================================================================
# 26a_estructura.R
# Solo imprime la estructura de cuatro archivos, para poder escribir el script
# que arma la variante de "rango estricto". No calcula ni escribe nada.
# =============================================================================

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

suppressPackageStartupMessages({library(readr); library(dplyr)})

base    <- RUTA_BASE
archivos <- c(
  file.path(base, "datos", "rejilla", "registros_fuera_rango.csv"),
  file.path(base, "datos", "rejilla", "rango_vs_registros_100km.csv"),
  file.path(base, "datos", "rejilla", "presencia_100km.csv"),
  file.path(base, "datos", "modelos", "base_modelos_100km.csv")
)

for (ruta in archivos) {
  cat("\n\n=====================================================\n")
  cat(basename(ruta), "\n")
  cat("=====================================================\n")
  if (!file.exists(ruta)) { cat("NO EXISTE en", dirname(ruta), "\n"); next }
  d <- read_csv(ruta, show_col_types = FALSE, progress = FALSE)
  cat("filas:", nrow(d), "   columnas:", ncol(d), "\n\n")
  cat("columnas:\n")
  print(data.frame(columna = names(d), tipo = vapply(d, function(x) class(x)[1], "")),
        row.names = FALSE)
  cat("\nprimeras 3 filas:\n")
  print(as.data.frame(head(d, 3)), row.names = FALSE)
}

cat("\n\nListo. Copia y pega todo lo que imprimió.\n")
