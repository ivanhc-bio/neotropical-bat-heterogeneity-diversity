# =============================================================================
# 31_tabla_s1.R
# Reúne los conteos que faltan para la Tabla S1 del material suplementario:
#   (1) cuántos registros quedan al final de cada paso de la curación
#   (2) cuántas celdas ocupa cada género en cada tamaño de celda
#
# Salidas (en <raiz del repositorio>/tablas):
#   tabla_s1_pasos.csv
#   tabla_s1_generos.csv
#
# El script solo lee; no modifica nada.
# =============================================================================

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

base    <- RUTA_BASE
dir_cur <- file.path(base, "datos", "curacion")
dir_rej <- file.path(base, "datos", "rejilla")
dir_out <- file.path(base, "tablas")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

col_de <- function(d, candidatas) {
  hit <- intersect(candidatas, names(d))
  if (length(hit) == 0) NA_character_ else hit[1]
}

# -----------------------------------------------------------------------------
# 1. Registros al final de cada paso de la curación
# -----------------------------------------------------------------------------

etapas <- c(
  "Records retrieved and standardized" = "especimenes_con_marcas.rds",
  "After the coordinate filters"       = "especimenes_curados_coord.rds",
  "With a country assigned"            = "especimenes_con_pais.rds",
  "After the country filter"           = "especimenes_curados_pais.rds",
  "After the altitudinal range filter" = "especimenes_curados_rango.rds"
)

pasos <- lapply(seq_along(etapas), function(i) {

  nombre <- names(etapas)[i]
  ruta   <- file.path(dir_cur, etapas[[i]])

  if (!file.exists(ruta)) {
    message("AVISO: no encontré ", ruta)
    return(NULL)
  }

  d <- readRDS(ruta)

  cat("\n--- ", nombre, " ---\n", sep = "")
  cat("filas: ", nrow(d), "   columnas: ", ncol(d), "\n", sep = "")
  cat("nombres: ", paste(names(d), collapse = ", "), "\n", sep = "")

  c_sp <- col_de(d, c("especie", "especie_final", "species", "scientificName"))
  c_ge <- col_de(d, c("genero_analisis", "genero", "genus"))

  tibble(
    paso      = nombre,
    archivo   = etapas[[i]],
    registros = nrow(d),
    especies  = if (!is.na(c_sp)) n_distinct(d[[c_sp]]) else NA_integer_,
    generos   = if (!is.na(c_ge)) n_distinct(d[[c_ge]]) else NA_integer_
  )
}) |> bind_rows()

if (nrow(pasos) > 0) {
  pasos <- pasos |>
    mutate(eliminados_en_este_paso = c(NA_integer_, -diff(registros)))

  cat("\n=== Registros por paso ===\n")
  print(as.data.frame(pasos))
  write_csv(pasos, file.path(dir_out, "tabla_s1_pasos.csv"))
}

# -----------------------------------------------------------------------------
# 2. Celdas ocupadas por cada género
# -----------------------------------------------------------------------------

esp <- read_csv(file.path(dir_rej, "especies_celdas_ocupadas.csv"),
                show_col_types = FALSE)

celdas_genero <- lapply(c(50, 100, 200), function(tam) {

  ruta <- file.path(dir_rej, paste0("presencia_", tam, "km.csv"))

  if (!file.exists(ruta)) {
    message("AVISO: no encontré ", ruta)
    return(NULL)
  }

  p <- read_csv(ruta, show_col_types = FALSE)

  cat("\n--- presencia_", tam, "km.csv ---\n", sep = "")
  cat("filas: ", nrow(p), "   columnas: ", ncol(p), "\n", sep = "")
  cat("nombres: ", paste(utils::head(names(p), 10), collapse = ", "),
      if (ncol(p) > 10) " ..." else "", "\n", sep = "")

  c_celda <- col_de(p, c("celda", "id_celda", "cell", "celda_id"))
  c_sp    <- col_de(p, c("especie", "especie_final", "species"))

  if (is.na(c_celda)) {
    message("No reconozco la columna de celda en ", ruta, "; lo salto.")
    return(NULL)
  }

  if (!is.na(c_sp)) {
    # formato largo: una fila por celda y especie
    largo <- p |> select(celda = all_of(c_celda), especie = all_of(c_sp))
  } else {
    # formato ancho: una columna por especie
    largo <- p |>
      pivot_longer(-all_of(c_celda), names_to = "especie", values_to = "pres") |>
      filter(!is.na(pres), pres > 0) |>
      select(celda = all_of(c_celda), especie)
  }

  largo |>
    left_join(select(esp, especie, genero_analisis), by = "especie") |>
    filter(!is.na(genero_analisis)) |>
    distinct(genero_analisis, celda) |>
    count(genero_analisis, name = "celdas") |>
    mutate(tam_km = tam)
}) |> bind_rows()

generos <- esp |>
  group_by(genero_analisis) |>
  summarise(
    especies  = n_distinct(especie),
    registros = sum(n_registros_100, na.rm = TRUE),
    .groups   = "drop"
  )

if (nrow(celdas_genero) > 0) {
  generos <- generos |>
    left_join(
      celdas_genero |>
        pivot_wider(names_from = tam_km, values_from = celdas,
                    names_prefix = "celdas_"),
      by = "genero_analisis"
    )
}

generos <- arrange(generos, genero_analisis)

cat("\n=== Especies, registros y celdas por género ===\n")
print(as.data.frame(generos))
write_csv(generos, file.path(dir_out, "tabla_s1_generos.csv"))

cat("\nListo. Revisa:\n",
    "  ", file.path(dir_out, "tabla_s1_pasos.csv"), "\n",
    "  ", file.path(dir_out, "tabla_s1_generos.csv"), "\n", sep = "")
