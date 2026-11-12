# =============================================================================
# config.R
# Rutas del proyecto. Todos los scripts lo cargan con source("config.R") y
# ninguno contiene rutas absolutas.
#
# Estructura esperada, con las carpetas de datos en la raíz del repositorio:
#
#   <raiz>/
#     scripts/      (esta carpeta, con config.R)
#     rejilla/  modelos/  predictores/  curacion/  filogenia/  ambiente/
#     figuras/  tablas/
#
# Por omisión asume que lo corres desde scripts/. Si tu estructura es otra,
# define la variable de entorno CAP4_BASE antes de empezar.
# =============================================================================

RUTA_BASE <- Sys.getenv("CAP4_BASE", unset = NA_character_)
if (is.na(RUTA_BASE) || !nzchar(RUTA_BASE)) {
  RUTA_BASE <- if (dir.exists("rejilla")) "."
               else if (dir.exists(file.path("..", "rejilla"))) ".."
               else "."
}

# alias: los scripts usan uno u otro nombre
base      <- RUTA_BASE
dir_base  <- RUTA_BASE
dir_datos <- RUTA_BASE          # las carpetas de datos cuelgan de la raíz

dir_rej   <- file.path(RUTA_BASE, "rejilla")
dir_mod   <- file.path(RUTA_BASE, "modelos")
dir_cur   <- file.path(RUTA_BASE, "curacion")
dir_amb   <- file.path(RUTA_BASE, "ambiente")
dir_pred  <- file.path(RUTA_BASE, "predictores")
dir_arb   <- file.path(RUTA_BASE, "filogenia")
dir_fig   <- file.path(RUTA_BASE, "figuras")
dir_out   <- file.path(RUTA_BASE, "tablas")

for (d in c(dir_rej, dir_mod, dir_cur, dir_amb, dir_pred,
            dir_arb, dir_fig, dir_out))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Entradas externas. NO están en el repositorio: son rásteres grandes de
# archivos públicos. Descárgalos (ver README, "Data sources") y apunta aquí,
# o define las variables de entorno correspondientes.
# -----------------------------------------------------------------------------

# Capas de CHELSA-TraCE21k ya recortadas a ~5 km: salida del script 00.
RUTA_TRACE <- Sys.getenv("CAP4_TRACE",
                         unset = file.path(dir_amb, "chelsa_trace_5km"))

# CHELSA v2.1 (clima actual, 1981-2010) recortada al Neotrópico.
RUTA_CLIMA_ACTUAL <- Sys.getenv("CAP4_CHELSA",
                                unset = file.path(dir_amb, "chelsa_v21"))

# Capas topográficas de Amatulli et al. (2018), script 20_topo_rugosidad_fina.R.
RUTA_TOPO_AMATULLI <- Sys.getenv("CAP4_TOPO", unset = dir_amb)

# Polígonos de rango de Marsh et al. (2022), 260 MB; no se distribuyen aquí.
RUTA_RANGOS <- Sys.getenv("CAP4_RANGOS", unset = file.path(RUTA_BASE, "rangos"))

# Carpeta temporal de terra. Conviene que esté en un disco con espacio.
RUTA_TMP_TERRA <- Sys.getenv("CAP4_TMP", unset = file.path(tempdir(), "terra"))
dir.create(RUTA_TMP_TERRA, showWarnings = FALSE, recursive = TRUE)
