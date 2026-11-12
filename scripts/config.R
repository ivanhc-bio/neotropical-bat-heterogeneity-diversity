# =============================================================================
# config.R
# Rutas del proyecto. Todos los scripts lo cargan con source("config.R") y
# ninguno contiene rutas absolutas.
#
# Por omisión asume que lo corres desde la carpeta scripts/, con los datos en
# ../datos. Si tu estructura es otra, define la variable de entorno CAP4_BASE
# antes de empezar, o edita RUTA_BASE aquí abajo.
# =============================================================================

RUTA_BASE <- Sys.getenv("CAP4_BASE", unset = NA_character_)
if (is.na(RUTA_BASE) || !nzchar(RUTA_BASE)) {
  RUTA_BASE <- if (dir.exists("datos")) "."
               else if (dir.exists(file.path("..", "datos"))) ".."
               else "."
}

# alias: los scripts usan uno u otro nombre
base     <- RUTA_BASE
dir_base <- RUTA_BASE

dir_datos  <- file.path(RUTA_BASE, "datos")
dir_rej    <- file.path(dir_datos, "rejilla")
dir_mod    <- file.path(dir_datos, "modelos")
dir_cur    <- file.path(dir_datos, "curacion")
dir_amb    <- file.path(dir_datos, "ambiente")
dir_pred   <- file.path(dir_datos, "predictores")
dir_arb    <- file.path(dir_datos, "filogenia")
dir_fig    <- file.path(RUTA_BASE, "figuras")
dir_out    <- file.path(RUTA_BASE, "tablas")

for (d in c(dir_datos, dir_rej, dir_mod, dir_cur, dir_amb, dir_pred,
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

# Capas topográficas de Amatulli et al. (2018) usadas por 20_topo_rugosidad_fina.R.
RUTA_TOPO_AMATULLI <- Sys.getenv("CAP4_TOPO", unset = dir_amb)

# Carpeta temporal de terra. Conviene que esté en un disco con espacio.
RUTA_TMP_TERRA <- Sys.getenv("CAP4_TMP", unset = file.path(tempdir(), "terra"))
dir.create(RUTA_TMP_TERRA, showWarnings = FALSE, recursive = TRUE)
