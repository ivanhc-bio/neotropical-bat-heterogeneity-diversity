# 16_sensibilidad_datos.R
# Fase 4 (parte A) del plan: robustez de las pruebas confirmatorias (H1, H2, H3) a decisiones sobre los DATOS, sin recorrer de
# nuevo el flujo de curación. Usa solo lo que ya existe en disco. Todo con el tratamiento espacial "actual" del script 15
# (suavizado s(x, y) por género con k = n / 8), que es la referencia de las pruebas confirmatorias.
#
# Variantes de riqueza (H1: topo_sd_elev -> riqueza, pendiente común; H2: interacción género x topo_sd_elev):
#   principal        criterio principal, taxonomía principal (debe reproducir el script 15: H1 = 0.115)
#   estricto, laxo   criterios de inclusión de celdas del script 06 (n_reg >= 20 / >= 5, cobertura 0.9 / 0.7, n_reg_género 10 / 3)
#   taxonomia_alt    riqueza con la taxonomía de sensibilidad (presencia_sens_{tam}km.csv: A. intermedius separada,
#                    S. lilium del norte de Sudamérica a S. giannae), mismas celdas del criterio principal
# Variantes de SES-MPD (H3: hist_sd_temp -> SES-MPD, pendiente común):
#   principal        SES con los árboles completados (debe reproducir el script 15: H3 = -0.049)
#   sin_hermanas     SES con árboles completados sin hermanas (especies injertadas al azar dentro del género)
#   solo_adn         SES solo con especies con ADN (sin imputadas ni injertadas); celdas con S < 3 tras quitarlas quedan fuera
# Ambas respuestas:
#   sin_<género>     se quita un género frugívoro a la vez: ¿el efecto común depende de uno solo (p. ej. Sturnira)?
#
# Parte B (riqueza, solo H1 y H2): la presencia celda x especie se rehace desde los registros con otro filtro, SIN volver a correr
# 04-09 (las celdas incluidas y los predictores se mantienen; solo cambia S):
#   control_registros                debe reproducir a "principal" (comprueba la asignación a celdas)
#   rango_estricto                   sin registros a más de 100 km del polígono de rango (especimenes_curados_rango_estricto.rds)
#   sin_capitales                    sin los registros marcados por la prueba de capitales (flag_cc_cap)
#   rango_estricto_sin_capitales     ambos filtros
# Lo que NO cubre: criterios estricto/laxo para el SES y árboles de Upham sin injertos (la variante solo_adn los cubre en gran parte).
#
# Salida (datos/modelos): robustez_datos_{tam}km.csv, con una fila por variante x respuesta.
# No probado: no tengo R en mi entorno. Prueba con tamanos_km <- 100 y revisa que la variante "principal" reproduzca el script 15.
# Requiere: mgcv, dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv)
library(dplyr)
library(readr)
library(tidyr)
library(sf)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_cur  <- file.path(dir_base, "datos", "curacion")
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_pred <- file.path(dir_base, "datos", "predictores")

tamanos_km   <- c(50, 100, 200)
frugivoros   <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
n_min_modelo <- 25
k_por_celdas <- 8
esf_origen   <- "esf_mam_log"     # esfuerzo independiente, igual que 09_bases_modelos.R
crs_laea     <- "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"   # igual que 04_rejilla_celdas.R

# Fuentes de presencia para la riqueza (parte B: se rehace celda x especie desde los registros, con las mismas celdas incluidas)
#   control_registros       registros de especimenes_curados_rango.rds (lo que usó 04): debe dar la misma presencia que presencia_{tam}km.csv
#   rango_estricto          especimenes_curados_rango_estricto.rds (04b): sin registros a mas de buffer_km del poligono
#   sin_capitales           registros curados sin los marcados por la prueba de capitales (flag_cc_cap)
#   estricto_sin_capitales  ambos filtros a la vez
fuentes_registros <- list(
  control_registros      = list(archivo = "especimenes_curados_rango.rds",          quitar_capitales = FALSE),
  rango_estricto         = list(archivo = "especimenes_curados_rango_estricto.rds", quitar_capitales = FALSE),
  sin_capitales          = list(archivo = "especimenes_curados_rango.rds",          quitar_capitales = TRUE),
  estricto_sin_capitales = list(archivo = "especimenes_curados_rango_estricto.rds", quitar_capitales = TRUE))

hip      <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp")
ctrl_riq <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2", "clima_bio01_media", "clima_bio12_media")
ctrl_ses <- c("ctrl_esfuerzo_log", "clima_bio01_media", "clima_bio12_media", "log_S")
respuestas <- list(
  riqueza = list(y = "cbind(S, S_fallo)", familia = "binomial", ctrl = ctrl_riq, h_conf = "topo_sd_elev", h_nombre = "H1"),
  ses_mpd = list(y = "ses_mpd", familia = "gaussian", ctrl = ctrl_ses, h_conf = "hist_sd_temp", h_nombre = "H3"))

# ---- 1. Funciones (las mismas que en el script 15, solo con el tratamiento "actual") ----------------
zscore <- function(x) as.numeric(scale(x))

preparar <- function(d, vars) {
  gens <- sort(unique(d$genero))
  d <- d %>% group_by(genero) %>% mutate(across(all_of(vars), zscore)) %>% ungroup()
  for (g in gens) d[[paste0("I_", g)]] <- as.numeric(d$genero == g)
  d$genero <- factor(d$genero)
  d
}

k_por_genero <- function(d) {
  n <- table(d$genero)
  setNames(pmax(5, pmin(40, floor(as.numeric(n) / k_por_celdas))), names(n))
}

formula_conjunta <- function(y, d, ctrl, separados, comunes_hip) {
  k <- k_por_genero(d)
  esp <- sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k), as.integer(k))
  partes <- c("genero", ctrl, comunes_hip, if (length(separados)) paste0("genero:", separados), esp)
  as.formula(paste(y, "~", paste(partes, collapse = " + ")))
}

ajusta <- function(f, d, fam) tryCatch(gam(f, family = fam, data = d, method = "REML"), error = function(e) NULL)

# Ajusta una variante: pendiente común de la hipótesis confirmatoria, pendientes por género y prueba de Wald de H2
ajustar_variante <- function(d0, rn, variante, excluir = character()) {
  R <- respuestas[[rn]]
  vars <- c(hip, R$ctrl)
  d0 <- d0 %>% filter(genero %in% setdiff(frugivoros, excluir)) %>% filter(if_all(all_of(vars), ~ !is.na(.x)))
  ok <- names(which(table(d0$genero) >= n_min_modelo))
  if (length(ok) < 2) return(NULL)
  d <- preparar(d0 %>% filter(genero %in% ok), vars)
  fam <- if (R$familia == "binomial") binomial() else gaussian()

  m_com <- ajusta(formula_conjunta(R$y, d, R$ctrl, character(), hip), d, fam)
  m_gen <- ajusta(formula_conjunta(R$y, d, R$ctrl, hip, character()), d, fam)
  m_wal <- ajusta(formula_conjunta(R$y, d, R$ctrl, R$h_conf, hip), d, fam)
  if (is.null(m_com)) return(NULL)

  tc <- as.data.frame(summary(m_com)$p.table)[R$h_conf, ]
  por_gen <- NA_character_
  if (!is.null(m_gen)) {
    tg <- as.data.frame(summary(m_gen)$p.table); tg$nombre <- rownames(tg)
    tg <- tg[grepl(paste0("^genero.*:", R$h_conf, "$"), tg$nombre), ]
    por_gen <- paste(sprintf("%s %.2f (p %.3f)", sub("^genero([^:]+):.*$", "\\1", tg$nombre), tg$Estimate, tg[["Pr(>|t|)"]] %||% tg[["Pr(>|z|)"]]),
                     collapse = "; ")
  }
  p_h2 <- NA_real_
  if (!is.null(m_wal)) {
    at <- tryCatch(anova(m_wal)$pTerms.table, error = function(e) NULL)
    fila <- if (is.null(at)) integer() else which(rownames(at) == paste0("genero:", R$h_conf))
    if (length(fila) == 1) p_h2 <- at[fila, "p-value"]
  }
  tibble(variante = variante, respuesta = rn, prueba = R$h_nombre, predictor = R$h_conf,
         generos = paste(ok, collapse = ","), n = nrow(d),
         estimate = tc[[1]], se = tc[[2]], p = tc[[4]],
         ic_inf = tc[[1]] - 1.96 * tc[[2]], ic_sup = tc[[1]] + 1.96 * tc[[2]],
         p_H2_wald = p_h2, por_genero = por_gen)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Presencia celda x genero x especie rehecha desde los registros (misma asignacion a celdas que 04: indices enteros en LAEA)
presencia_desde_registros <- function(tam_km, archivo, quitar_capitales) {
  x <- readRDS(file.path(dir_cur, archivo))
  req <- c("decimalLatitude", "decimalLongitude", "nombre_auditado", "genero_analisis", if (quitar_capitales) "flag_cc_cap")
  if (length(setdiff(req, names(x)))) stop(archivo, ": faltan columnas ", paste(setdiff(req, names(x)), collapse = ", "))
  x <- x %>% filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(nombre_auditado))
  if (quitar_capitales) x <- x %>% filter(!(flag_cc_cap %in% TRUE))
  xy <- st_coordinates(st_transform(st_as_sf(x, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE), crs_laea))
  tam <- tam_km * 1000
  x$celda <- sprintf("c%d_%d_%d", as.integer(tam_km), as.integer(floor(xy[, 1] / tam)), as.integer(floor(xy[, 2] / tam)))
  x %>% count(celda, genero_analisis, especie = nombre_auditado, name = "n_reg")
}

leer_presencia <- function(tam_km, fuente) {
  if (fuente %in% names(fuentes_registros)) {
    f <- fuentes_registros[[fuente]]
    return(presencia_desde_registros(tam_km, f$archivo, f$quitar_capitales))
  }
  read_csv(file.path(dir_rej, sprintf(if (fuente == "alt") "presencia_sens_%dkm.csv" else "presencia_%dkm.csv", tam_km)),
           show_col_types = FALSE)
}

# Datos de riqueza para un criterio de celdas y una fuente de presencia ("principal", "alt" o una de fuentes_registros)
datos_riqueza <- function(tam_km, criterio, taxonomia) {
  inc  <- read_csv(file.path(dir_mod, sprintf("celdas_incluidas_%dkm.csv", tam_km)), show_col_types = FALSE) %>%
    select(celda, starts_with("incl_"))
  cel  <- read_csv(file.path(dir_rej, sprintf("celdas_%dkm.csv", tam_km)), show_col_types = FALSE) %>% select(celda, x_centro_km, y_centro_km)
  pres <- leer_presencia(tam_km, taxonomia)
  vars_pred <- unique(c(hip, setdiff(ctrl_riq, "ctrl_esfuerzo_log")))
  pred <- read_csv(file.path(dir_pred, sprintf("predictores_%dkm.csv", tam_km)), show_col_types = FALSE) %>% select(celda, all_of(vars_pred))
  esf  <- read_csv(file.path(dir_pred, sprintf("esfuerzo_independiente_%dkm.csv", tam_km)), show_col_types = FALSE) %>%
    select(celda, ctrl_esfuerzo_log = all_of(esf_origen))
  todo <- inc %>% left_join(cel, by = "celda") %>% left_join(pred, by = "celda") %>% left_join(esf, by = "celda")
  bind_rows(lapply(frugivoros, function(g) {
    d <- todo[todo[[paste0("incl_", criterio, "_", g)]] %in% TRUE, ]
    p <- pres %>% filter(genero_analisis == g, celda %in% d$celda)
    Sg <- p %>% group_by(celda) %>% summarise(S = n_distinct(especie), .groups = "drop")
    pool <- n_distinct(p$especie)
    d %>% select(celda, x_centro_km, y_centro_km, all_of(c(vars_pred, "ctrl_esfuerzo_log"))) %>%
      left_join(Sg, by = "celda") %>% mutate(S = coalesce(S, 0L), S_fallo = pool - S, genero = g) %>% filter(S_fallo >= 0)
  }))
}

# ---- 2. Por tamaño de celda --------------------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE) %>%
    filter(genero %in% frugivoros)
  res <- list()

  # --- riqueza: criterios, taxonomía y quitar un género
  variantes_riq <- list(principal = c("principal", "principal"), estricto = c("estricto", "principal"),
                        laxo = c("laxo", "principal"), taxonomia_alt = c("principal", "alt"),
                        control_registros = c("principal", "control_registros"),
                        rango_estricto = c("principal", "rango_estricto"),
                        sin_capitales = c("principal", "sin_capitales"),
                        rango_estricto_sin_capitales = c("principal", "estricto_sin_capitales"))
  d_riq <- lapply(variantes_riq, function(v) datos_riqueza(tam_km, v[1], v[2]))
  # control: rehacer la presencia desde los registros debe dar lo mismo que presencia_{tam}km.csv
  p_ref <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE) %>% filter(genero_analisis %in% frugivoros)
  p_reg <- leer_presencia(tam_km, "control_registros") %>% filter(genero_analisis %in% frugivoros)
  cat(sprintf("Control de la asignacion a celdas: combinaciones celda-especie en presencia_%dkm.csv = %d; rehechas desde registros = %d; distintas = %d\n",
              tam_km, nrow(p_ref), nrow(p_reg), nrow(anti_join(p_ref, p_reg, by = c("celda", "especie"))) + nrow(anti_join(p_reg, p_ref, by = c("celda", "especie")))))
  for (id in c("rango_estricto", "sin_capitales", "estricto_sin_capitales")) {
    n_id <- nrow(leer_presencia(tam_km, id) %>% filter(genero_analisis %in% frugivoros))
    cat(sprintf("  %-28s combinaciones celda-especie: %d (%.1f %% de las de control)\n", id, n_id, 100 * n_id / nrow(p_reg)))
  }
  for (id in names(d_riq)) res[[paste("riqueza", id)]] <- ajustar_variante(d_riq[[id]], "riqueza", id)
  for (g in frugivoros) res[[paste("riqueza sin", g)]] <- ajustar_variante(d_riq$principal, "riqueza", paste0("sin_", g), excluir = g)
  n_base <- sum(base$genero %in% frugivoros)
  cat(sprintf("Comprobación: celdas de riqueza en 'principal' = %d; en base_modelos = %d (deben coincidir salvo NA en predictores)\n",
              nrow(d_riq$principal), n_base))

  # --- SES-MPD: árboles alternativos y quitar un género
  col_ses <- c(principal = "ses_mpd", sin_hermanas = "ses_mpd_sin_hermanas", solo_adn = "ses_mpd_solo_adn")
  for (id in names(col_ses)) {
    if (!col_ses[[id]] %in% names(base)) { cat("  Falta la columna", col_ses[[id]], "; se omite", id, "\n"); next }
    d <- base %>% mutate(ses_mpd = .data[[col_ses[[id]]]]) %>% filter(!is.na(ses_mpd))
    res[[paste("ses", id)]] <- ajustar_variante(d, "ses_mpd", id)
  }
  d_ses <- base %>% filter(filo_valida, !is.na(ses_mpd))
  for (g in frugivoros) res[[paste("ses sin", g)]] <- ajustar_variante(d_ses, "ses_mpd", paste0("sin_", g), excluir = g)

  res <- bind_rows(res) %>% mutate(tam_km = tam_km)
  write_csv(res, file.path(dir_mod, sprintf("robustez_datos_%dkm.csv", tam_km)))

  # ---- 3. Resumen legible ----
  r3 <- function(x) round(x, 3)
  for (rn in names(respuestas)) {
    rr <- res %>% filter(respuesta == rn)
    if (!nrow(rr)) next
    cat(sprintf("\n--- %s: pendiente común de %s (%s) por variante de datos ---\n", rn, respuestas[[rn]]$h_conf, respuestas[[rn]]$h_nombre))
    print(as.data.frame(rr %>% mutate(estimado = sprintf("%.3f (%.3f)", estimate, se),
                                      ic95 = sprintf("[%.3f, %.3f]", ic_inf, ic_sup), p = r3(p), p_H2_wald = r3(p_H2_wald)) %>%
            select(variante, n, estimado, ic95, p, p_H2_wald)), row.names = FALSE)
    cat(sprintf("\n--- %s: pendientes por género de %s ---\n", rn, respuestas[[rn]]$h_conf))
    print(as.data.frame(rr %>% select(variante, por_genero)), row.names = FALSE)
  }
  cat("\nGuardado:", file.path(dir_mod, sprintf("robustez_datos_%dkm.csv", tam_km)), "\n")
  cat("Referencia (script 15, tratamiento actual): H1 = 0.115 (0.037); H3 = -0.049 (0.045); H2 riqueza Wald p = .001\n")
}
