# =============================================================================
# 34_carollia.R
# Examina las dos rarezas de Carollia que quedaron sin probar:
#
#   (A) Autocorrelación espacial residual en el SES-MPD. El diagnóstico del
#       script 27 da una sola I de Moran por género. Aquí se desglosa por
#       bandas de distancia, para ver A QUÉ ESCALA vive la autocorrelación, y
#       se vuelve a ajustar con una base espacial más grande para distinguir
#       entre "falta estructura en el modelo" y "el suavizado quedó corto".
#
#   (B) El efecto histórico sobre la riqueza. La pendiente común es nula, pero
#       la de Carollia no se examinó por separado. Aquí se estima su pendiente
#       en seis especificaciones y se prueba su desviación respecto de la común.
#
# IMPORTANTE: todo lo de (B) es EXPLORATORIO. No entró en las 16 pruebas de
# Holm. Si algo sale, se reporta como hallazgo post hoc, no como prueba.
#
# Salidas (datos/modelos):
#   carollia_moran_bandas.csv
#   carollia_ses_base_grande.csv
#   carollia_riqueza_historica.csv
#
# Requiere: mgcv, dplyr, readr, tibble, ape. Usa las funciones del script 20.
# =============================================================================

solo_funciones <- TRUE
source("20_alfa_pruebas.R")     # funciones, rutas, hip, ctrl_*, generos, tamanos_km

suppressPackageStartupMessages({library(tibble)})
if (!requireNamespace("ape", quietly = TRUE))
  stop("Falta el paquete ape. Instálalo con install.packages('ape').")

tamanos_km <- c(50, 100, 200)
FOCAL      <- "Carollia"

# ---- datos de una escala: copia literal de la del script 27 -----------------

datos_alfa <- function(base, pres, tam_km, generos) {
  h <- unname(hip)
  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g); if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>%
      summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    d
  }))
  vr <- c(h, ctrl_riq)
  d0 <- base %>% filter(if_all(all_of(vr), ~ !is.na(.x)))
  ok <- names(which(table(d0$genero) >= n_min_modelo))
  dr <- preparar(d0 %>% filter(genero %in% ok), vr)

  vs <- c(h, ctrl_ses)
  d1 <- base %>% filter(filo_valida, !is.na(ses_mpd), if_all(all_of(vs), ~ !is.na(.x)))
  ok <- names(which(table(d1$genero) >= n_min_modelo))
  ds <- preparar(d1 %>% filter(genero %in% ok), vs)

  list(riqueza = list(d = dr, ctrl = ctrl_riq, y = "cbind(S, S_fallo)", fam = binomial()),
       ses_mpd = list(d = ds, ctrl = ctrl_ses, y = "ses_mpd",          fam = gaussian()))
}

# ---- suavizados con base ajustable (el del script 20 fija divisor 8, tope 40)

suaves_k <- function(d, divisor, kmax) {
  n <- table(d$genero)
  k <- setNames(pmax(5, pmin(kmax, floor(as.numeric(n) / divisor))), names(n))
  sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)",
          names(k), as.integer(k))
}

fml_k <- function(y, comunes, separados, d, divisor, kmax, espacial = TRUE) {
  piezas <- c("genero", comunes,
              if (length(separados)) paste0("genero:", separados),
              if (espacial) suaves_k(d, divisor, kmax))
  as.formula(paste(y, "~", paste(piezas, collapse = " + ")))
}

# =============================================================================
# (A) AUTOCORRELACIÓN RESIDUAL DEL SES-MPD, POR BANDAS DE DISTANCIA
# =============================================================================

# Una sola I de Moran mezcla todas las distancias y no dice si el problema está
# en celdas vecinas (suavizado corto) o a escala continental (falta un proceso).

moran_bandas <- function(d, resid,
                         cortes = c(0, 100, 200, 400, 800, 1600, Inf)) {
  bind_rows(lapply(levels(d$genero), function(g) {
    i <- which(d$genero == g)
    if (length(i) < 20) return(NULL)
    xy <- as.matrix(d[i, c("x_centro_km", "y_centro_km")])
    D  <- as.matrix(dist(xy))           # km, porque las coordenadas ya están en km
    r  <- resid[i]
    bind_rows(lapply(seq_len(length(cortes) - 1), function(b) {
      w <- (D > cortes[b] & D <= cortes[b + 1]) * 1
      diag(w) <- 0
      etiqueta <- sprintf("%g-%g", cortes[b], cortes[b + 1])
      if (sum(w) < 20)
        return(tibble(genero = g, n_g = length(i), banda_km = etiqueta,
                      pares = sum(w) / 2, moran_I = NA_real_, moran_p = NA_real_))
      mi <- tryCatch(ape::Moran.I(r, w, scaled = TRUE),
                     error = function(e) list(observed = NA_real_, p.value = NA_real_))
      tibble(genero = g, n_g = length(i), banda_km = etiqueta,
             pares = sum(w) / 2, moran_I = mi$observed, moran_p = mi$p.value)
    }))
  }))
}

bandas_todo <- list(); base_grande_todo <- list()

for (tam_km in tamanos_km) {

  cat("\n\n==================  (A) SES-MPD,", tam_km, "km  ==================\n")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  dd   <- datos_alfa(base, pres, tam_km, generos)
  d    <- dd$ses_mpd$d
  if (is.null(d) || !nrow(d)) { message("sin datos de SES-MPD a ", tam_km, " km"); next }

  h <- unname(hip)

  # modelo principal, igual que el del script 27
  m0 <- ajusta(fml(dd$ses_mpd$y, c(dd$ses_mpd$ctrl, h), character(), d), d, dd$ses_mpd$fam)
  r0 <- residuals(m0, type = "deviance")

  cat("\n-- I de Moran por banda, modelo principal --\n")
  b0 <- moran_bandas(d, r0) %>% mutate(modelo = "principal", tam_km = tam_km)
  print(as.data.frame(b0 %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

  # base espacial más grande: divisor 4 en vez de 8, tope 100 en vez de 40
  mG <- ajusta(fml_k(dd$ses_mpd$y, c(dd$ses_mpd$ctrl, h), character(), d,
                     divisor = 4, kmax = 100), d, dd$ses_mpd$fam)
  rG <- residuals(mG, type = "deviance")

  cat("\n-- I de Moran por banda, base espacial grande --\n")
  bG <- moran_bandas(d, rG) %>% mutate(modelo = "base_grande", tam_km = tam_km)
  print(as.data.frame(bG %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

  bandas_todo[[as.character(tam_km)]] <- bind_rows(b0, bG)

  # ¿cambian las conclusiones del SES-MPD con la base grande?
  cmp <- bind_rows(lapply(c("principal", "base_grande"), function(etq) {
    m <- if (etq == "principal") m0 else mG
    tt <- as.data.frame(summary(m)$p.table)
    tibble(tam_km = tam_km, modelo = etq, predictor = h,
           estimado = tt[h, 1], se = tt[h, 2], p = tt[h, 4],
           ic_inf = tt[h, 1] - 1.96 * tt[h, 2], ic_sup = tt[h, 1] + 1.96 * tt[h, 2],
           dev_expl = summary(m)$dev.expl, edf_total = sum(m$edf))
  }))
  cat("\n-- Pendientes comunes del SES-MPD bajo las dos bases --\n")
  print(as.data.frame(cmp %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
  base_grande_todo[[as.character(tam_km)]] <- cmp
}

bandas_todo      <- bind_rows(bandas_todo)
base_grande_todo <- bind_rows(base_grande_todo)
write_csv(bandas_todo,      file.path(dir_mod, "carollia_moran_bandas.csv"))
write_csv(base_grande_todo, file.path(dir_mod, "carollia_ses_base_grande.csv"))

cat("\n\n>>> Cómo leer (A):\n",
    "  - Si la I de Carollia es alta SOLO en la banda 0-100 o 0-200 km y baja con la\n",
    "    base grande, el suavizado quedó corto y no hay nada sustantivo que reportar.\n",
    "  - Si sigue alta a 400 km o más, y con la base grande, hay estructura espacial\n",
    "    que el modelo no captura: eso sí va a las limitaciones.\n",
    "  - Compara también si las pendientes del SES-MPD cambian de signo o de\n",
    "    significancia entre las dos bases.\n", sep = "")

# =============================================================================
# (B) EFECTO HISTÓRICO SOBRE LA RIQUEZA EN CAROLLIA
# =============================================================================

pred_hist <- unname(hip["historica"])   # "hist_sd_temp"
riq_todo  <- list()

saca <- function(m, termino, etq, tam_km, extra = list()) {
  tt <- as.data.frame(summary(m)$p.table)
  if (!(termino %in% rownames(tt)))
    return(tibble(tam_km = tam_km, especificacion = etq, termino = termino,
                  estimado = NA_real_, se = NA_real_, p = NA_real_,
                  ic_inf = NA_real_, ic_sup = NA_real_))
  cw <- tryCatch({
    cc <- concurvity(m, full = TRUE)
    j  <- grep(pred_hist, colnames(cc), fixed = TRUE)
    if (length(j)) max(cc["worst", j]) else NA_real_
  }, error = function(e) NA_real_)
  tibble(tam_km = tam_km, especificacion = etq, termino = termino,
         estimado = tt[termino, 1], se = tt[termino, 2], p = tt[termino, 4],
         ic_inf = tt[termino, 1] - 1.96 * tt[termino, 2],
         ic_sup = tt[termino, 1] + 1.96 * tt[termino, 2],
         concurvidad_peor = cw, n = nrow(m$model))
}

for (tam_km in tamanos_km) {

  cat("\n\n==================  (B) Riqueza,", tam_km, "km  ==================\n")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  dd   <- datos_alfa(base, pres, tam_km, generos)
  d    <- dd$riqueza$d
  if (is.null(d) || !nrow(d)) { message("sin datos de riqueza a ", tam_km, " km"); next }
  if (!(FOCAL %in% levels(d$genero))) { message(FOCAL, " no entra a ", tam_km, " km"); next }

  y <- dd$riqueza$y; fam <- dd$riqueza$fam; h <- unname(hip)
  trm <- paste0("genero", FOCAL, ":", pred_hist)

  filas <- list()

  # 1. principal: pendiente por género, con suavizado espacial
  m <- ajusta(fml(y, ctrl_riq, h, d), d, fam)
  filas[[1]] <- saca(m, trm, "principal", tam_km)

  # 2. sin término espacial
  m <- ajusta(fml_k(y, ctrl_riq, h, d, divisor = 8, kmax = 40, espacial = FALSE), d, fam)
  filas[[2]] <- saca(m, trm, "sin_espacio", tam_km)

  # 3. suavizado rígido (k = 7 para todos)
  m <- ajusta(fml_k(y, ctrl_riq, h, d, divisor = 1e9, kmax = 7), d, fam)
  filas[[3]] <- saca(m, trm, "suave_rigido_k7", tam_km)

  # 4. esfuerzo medido con los registros del propio género
  if ("esfuerzo_focal_log" %in% names(d)) {
    ctrl_f <- sub("ctrl_esfuerzo_log", "esfuerzo_focal_log", ctrl_riq)
    df <- d %>% filter(!is.na(esfuerzo_focal_log)) %>%
      group_by(genero) %>% mutate(esfuerzo_focal_log = zscore(esfuerzo_focal_log)) %>% ungroup()
    m <- ajusta(fml(y, ctrl_f, h, df), df, fam)
    filas[[4]] <- saca(m, trm, "esfuerzo_del_genero", tam_km)
  }

  # 5. base espacial grande
  m <- ajusta(fml_k(y, ctrl_riq, h, d, divisor = 4, kmax = 100), d, fam)
  filas[[5]] <- saca(m, trm, "base_espacial_grande", tam_km)

  # 6. solo Carollia, sin términos de género
  dc <- d %>% filter(genero == FOCAL) %>% droplevels()
  if (nrow(dc) >= n_min_modelo) {
    kc <- max(5, min(40, floor(nrow(dc) / 8)))
    fc <- as.formula(sprintf("%s ~ %s + s(x_centro_km, y_centro_km, bs = 'tp', k = %d)",
                             y, paste(c(ctrl_riq, h), collapse = " + "), kc))
    mc <- ajusta(fc, dc, fam)
    filas[[6]] <- saca(mc, pred_hist, "solo_carollia", tam_km)
  }

  # 7. desviación respecto de la pendiente común (modelo M1 del script 20)
  m1 <- ajusta(fml(y, c(ctrl_riq, h), character(), d, desv = h), d, fam)
  filas[[7]] <- saca(m1, trm, "desviacion_vs_comun", tam_km)

  tab <- bind_rows(filas)
  print(as.data.frame(tab %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
  riq_todo[[as.character(tam_km)]] <- tab
}

riq_todo <- bind_rows(riq_todo)
write_csv(riq_todo, file.path(dir_mod, "carollia_riqueza_historica.csv"))

cat("\n\n>>> Cómo leer (B):\n",
    "  - 'principal' y 'solo_carollia' deberían coincidir en signo. Si no, el efecto\n",
    "    depende de compartir controles con los otros géneros.\n",
    "  - 'desviacion_vs_comun' es la diferencia respecto de la pendiente compartida,\n",
    "    no la pendiente de Carollia: son cosas distintas y conviene reportar ambas.\n",
    "  - 'concurvidad_peor' manda. Si pasa de 0.9, el coeficiente no está\n",
    "    identificado y no se puede interpretar su magnitud, igual que en la Tabla S6.\n",
    "  - Si el signo cambia entre especificaciones, el hallazgo es ruido y lo mejor\n",
    "    es no reportarlo.\n", sep = "")

cat("\nListo:\n  ", file.path(dir_mod, "carollia_moran_bandas.csv"),
    "\n  ", file.path(dir_mod, "carollia_ses_base_grande.csv"),
    "\n  ", file.path(dir_mod, "carollia_riqueza_historica.csv"), "\n", sep = "")
