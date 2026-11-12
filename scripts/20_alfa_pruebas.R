# 20_alfa_pruebas.R
# Diseño definitivo, parte alfa. Pruebas formales de heterogeneidad presente (topo_sd_elev) e histórica (hist_sd_temp)
# contra la riqueza y el SES-MPD de los 4 géneros frugívoros, en un modelo conjunto (GAM, REML).
#
# Por cada respuesta (riqueza, SES-MPD) y cada predictor de hipótesis se obtienen 2 pruebas:
#   (a) pendiente COMÚN a los géneros: Wald del predictor en el modelo M0 (ambos predictores con pendiente común, sin interacción).
#   (b) INTERACCIÓN género x predictor: Wald de 3 gl en el modelo M1 (pendiente común + desviación por género de AMBOS predictores).
# Son 2 predictores x 2 respuestas x 2 pruebas = 8 pruebas por tamaño de celda; las otras 8 vienen del script de beta.
# Holm sobre las 16 pruebas del tamaño de celda: lo hace el script 22 (junta alfa y beta).
#
# Modelo (todas las respuestas): género + controles comunes + 2 predictores de hipótesis + s(x, y) por género (thin-plate,
#   k = n_género / 8, entre 5 y 40). Predictores estandarizados dentro de género y tamaño de celda.
#   riqueza : binomial cbind(S, S_fallo); controles = esfuerzo (Mammalia, log), área de tierra, temperatura media, precipitación media.
#   SES-MPD : gaussiano; controles = esfuerzo, temperatura media, precipitación media, log_S (sin área, como antes).
# Nota: en M1 el coeficiente "topo_sd_elev" es la pendiente del género de referencia, no la común; para inferir la pendiente
# común se usa M0. Las pendientes por género (con IC) salen de M_sep (predictor solo por género).
#
# Salidas (datos/modelos): alfa_pruebas_{tam}km.csv (8 filas) y alfa_pendientes_genero_{tam}km.csv
# No probado: no hay R en mi entorno. Prueba con tamanos_km <- 100.
# Requiere: mgcv, dplyr, readr
# Para reutilizar las funciones (robustez): solo_funciones <- TRUE; source("20_alfa_pruebas.R")

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv); library(dplyr); library(readr)

dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")

tamanos_km   <- c(100, 50, 200)
generos      <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
hip          <- c(presente = "topo_sd_elev", historica = "hist_sd_temp")
ctrl_riq     <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2", "clima_bio01_media", "clima_bio12_media")
ctrl_ses     <- c("ctrl_esfuerzo_log", "clima_bio01_media", "clima_bio12_media", "log_S")
n_min_modelo <- 25
k_por_celdas <- 8

zscore <- function(x) as.numeric(scale(x))
preparar <- function(d, vars) {
  gens <- sort(unique(d$genero))
  d <- d %>% group_by(genero) %>% mutate(across(all_of(vars), zscore)) %>% ungroup()
  for (g in gens) d[[paste0("I_", g)]] <- as.numeric(d$genero == g)
  d$genero <- factor(d$genero); d
}
k_por_genero <- function(d) { n <- table(d$genero); setNames(pmax(5, pmin(40, floor(as.numeric(n) / k_por_celdas))), names(n)) }
suaves <- function(d) { k <- k_por_genero(d)
  sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k), as.integer(k)) }
fml <- function(y, comunes, separados, d, desv = character()) {  # desv: desviación por género ADEMÁS de la pendiente común
  sep <- c(separados, desv)
  # ojo: paste0("genero:", character(0)) devuelve "genero:", por eso se protege
  as.formula(paste(y, "~", paste(c("genero", comunes, if (length(sep)) paste0("genero:", sep), suaves(d)), collapse = " + ")))
}
ajusta <- function(f, d, fam) tryCatch(gam(f, family = fam, data = d, method = "REML"),
                                       error = function(e) gam(f, family = fam, data = d, method = "GCV.Cp"))
wald <- function(m, termino) {            # Wald de un término paramétrico (gl, p)
  at <- tryCatch(anova(m)$pTerms.table, error = function(e) NULL)
  if (is.null(at) || !(termino %in% rownames(at))) return(c(gl = NA_real_, p = NA_real_))
  c(gl = at[termino, "df"], p = at[termino, "p-value"])
}

# Una respuesta: devuelve (pruebas, pendientes por género)
alfa_respuesta <- function(d, y, fam, ctrl, etiqueta, tam_km, hipotesis = hip) {
  h <- unname(hipotesis)
  m0 <- ajusta(fml(y, c(ctrl, h), character(), d), d, fam)                 # pendientes comunes
  m1 <- ajusta(fml(y, c(ctrl, h), character(), d, desv = h), d, fam)       # común + desviación por género
  ms <- ajusta(fml(y, ctrl, h, d), d, fam)                                 # pendiente por género (para IC por género)
  tab0 <- as.data.frame(summary(m0)$p.table)
  pruebas <- bind_rows(lapply(names(hipotesis), function(tipo) {
    p <- hipotesis[[tipo]]
    wi <- wald(m1, paste0("genero:", p))
    tibble(tam_km = tam_km, respuesta = etiqueta, tipo = tipo, predictor = p, n = nrow(d),
           generos = paste(levels(d$genero), collapse = ","),
           est_comun = tab0[p, 1], se_comun = tab0[p, 2], ic_inf = tab0[p, 1] - 1.96 * tab0[p, 2],
           ic_sup = tab0[p, 1] + 1.96 * tab0[p, 2], p_comun = tab0[p, 4],
           gl_int = wi[["gl"]], p_int = wi[["p"]])
  }))
  ts <- as.data.frame(summary(ms)$p.table); ts$term <- rownames(ts); names(ts)[1:4] <- c("estimate", "se", "stat", "p")
  pend <- ts %>% filter(grepl("^genero[^:]+:", term)) %>%
    mutate(genero = sub("^genero([^:]+):.*$", "\\1", term), predictor = sub("^genero[^:]+:", "", term),
           ic_inf = estimate - 1.96 * se, ic_sup = estimate + 1.96 * se, respuesta = etiqueta, tam_km = tam_km) %>%
    filter(predictor %in% h) %>% select(tam_km, respuesta, genero, predictor, estimate, se, p, ic_inf, ic_sup)
  list(pruebas = pruebas, pendientes = pend)
}

# Prepara los datos de un tamaño de celda y corre las dos respuestas
alfa_tamano <- function(base, pres, tam_km, generos, hipotesis = hip, ctrl_r = ctrl_riq, ctrl_s = ctrl_ses) {
  h <- unname(hipotesis)
  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g); if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    if (any(d$S_fallo < 0)) stop("Celdas con más especies que el grupo de especies: ", g)
    d
  }))
  out <- list()
  vr <- c(h, ctrl_r)
  d0 <- base %>% filter(if_all(all_of(vr), ~ !is.na(.x)))
  ok <- names(which(table(d0$genero) >= n_min_modelo))
  dr <- preparar(d0 %>% filter(genero %in% ok), vr)
  cat(sprintf("  riqueza: %s; n = %d\n", paste(ok, collapse = ","), nrow(dr)))
  out$riq <- alfa_respuesta(dr, "cbind(S, S_fallo)", binomial(), ctrl_r, "riqueza", tam_km, hipotesis)
  vs <- c(h, ctrl_s)
  d1 <- base %>% filter(filo_valida, !is.na(ses_mpd), if_all(all_of(vs), ~ !is.na(.x)))
  ok <- names(which(table(d1$genero) >= n_min_modelo))
  ds <- preparar(d1 %>% filter(genero %in% ok), vs)
  cat(sprintf("  SES-MPD: %s; n = %d\n", paste(ok, collapse = ","), nrow(ds)))
  out$ses <- alfa_respuesta(ds, "ses_mpd", gaussian(), ctrl_s, "ses_mpd", tam_km, hipotesis)
  list(pruebas = bind_rows(out$riq$pruebas, out$ses$pruebas), pendientes = bind_rows(out$riq$pendientes, out$ses$pendientes))
}

if (!exists("solo_funciones") || !solo_funciones) {
  for (tam_km in tamanos_km) {
    cat("\n=== Celdas de", tam_km, "km ===\n")
    base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
    pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
    r <- alfa_tamano(base, pres, tam_km, generos)
    write_csv(r$pruebas, file.path(dir_mod, sprintf("alfa_pruebas_%dkm.csv", tam_km)))
    write_csv(r$pendientes, file.path(dir_mod, sprintf("alfa_pendientes_genero_%dkm.csv", tam_km)))
    cat("\nPruebas (est_comun = pendiente común por DE; p_comun = Wald; p_int = Wald de la interacción, 3 gl):\n")
    print(as.data.frame(r$pruebas %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
                          select(respuesta, tipo, predictor, n, est_comun, ic_inf, ic_sup, p_comun, gl_int, p_int)), row.names = FALSE)
  }
}
