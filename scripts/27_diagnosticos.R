# =============================================================================
# 27_diagnosticos.R
# Diagnósticos de los modelos alfa principales, para la tabla suplementaria:
#   - concurvidad del suavizado espacial de cada género (mgcv::concurvity)
#   - I de Moran de los residuos, por género, con pesos de distancia inversa
#   - VIF de los predictores paramétricos
#   - devianza explicada y dispersión de cada modelo
#
# Salidas (datos/modelos):
#   diagnosticos_alfa.csv    una fila por tamaño, respuesta y género
#   diagnosticos_vif.csv     una fila por tamaño, respuesta y predictor
#
# Usa exactamente los mismos datos y la misma fórmula que el script 20.
# Requiere: mgcv, dplyr, readr, tibble y ape (para la I de Moran).
# =============================================================================

solo_funciones <- TRUE
source("20_alfa_pruebas.R")     # funciones, rutas y parámetros

suppressPackageStartupMessages({library(tibble)})
if (!requireNamespace("ape", quietly = TRUE))
  stop("Falta el paquete ape. Instálalo con install.packages('ape').")

tamanos_km <- c(50, 100, 200)

# ---- datos de una escala, igual que alfa_tamano ------------------------------

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

# ---- I de Moran por género, pesos de distancia inversa -----------------------

moran_generos <- function(d, resid) {
  bind_rows(lapply(levels(d$genero), function(g) {
    i <- which(d$genero == g)
    if (length(i) < 10) return(tibble(genero = g, n_g = length(i),
                                      moran_I = NA_real_, moran_p = NA_real_))
    xy <- as.matrix(d[i, c("x_centro_km", "y_centro_km")])
    w  <- 1 / as.matrix(dist(xy))
    diag(w) <- 0
    w[!is.finite(w)] <- 0
    mi <- tryCatch(ape::Moran.I(resid[i], w, scaled = TRUE),
                   error = function(e) list(observed = NA_real_, p.value = NA_real_))
    tibble(genero = g, n_g = length(i), moran_I = mi$observed, moran_p = mi$p.value)
  }))
}

# ---- VIF simple sobre los predictores paramétricos ---------------------------

vif_simple <- function(d, vars) {
  vapply(vars, function(v) {
    otras <- setdiff(vars, v)
    r2 <- summary(lm(as.formula(paste(v, "~", paste(otras, collapse = " + "))), data = d))$r.squared
    1 / (1 - r2)
  }, numeric(1))
}

# ---- bucle -------------------------------------------------------------------

diag_todo <- list(); vif_todo <- list()

for (tam_km in tamanos_km) {

  cat("\n\n=====================  ", tam_km, " km  =====================\n", sep = "")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  dd <- datos_alfa(base, pres, tam_km, generos)

  for (resp in names(dd)) {

    d <- dd[[resp]]$d
    if (is.null(d) || !nrow(d)) { cat("\n--", resp, ": sin datos\n"); next }
    cat("\n--", resp, ": n =", nrow(d), "\n")

    m0 <- ajusta(fml(dd[[resp]]$y, c(dd[[resp]]$ctrl, unname(hip)), character(), d),
                 d, dd[[resp]]$fam)
    s0 <- summary(m0)

    # concurvidad: una columna por término, filas worst / observed / estimate
    cc <- concurvity(m0, full = TRUE)
    cc <- as.data.frame(t(cc)) %>% rownames_to_column("termino")
    names(cc) <- make.names(names(cc))
    espaciales <- cc %>% filter(grepl("^s\\(x_centro_km", termino))
    espaciales$genero <- sub(".*I_", "", espaciales$termino)

    mo <- moran_generos(d, residuals(m0, type = "deviance"))

    fila <- espaciales %>%
      select(genero, concurvidad_peor = worst, concurvidad_obs = observed) %>%
      full_join(mo, by = "genero") %>%
      mutate(tam_km = tam_km, respuesta = resp, n_total = nrow(d),
             dev_expl = s0$dev.expl, dispersion = s0$scale) %>%
      select(tam_km, respuesta, genero, n_g, n_total, concurvidad_peor,
             concurvidad_obs, moran_I, moran_p, dev_expl, dispersion)

    print(as.data.frame(fila %>% mutate(across(where(is.numeric), ~ round(.x, 4)))),
          row.names = FALSE)
    diag_todo[[paste(tam_km, resp)]] <- fila

    vv <- c(unname(hip), dd[[resp]]$ctrl)
    vv <- vv[vapply(vv, function(v) is.numeric(d[[v]]) && var(d[[v]], na.rm = TRUE) > 0, TRUE)]
    v  <- vif_simple(d, vv)
    cat("\n   VIF:", paste(sprintf("%s = %.2f", names(v), v), collapse = "; "), "\n")
    vif_todo[[paste(tam_km, resp)]] <-
      tibble(tam_km = tam_km, respuesta = resp, predictor = names(v), vif = as.numeric(v))
  }
}

diag_todo <- bind_rows(diag_todo); vif_todo <- bind_rows(vif_todo)
write_csv(diag_todo, file.path(dir_mod, "diagnosticos_alfa.csv"))
write_csv(vif_todo,  file.path(dir_mod, "diagnosticos_vif.csv"))

cat("\n\nListo:\n  ", file.path(dir_mod, "diagnosticos_alfa.csv"),
    "\n  ", file.path(dir_mod, "diagnosticos_vif.csv"), "\n", sep = "")
