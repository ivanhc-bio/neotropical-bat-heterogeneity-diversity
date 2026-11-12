# 15_robustez_espacial.R
# Fase 3 del plan: robustez de las pruebas confirmatorias al tratamiento del espacio y validación cruzada en bloques espaciales.
# Pruebas confirmatorias (fase 0): H1 topo_sd_elev -> riqueza (pendiente común), H2 interacción género x topo_sd_elev,
# H3 hist_sd_temp -> SES-MPD (pendiente común). Modelo conjunto de los 4 géneros frugívoros (Artibeus, Carollia, Platyrrhinus,
# Sturnira); Anoura se analiza aparte como contraste de H1.
#
# Predictores (modelo principal de 7): hipótesis = topo_sd_elev, clima_bio12_sd (secundaria), hist_sd_temp;
# controles = esfuerzo independiente, área terrestre, temperatura media y precipitación media (en SES: sin área y con log_S).
# Estandarizados dentro de cada género, así que una pendiente es "por DE dentro del género".
#
# Tratamientos del espacio (columna `tratamiento`):
#   ninguno   sin término espacial
#   rigido    suavizado espacial s(x, y) por género con k fijo = k_rigido (curva rígida: solo estructura regional)
#   actual    s(x, y) por género con k = n / 8 (entre 5 y 40), como el script 11
#   bloques   efecto aleatorio de bloques de bloque_km x bloque_km (estructura regional sin suavizado fino)
# Para cada uno: pendiente común y pendientes por género de cada predictor de hipótesis, prueba global de H2 (razón de
# verosimilitud con ML; si falla, solo dAIC) y el contraste de Anoura (pendiente de Anoura frente a la común de los frugívoros).
# Holm: se calcula sobre las 3 pruebas confirmatorias del tratamiento "actual" (p de H1, H2 y H3).
#
# Validación cruzada en bloques espaciales (sin término espacial: un suavizado no puede predecir bloques no vistos): se divide
# el área en bloques de bloque_km, se asignan al azar a n_folds grupos y se compara la devianza en los bloques retenidos del
# modelo solo con controles frente al modelo con los 3 predictores de hipótesis. Ganancia = 1 - devianza_completo / devianza_controles.
#
# Salidas (datos/modelos): robustez_espacial_{tam}km.csv, robustez_h2_{tam}km.csv, robustez_contraste_anoura_{tam}km.csv,
#                          validacion_bloques_{tam}km.csv
# No probado: no tengo R en mi entorno. Prueba con tamanos_km <- 100 y revisa lo que imprime.
# Requiere: mgcv, dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(mgcv)
library(dplyr)
library(readr)
library(tidyr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")

tamanos_km   <- c(50, 200)
frugivoros   <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
contraste    <- "Anoura"
n_min_modelo <- 25
k_por_celdas <- 8
k_rigido     <- 7
bloque_km    <- 500
n_folds      <- 5
semilla      <- 1
tratamientos <- c("ninguno", "rigido", "actual", "bloques")

hip      <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp")
ctrl_riq <- c("ctrl_esfuerzo_log", "ctrl_area_tierra_km2", "clima_bio01_media", "clima_bio12_media")
ctrl_ses <- c("ctrl_esfuerzo_log", "clima_bio01_media", "clima_bio12_media", "log_S")
respuestas <- list(
  riqueza = list(y = "cbind(S, S_fallo)", familia = "binomial", ctrl = ctrl_riq, solo_filo = FALSE, h_conf = "topo_sd_elev", h_nombre = "H1"),
  ses_mpd = list(y = "ses_mpd", familia = "gaussian", ctrl = ctrl_ses, solo_filo = TRUE, h_conf = "hist_sd_temp", h_nombre = "H3"))

# ---- 1. Funciones -------------------------------------------------------------------
zscore <- function(x) as.numeric(scale(x))

preparar <- function(d, vars) {
  gens <- sort(unique(d$genero))
  d <- d %>% group_by(genero) %>% mutate(across(all_of(vars), zscore)) %>% ungroup()
  d$bloque <- factor(paste(floor(d$x_centro_km / bloque_km), floor(d$y_centro_km / bloque_km)))
  for (g in gens) d[[paste0("I_", g)]] <- as.numeric(d$genero == g)
  d$genero <- factor(d$genero)
  d
}

k_por_genero <- function(d) {
  n <- table(d$genero)
  setNames(pmax(5, pmin(40, floor(as.numeric(n) / k_por_celdas))), names(n))
}

terminos_espacio <- function(trat, d) {
  gens <- levels(d$genero)
  switch(trat,
    ninguno = character(),
    rigido  = sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", gens, k_rigido),
    actual  = { k <- k_por_genero(d); sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k), as.integer(k)) },
    bloques = "s(bloque, bs = 're')")
}

formula_conjunta <- function(y, trat, d, ctrl, separados, comunes_hip) {
  partes <- c("genero", ctrl, comunes_hip, if (length(separados)) paste0("genero:", separados), terminos_espacio(trat, d))
  as.formula(paste(y, "~", paste(partes, collapse = " + ")))
}

ajusta <- function(f, d, fam, metodo) {
  tryCatch(gam(f, family = fam, data = d, method = metodo), error = function(e) NULL)
}

# ---- 2. Por tamaño de celda ----------------------------------------------------------------
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  base <- bind_rows(lapply(c(frugivoros, contraste), function(g) {
    d <- base %>% filter(genero == g)
    if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    d
  }))

  coefs <- list(); h2 <- list(); anou <- list(); cvs <- list()
  for (rn in names(respuestas)) {
    R <- respuestas[[rn]]
    vars <- c(hip, R$ctrl)
    d0 <- base %>% filter(genero %in% c(frugivoros, contraste))
    if (R$solo_filo) d0 <- d0 %>% filter(filo_valida, !is.na(ses_mpd))
    d0 <- d0 %>% filter(if_all(all_of(vars), ~ !is.na(.x)))
    gens_ok <- names(which(table(d0$genero) >= n_min_modelo))
    dfr <- preparar(d0 %>% filter(genero %in% intersect(frugivoros, gens_ok)), vars)
    dan <- if (contraste %in% gens_ok) preparar(d0 %>% filter(genero == contraste), vars) else NULL
    fam <- if (R$familia == "binomial") binomial() else gaussian()
    cat(sprintf("\n[%s] frugívoros: %s (n = %d)%s\n", rn, paste(levels(dfr$genero), collapse = ", "), nrow(dfr),
                if (is.null(dan)) "; Anoura omitida (pocas celdas)" else sprintf("; Anoura n = %d", nrow(dan))))

    for (trat in tratamientos) {
      # pendiente común y por género de los predictores de hipótesis (REML)
      m_com <- ajusta(formula_conjunta(R$y, trat, dfr, R$ctrl, character(), hip), dfr, fam, "REML")
      m_gen <- ajusta(formula_conjunta(R$y, trat, dfr, R$ctrl, hip, character()), dfr, fam, "REML")
      if (!is.null(m_com)) {
        tc <- as.data.frame(summary(m_com)$p.table); names(tc) <- c("estimate", "se", "stat", "p"); tc$termino <- rownames(tc)
        coefs[[paste(rn, trat, "comun")]] <- tc %>% filter(termino %in% hip) %>%
          mutate(respuesta = rn, tratamiento = trat, modelo = "comun", genero = "todos")
      }
      if (!is.null(m_gen)) {
        tg <- as.data.frame(summary(m_gen)$p.table); names(tg) <- c("estimate", "se", "stat", "p"); tg$nombre <- rownames(tg)
        coefs[[paste(rn, trat, "genero")]] <- tg %>% filter(grepl("^genero.*:", nombre)) %>%
          mutate(genero = sub("^genero([^:]+):.*$", "\\1", nombre), termino = sub("^genero[^:]+:", "", nombre),
                 respuesta = rn, tratamiento = trat, modelo = "por_genero") %>% select(-nombre)
      }
      # H2: interacción genero x predictor de la hipótesis confirmatoria (ML)
      f_full <- formula_conjunta(R$y, trat, dfr, R$ctrl, R$h_conf, setdiff(hip, R$h_conf))
      f_nulo <- formula_conjunta(R$y, trat, dfr, R$ctrl, character(), hip)
      mf <- ajusta(f_full, dfr, fam, "ML"); mn <- ajusta(f_nulo, dfr, fam, "ML")
      if (!is.null(mf) && !is.null(mn)) {
        a <- tryCatch(anova(mn, mf, test = "Chisq"), error = function(e) NULL)
        dv <- if (is.null(a)) NA_real_ else a$Deviance[2]; pv <- if (is.null(a)) NA_real_ else a[["Pr(>Chi)"]][2]
        ok <- !is.na(dv) && dv > 0
        # prueba de Wald conjunta de las interacciones genero:predictor (REML; valida aunque la LR falle)
        mw <- ajusta(formula_conjunta(R$y, trat, dfr, R$ctrl, R$h_conf, hip), dfr, fam, "REML")
        pw <- NA_real_; glw <- NA_real_
        if (!is.null(mw)) {
          at <- tryCatch(anova(mw)$pTerms.table, error = function(e) NULL)
          fila <- if (is.null(at)) integer() else which(rownames(at) == paste0("genero:", R$h_conf))
          if (length(fila) == 1) { pw <- at[fila, "p-value"]; glw <- at[fila, "df"] }
        }
        h2[[paste(rn, trat)]] <- tibble(respuesta = rn, tratamiento = trat, predictor = R$h_conf, gl = if (is.null(a)) NA_real_ else a$Df[2],
                                        devianza = dv, p_LR = if (ok) pv else NA_real_, gl_wald = glw, p_wald = pw,
                                        dAIC = AIC(mf) - AIC(mn),
                                        nota = if (ok) "" else "LR no valida; usar p_wald y dAIC")
      }
      # contraste de Anoura: pendiente de Anoura frente a la común de los frugívoros
      if (!is.null(dan) && !is.null(m_com)) {
        t_an <- switch(trat, ninguno = character(),
          rigido = sprintf("s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", k_rigido),
          actual = sprintf("s(x_centro_km, y_centro_km, bs = 'tp', k = %d)", max(5, min(40, floor(nrow(dan) / k_por_celdas)))),
          bloques = "s(bloque, bs = 're')")
        f_an <- as.formula(paste(R$y, "~", paste(c(R$ctrl, hip, t_an), collapse = " + ")))
        m_an <- ajusta(f_an, dan, fam, "REML")
        if (!is.null(m_an)) {
          ta <- as.data.frame(summary(m_an)$p.table); names(ta) <- c("estimate", "se", "stat", "p"); ta$termino <- rownames(ta)
          ta <- ta %>% filter(termino %in% hip)
          tc <- coefs[[paste(rn, trat, "comun")]]
          anou[[paste(rn, trat)]] <- ta %>% select(termino, est_anoura = estimate, se_anoura = se) %>%
            left_join(tc %>% select(termino, est_comun = estimate, se_comun = se), by = "termino") %>%
            mutate(z = (est_anoura - est_comun) / sqrt(se_anoura^2 + se_comun^2), p_contraste = 2 * pnorm(-abs(z)),
                   respuesta = rn, tratamiento = trat)
        }
      }
    }

    # ---- validación cruzada en bloques espaciales (sin término espacial) ----
    set.seed(semilla)
    bl <- levels(dfr$bloque); fold_b <- setNames(sample(rep(seq_len(n_folds), length.out = length(bl))), bl)
    dfr$fold <- fold_b[as.character(dfr$bloque)]
    f_c <- as.formula(paste(R$y, "~ genero +", paste(R$ctrl, collapse = " + ")))
    f_f <- as.formula(paste(R$y, "~ genero +", paste(c(R$ctrl, hip), collapse = " + ")))
    devianza <- function(m, te) {
      mu <- as.numeric(predict(m, te, type = "response"))
      if (R$familia == "binomial") { n <- te$S + te$S_fallo; sum(binomial()$dev.resids(te$S / n, mu, n)) }
      else sum((te$ses_mpd - mu)^2)
    }
    cv <- lapply(seq_len(n_folds), function(k) {
      tr <- dfr[dfr$fold != k, ]; te <- dfr[dfr$fold == k, ]
      if (nrow(te) < 5) return(NULL)
      m1 <- ajusta(f_c, tr, fam, "REML"); m2 <- ajusta(f_f, tr, fam, "REML")
      if (is.null(m1) || is.null(m2)) return(NULL)
      tibble(fold = k, n_test = nrow(te), dev_controles = devianza(m1, te), dev_completo = devianza(m2, te))
    })
    cvs[[rn]] <- bind_rows(cv) %>% mutate(respuesta = rn, ganancia = 1 - dev_completo / dev_controles)
  }

  coefs <- bind_rows(coefs) %>% mutate(ic_inf = estimate - 1.96 * se, ic_sup = estimate + 1.96 * se, tam_km = tam_km) %>%
    select(tam_km, respuesta, tratamiento, modelo, genero, termino, estimate, se, p, ic_inf, ic_sup)
  h2 <- bind_rows(h2) %>% mutate(tam_km = tam_km); anou <- bind_rows(anou) %>% mutate(tam_km = tam_km)
  cvs <- bind_rows(cvs) %>% mutate(tam_km = tam_km)
  write_csv(coefs, file.path(dir_mod, sprintf("robustez_espacial_%dkm.csv", tam_km)))
  write_csv(h2, file.path(dir_mod, sprintf("robustez_h2_%dkm.csv", tam_km)))
  write_csv(anou, file.path(dir_mod, sprintf("robustez_contraste_anoura_%dkm.csv", tam_km)))
  write_csv(cvs, file.path(dir_mod, sprintf("validacion_bloques_%dkm.csv", tam_km)))

  # ---- 3. Resumen ----------------------------------------------------------------------------
  r3 <- function(x) round(x, 3)
  for (rn in names(respuestas)) {
    R <- respuestas[[rn]]
    cat(sprintf("\n--- %s: pendiente común de los predictores de hipótesis (por DE; ¿se mantiene al cambiar el espacio?) ---\n", rn))
    print(as.data.frame(coefs %>% filter(respuesta == rn, modelo == "comun") %>%
            mutate(estimado = sprintf("%.3f (%.3f)", estimate, se), p = r3(p)) %>% select(tratamiento, termino, estimado, p)), row.names = FALSE)
    cat(sprintf("\n--- %s: pendientes por género de %s ---\n", rn, R$h_conf))
    print(as.data.frame(coefs %>% filter(respuesta == rn, modelo == "por_genero", termino == R$h_conf) %>%
            mutate(estimado = sprintf("%.2f (%.2f)", estimate, se)) %>% select(tratamiento, genero, estimado, p) %>%
            mutate(p = r3(p)) %>% pivot_wider(names_from = genero, values_from = c(estimado, p))), row.names = FALSE)
  }
  cat("\n--- H2: prueba global de interacción genero x predictor (p pequeña = el efecto difiere entre géneros) ---\n")
  print(as.data.frame(h2 %>% mutate(across(where(is.numeric), r3)) %>% select(respuesta, tratamiento, predictor, p_LR, gl_wald, p_wald, dAIC, nota)), row.names = FALSE)
  cat("\n--- Contraste de Anoura frente a la pendiente común de los frugívoros ---\n")
  print(as.data.frame(anou %>% filter(termino %in% c("topo_sd_elev", "hist_sd_temp")) %>%
          mutate(across(where(is.numeric), r3)) %>% select(respuesta, tratamiento, termino, est_anoura, est_comun, p_contraste)), row.names = FALSE)
  cat("\n--- Validación cruzada en bloques (ganancia > 0: los predictores de hipótesis mejoran la predicción en bloques no vistos) ---\n")
  print(as.data.frame(cvs %>% group_by(respuesta) %>%
          summarise(folds = n(), ganancia_media = r3(mean(ganancia)), folds_con_ganancia = sum(ganancia > 0), .groups = "drop")), row.names = FALSE)

  # Holm sobre las 3 pruebas confirmatorias (tratamiento "actual")
  conf <- bind_rows(
    coefs %>% filter(tratamiento == "actual", modelo == "comun", respuesta == "riqueza", termino == "topo_sd_elev") %>% transmute(prueba = "H1", p),
    h2 %>% filter(tratamiento == "actual", respuesta == "riqueza") %>% transmute(prueba = "H2", p = p_wald),
    coefs %>% filter(tratamiento == "actual", modelo == "comun", respuesta == "ses_mpd", termino == "hist_sd_temp") %>% transmute(prueba = "H3", p))
  conf$p_holm <- p.adjust(conf$p, method = "holm")
  cat("\n--- Pruebas confirmatorias (tratamiento actual), con corrección de Holm ---\n")
  print(as.data.frame(conf %>% mutate(across(where(is.numeric), r3))), row.names = FALSE)
}
