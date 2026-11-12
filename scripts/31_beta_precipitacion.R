# =============================================================================
# 31_beta_precipitacion.R
#
# Completa el 2x2 de dimensiones de heterogeneidad en la parte BETA:
#   variable (temperatura / precipitación) x marco (contemporáneo / histórico)
#
# Reutiliza toda la maquinaria del script 21 (recambio de Simpson, GLM
# cuasibinomial, covarianza robusta diádica) y solo cambia los predictores.
#
# Escenarios:
#   lluvia_sin_topo   |Δ CV de bio12| y |Δ hist_sd_prec_log|, controles de 21
#   lluvia_con_topo   lo mismo + |Δ elevación media| como control
#                     (cuánto añade la lluvia por encima del eje topográfico)
#   topo_sd           |Δ SD de elevación| y |Δ hist_sd_temp|
#                     (coherencia con alfa: el modelo principal de beta usa la
#                      diferencia en elevación MEDIA, no en heterogeneidad)
#
# La tabla de pares se construye UNA vez por tamaño de celda y se reutiliza en
# los tres escenarios, que es donde está casi todo el cómputo.
#
# Salidas (datos/modelos):
#   beta_precipitacion.csv
#   beta_precipitacion_generos.csv
#
# Ponlo junto a 21_beta_pruebas.R. EXPLORATORIO: sin corrección por multiplicidad.
# =============================================================================

solo_funciones <- TRUE
source("21_beta_pruebas.R")

# terra (script 30a) enmascara select(); nos aseguramos de usar el de dplyr
select <- dplyr::select

tamanos_km <- c(100, 50, 200)

# ---- variables de celda que hay que arrastrar al par ------------------------
vars_celda <<- c("topo_elev_media", "topo_sd_elev", "hist_sd_temp",
                 "clima_bio01_media", "clima_bio12_media", "ctrl_esfuerzo_log",
                 "clima_bio12_cv", "hist_sd_prec_log")

# ---- construir_pares ampliado (sustituye al del script 21) ------------------
construir_pares <<- function(g, d, X, arboles) {
  n <- nrow(X); esp <- colnames(X)
  Mt <- beta_taxonomico(X)
  Mf <- Reduce(`+`, lapply(arboles, function(tr) {
    tr <- keep.tip(tr, esp); inc <- incidencia(tr)
    beta_filogenetico(X[, tr$tip.label, drop = FALSE], inc$A, inc$len)
  })) / length(arboles)
  ij <- which(lower.tri(Mt), arr.ind = TRUE); i <- ij[, 1]; j <- ij[, 2]
  dif <- function(x) abs(x[i] - x[j])
  lS <- log(d$S)
  tibble(genero = g, ci = as.character(d$celda[i]), cj = as.character(d$celda[j]),
         beta_tax = Mt[ij], beta_filo = Mf[ij],
         dist_km = sqrt((d$x_centro_km[i] - d$x_centro_km[j])^2 +
                        (d$y_centro_km[i] - d$y_centro_km[j])^2),
         d_topo_elev_media  = dif(d$topo_elev_media),
         d_topo_sd_elev     = dif(d$topo_sd_elev),
         d_hist_sd_temp     = dif(d$hist_sd_temp),
         d_bio12_cv         = dif(d$clima_bio12_cv),
         d_hist_sd_prec_log = dif(d$hist_sd_prec_log),
         d_bio01 = dif(d$clima_bio01_media), d_bio12 = dif(d$clima_bio12_media),
         d_esf = dif(d$ctrl_esfuerzo_log),
         d_logS = dif(lS), logS_min = pmin(lS[i], lS[j])) %>%
    filter(!is.na(beta_tax), !is.na(beta_filo))
}

ctrl_base  <- c("d_bio01", "d_bio12", "d_esf", "d_logS", "logS_min")
ctrl_union <- c(ctrl_base, "d_topo_elev_media")

escenarios <- list(
  lluvia_sin_topo = list(hip  = c(actual = "d_bio12_cv", historica = "d_hist_sd_prec_log"),
                         ctrl = ctrl_base),
  lluvia_con_topo = list(hip  = c(actual = "d_bio12_cv", historica = "d_hist_sd_prec_log"),
                         ctrl = ctrl_union),
  topo_sd         = list(hip  = c(presente = "d_topo_sd_elev", historica = "d_hist_sd_temp"),
                         ctrl = ctrl_base),
  # PRUEBA DECISIVA: lo mismo, pero con |Δ elevación media| como control.
  # Si la señal histórica sobrevive aquí, es real; si desaparece, lo que vimos
  # en 'topo_sd' fue solo haber quitado ese control.
  topo_sd_con_elev = list(hip  = c(presente = "d_topo_sd_elev", historica = "d_hist_sd_temp"),
                          ctrl = ctrl_union))

res_todo <- list(); pend_todo <- list()

arboles <- readRDS(archivo_arboles)
arboles <- arboles[seq_len(min(n_arboles_beta, length(arboles)))]

for (tam_km in tamanos_km) {

  cat("\n\n=====================  ", tam_km, " km  =====================\n", sep = "")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)),
                   show_col_types = FALSE) %>%
    mutate(clima_bio12_cv = clima_bio12_sd / clima_bio12_media)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)),
                   show_col_types = FALSE)

  # una sola construcción de pares; se estandariza todo lo que puede usarse
  ctrl <<- c(ctrl_union, "d_topo_sd_elev", "d_hist_sd_temp")
  cat("\nConstruyendo la tabla de pares...\n")
  P <- pares_tamano(base, pres, arboles, tam_km, generos,
                    hipotesis = c(a = "d_bio12_cv", b = "d_hist_sd_prec_log"))
  if (!nrow(P)) { message("Sin pares en ", tam_km, " km; se omite."); next }
  cat("Pares totales:", nrow(P), "\n")

  # dos agrupamientos: bloques de 500 km (inferencia PRIMARIA del manuscrito)
  # y celda (la del script 21, que va al suplementario)
  Pb <- a_bloques(P, base, 500)

  for (nm in names(escenarios)) {
    ctrl <<- escenarios[[nm]]$ctrl
    for (ag in c("bloque", "celda")) {
      cat("\n--", nm, "/", ag, "--\n")
      Pu <- if (ag == "bloque") Pb else P
      r <- pruebas_pares(Pu, tam_km, hipotesis = escenarios[[nm]]$hip)
      res_todo[[paste(tam_km, nm, ag)]]  <- r$pruebas    %>%
        mutate(escenario = nm, agrupamiento = ag)
      pend_todo[[paste(tam_km, nm, ag)]] <- r$pendientes %>%
        mutate(escenario = nm, agrupamiento = ag)
    }
  }
}

res_todo  <- bind_rows(res_todo)
pend_todo <- bind_rows(pend_todo)
write_csv(res_todo,  file.path(dir_mod, "beta_precipitacion.csv"))
write_csv(pend_todo, file.path(dir_mod, "beta_precipitacion_generos.csv"))

cat("\n\n=== PRIMARIA (bloques de 500 km): pendiente común [IC] (p) | interacción p ===\n")
print(as.data.frame(
  res_todo %>% filter(agrupamiento == "bloque") %>%
    mutate(resumen = sprintf("%+.3f [%+.3f, %+.3f] (p=%.3f) | int p=%.4f",
                             est_comun, ic_inf, ic_sup, p_comun, p_int)) %>%
    select(tam_km, escenario, respuesta, tipo, predictor, n_pares, resumen) %>%
    arrange(escenario, respuesta, tipo, tam_km)), row.names = FALSE)

cat("\n\n=== Bloques frente a celda: p de la pendiente común ===\n")
print(as.data.frame(
  res_todo %>%
    select(tam_km, escenario, respuesta, tipo, agrupamiento, est_comun, p_comun) %>%
    tidyr::pivot_wider(names_from = agrupamiento,
                       values_from = c(est_comun, p_comun)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(escenario, respuesta, tipo, tam_km)), row.names = FALSE)

cat("\n=== Pendientes por género (bloques) ===\n")
print(as.data.frame(
  pend_todo %>% filter(agrupamiento == "bloque") %>%
    mutate(pend = sprintf("%+.3f (p=%.3f)", estimate, p)) %>%
    select(tam_km, escenario, respuesta, predictor, genero, pend) %>%
    arrange(escenario, respuesta, predictor, tam_km, genero)), row.names = FALSE)

cat("\nExploratorio: sin corrección por multiplicidad.\n")
