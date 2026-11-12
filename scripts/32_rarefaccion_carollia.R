# =============================================================================
# 32_rarefaccion_carollia.R
#
# PREGUNTA: Carollia se comporta distinto de los otros tres generos en varios
# analisis. Es biologia o es que su riqueza tiene techo?
#   Carollia:     S max por celda = 5, S mediana = 2, pool regional 8 especies
#   los otros:    S max 8-10, pool 18-24
# Con 2 especies en la celda mediana, el recambio de Simpson solo puede tomar
# unos pocos valores (0, 0.5, 1). Esa resolucion mas gruesa basta para explicar
# que Carollia salga distinto, sin necesidad de ninguna hipotesis ecologica.
#
# PRUEBA: rarificar los otros tres generos al techo de Carollia (5 especies por
# celda, al azar) y ver si sus pendientes se mueven HACIA la de Carollia.
#   - se mueven  -> artefacto de resolucion; no hay nada que interpretar
#   - no se mueven -> queda algo biologico y merece un parrafo propio
# Carollia nunca pasa de 5, asi que la rarefaccion no lo toca: es el control interno.
#
# La rarefaccion se hace sobre la tabla de presencia ANTES de llamar a
# pares_tamano, asi que no se toca nada del script 21.
#
# Salidas (datos/modelos):
#   rarefaccion_pruebas.csv     pendiente comun y p de interaccion
#   rarefaccion_generos.csv     pendientes por genero  <- ESTA es la que importa
#
# Ponlo junto a 21_beta_pruebas.R. EXPLORATORIO: sin correccion por multiplicidad.
# TARDA: n_sorteos x 2 escenarios x 2 respuestas. Si va lento baja n_arboles_beta
# a 5 y n_sorteos a 5; para una prueba rapida usa tamanos_km <- 100 solamente.
# =============================================================================

solo_funciones <- TRUE
source("21_beta_pruebas.R")
select <- dplyr::select

tamanos_km <- c(100, 50, 200)   # 100 primero: ahi no hay submuestreo de celdas
techo      <- 5                 # S max de Carollia
referente  <- "Carollia"        # genero que NO se rarifica: es el control interno
n_sorteos  <- 10                # repeticiones de la rarefaccion

# ---- variables de celda y pares: igual que en el script 31 ------------------
vars_celda <<- c("topo_elev_media", "topo_sd_elev", "hist_sd_temp",
                 "clima_bio01_media", "clima_bio12_media", "ctrl_esfuerzo_log",
                 "clima_bio12_cv", "hist_sd_prec_log")

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

# ---- la rarefaccion ---------------------------------------------------------
# Deja como mucho `techo` especies por celda y genero, elegidas al azar.
# Las celdas que ya tienen <= techo quedan intactas (slice_sample devuelve el grupo entero).
# (a) TECHO: deja como mucho `techo` especies por celda. Intervencion suave:
#     a 100 km solo recorta 83 celdas de Artibeus y ~12-15 de los otros dos.
rarificar_techo <- function(pres, techo, ref) {
  tab <- pres %>% distinct(genero_analisis, celda, especie)
  bind_rows(
    tab %>% filter(genero_analisis == ref),                     # el referente no se toca
    tab %>% filter(genero_analisis != ref) %>%
      group_by(genero_analisis, celda) %>% slice_sample(n = techo) %>% ungroup())
}

# (b) DISTRIBUCION: a cada celda se le sortea un S objetivo de la distribucion de S
#     del genero de referencia (entre sus celdas con S >= s_min) y se recorta a ese
#     numero. Iguala la mediana y no solo el maximo, asi que es la prueba con fuerza.
rarificar_dist <- function(pres, ref) {
  tab  <- pres %>% distinct(genero_analisis, celda, especie)
  Sref <- tab %>% filter(genero_analisis == ref) %>% count(celda, name = "S") %>%
            filter(S >= s_min) %>% pull(S)
  obj  <- tab %>% filter(genero_analisis != ref) %>% count(genero_analisis, celda, name = "S") %>%
            mutate(objetivo = pmin(S, sample(Sref, n(), replace = TRUE)))
  bind_rows(
    tab %>% filter(genero_analisis == ref),
    tab %>% filter(genero_analisis != ref) %>%
      left_join(obj, by = c("genero_analisis", "celda")) %>%
      group_by(genero_analisis, celda) %>%
      mutate(.orden = sample(n())) %>%        # slice_sample() exige n constante: se baraja y se corta
      filter(.orden <= objetivo) %>%
      ungroup() %>% select(genero_analisis, celda, especie))
}

ctrl_base <- c("d_bio01", "d_bio12", "d_esf", "d_logS", "logS_min")
escenarios <- list(
  lluvia = list(hip = c(actual = "d_bio12_cv", historica = "d_hist_sd_prec_log"), ctrl = ctrl_base),
  topo   = list(hip = c(presente = "d_topo_sd_elev", historica = "d_hist_sd_temp"), ctrl = ctrl_base))

arboles <- readRDS(archivo_arboles)
arboles <- arboles[seq_len(min(n_arboles_beta, length(arboles)))]

res_todo <- list(); pend_todo <- list()

for (tam_km in tamanos_km) {
  cat("\n\n=====================  ", tam_km, " km  =====================\n", sep = "")

  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)),
                   show_col_types = FALSE) %>%
    mutate(clima_bio12_cv = clima_bio12_sd / clima_bio12_media)
  pres0 <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)),
                    show_col_types = FALSE)

  # cuanto recorta el techo en cada genero (control de que la prueba hace algo)
  cat("\nCeldas por genero con S > ", techo, " (las que el techo recorta):\n", sep = "")
  print(pres0 %>% distinct(genero_analisis, celda, especie) %>%
          count(genero_analisis, celda, name = "S") %>%
          group_by(genero_analisis) %>%
          summarise(celdas = n(), S_max = max(S), S_mediana = median(S),
                    recortadas = sum(S > techo), .groups = "drop") %>%
          as.data.frame(), row.names = FALSE)

  # tipo "completo" = referencia sin tocar; "techo" y "dist" = las dos rarefacciones
  corridas <- c(list(list(tipo = "completo", s = 0)),
                lapply(1:n_sorteos, function(s) list(tipo = "techo", s = s)),
                lapply(1:n_sorteos, function(s) list(tipo = "dist",  s = s)))

  for (cc in corridas) {
    etiqueta <- if (cc$tipo == "completo") "completo" else sprintf("%s_s%02d", cc$tipo, cc$s)
    cat("\n--", etiqueta, "--\n")

    # semilla propia de la rarefaccion; pares_tamano re-siembra despues, asi que el
    # submuestreo de celdas es el mismo en todos los sorteos y solo varia la rarefaccion
    set.seed(1000 + cc$s)
    pres <- switch(cc$tipo,
                   completo = pres0,
                   techo    = rarificar_techo(pres0, techo, referente),
                   dist     = rarificar_dist(pres0, referente))

    # Diagnostico que SI es sensible: la S por celda despues de rarificar.
    # (celdas, especies del genero y pares no cambian aunque la rarefaccion funcione.)
    cat("   S por celda (celdas con S >= ", s_min, "): ", sep = "")
    cat(paste(pres %>% distinct(genero_analisis, celda, especie) %>%
                count(genero_analisis, celda, name = "S") %>% filter(S >= s_min) %>%
                group_by(genero_analisis) %>%
                summarise(t = sprintf("%s media %.2f max %d", genero_analisis[1], mean(S), max(S)),
                          .groups = "drop") %>% pull(t), collapse = " | "), "\n")

    # OJO: ctrl NO debe repetir las variables que van en `hipotesis`: el script 21 hace
    # est <- c(unname(hipotesis), ctrl) y luego pr[est] <- ..., que falla con indices duplicados.
    ctrl <<- c(ctrl_base, "d_topo_sd_elev", "d_hist_sd_temp")
    P <- pares_tamano(base, pres, arboles, tam_km, generos,
                      hipotesis = c(a = "d_bio12_cv", b = "d_hist_sd_prec_log"))
    if (!nrow(P)) { message("Sin pares; se omite."); next }
    Pb <- a_bloques(P, base, 500)    # inferencia primaria: bloques de 500 km

    for (nm in names(escenarios)) {
      ctrl <<- escenarios[[nm]]$ctrl
      r <- pruebas_pares(Pb, tam_km, hipotesis = escenarios[[nm]]$hip)
      res_todo[[paste(tam_km, etiqueta, nm)]]  <- r$pruebas %>%
        mutate(escenario = nm, variante = etiqueta, tipo_rar = cc$tipo, sorteo = cc$s)
      pend_todo[[paste(tam_km, etiqueta, nm)]] <- r$pendientes %>%
        mutate(escenario = nm, variante = etiqueta, tipo_rar = cc$tipo, sorteo = cc$s)
    }
  }
}

res_todo  <- bind_rows(res_todo)
pend_todo <- bind_rows(pend_todo)
write_csv(res_todo,  file.path(dir_mod, "rarefaccion_pruebas.csv"))
write_csv(pend_todo, file.path(dir_mod, "rarefaccion_generos.csv"))

# ---- LA COMPARACION QUE IMPORTA --------------------------------------------
# Pendiente por genero sin rarificar, frente a la media de los sorteos rarificados.
# Si las de Artibeus, Platyrrhinus y Sturnira se acercan a la de Carollia, es artefacto.
cat("\n\n=== Pendientes por genero: completo frente a rarificado (media de", n_sorteos, "sorteos) ===\n")
comp <- pend_todo %>%
  group_by(tam_km, escenario, respuesta, predictor, genero, grupo = tipo_rar) %>%
  summarise(b = mean(estimate), de = sd(estimate), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = grupo, values_from = c(b, de)) %>%
  mutate(cambio_techo = b_techo - b_completo, cambio_dist = b_dist - b_completo) %>%
  arrange(escenario, respuesta, predictor, tam_km, genero)
print(as.data.frame(comp %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

cat("\nLectura: 'cambio_techo' y 'cambio_dist' son cuanto se mueve la pendiente con cada rarefaccion.\n")
cat(referente, "no se rarifica nunca: es el control interno y no deberia moverse.\n")
cat("Si los otros tres se mueven hacia el valor de Carollia, la diferencia era de resolucion.\n")
cat("de_techo y de_dist son la dispersion entre sorteos; si son grandes, sube n_sorteos.\n")
cat("La variante 'dist' es la que tiene fuerza: iguala la mediana de S, no solo el maximo.\n")
