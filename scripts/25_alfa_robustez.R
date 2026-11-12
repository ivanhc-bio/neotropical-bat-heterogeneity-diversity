# 25_alfa_robustez.R
# Robustez (nivel A) de las 8 pruebas de alfa del script 20: riqueza y SES-MPD frente a topografía (presente) e inestabilidad térmica (histórica).
# Se repiten las pruebas bajo variantes y se compara el signo y el intervalo de la pendiente común y el p de la interacción con la principal.
# Variantes (en el tamaño de celda principal; cambia tamanos_km para otros):
#   principal               script 20
#   espacio_ninguno         sin término espacial
#   espacio_rigido          suavizado espacial rígido (k = 7 por género)
#   espacio_bloque          efecto aleatorio de bloques de 500 x 500 km, sin suavizado
#   controles_por_genero    todos los controles con pendiente propia de cada género
#   esfuerzo_genero         esfuerzo = registros del género (esfuerzo_genero_log) en lugar de Mammalia
#   historia_anomalia       anomalía térmica desde el Último Máximo Glacial (hist_anom_temp) en lugar de hist_sd_temp
#   topo_rango              rango de elevación en lugar de la desviación estándar (si existe la columna; ver candidatos_rango)
#   solo_presente / solo_historica   un predictor de hipótesis por modelo
#   ses_sin_hermanas, ses_solo_adn, ses_pd   variantes de la filogenia (solo SES-MPD / SES-PD)
#   sin_<género>            exclusión de un género a la vez (la interacción pasa a 2 gl)
#   y las variantes de DATOS (criterios de inclusión estricto y laxo, taxonomía alternativa, filtro de rango estricto), si existen los
#   archivos base_modelos_{tam}km{sufijo}.csv (y presencia_{tam}km{sufijo}.csv): ajusta 'variantes_datos' al nombre real de tus archivos.
# Salidas (datos/modelos): alfa_robustez_{tam}km.csv (formato largo) y alfa_robustez_resumen.csv.
# Requiere el script 20 en la misma carpeta (se carga con source). No probado: no hay R en mi entorno.
# Prueba con las variantes espaciales primero (variantes_a_correr).

solo_funciones <- TRUE
source("20_alfa_pruebas.R")     # funciones, rutas y parámetros (hip, ctrl_riq, ctrl_ses, generos, n_min_modelo, ...)

tamanos_km <- c(50, 100, 200)
bloque_km  <- 500
candidatos_rango <- c("topo_rango_elev", "topo_range_elev", "topo_rango", "topo_rg_elev")
# Para correr solo algunas variantes (p. ej. solo_variantes <- "topo_rango" después de crear topo_rango_{tam}km.csv con el script 05c);
# NULL = todas. Con solo_variantes el resultado se guarda con el sufijo _parcial y no pisa a alfa_robustez_{tam}km.csv.
solo_variantes <- "rango_estricto"
# sufijo de los archivos de datos alternativos: AJUSTA estos nombres a los reales
variantes_datos <- c(criterio_estricto = "_estricto", criterio_laxo = "_laxo",
                     taxonomia_alternativa = "_taxonomia_alt", rango_estricto = "_rango_estricto")

# ---- funciones que sustituyen a las del script 20 (opciones nuevas) -------------------------------------------------------------
espacio_modo <- "principal"            # "principal", "ninguno", "rigido" o "bloque"
suaves <- function(d) {
  if (espacio_modo == "ninguno") return(character())
  if (espacio_modo == "bloque") return("s(bloque, bs = 're')")
  k <- k_por_genero(d); if (espacio_modo == "rigido") k[] <- 7
  sprintf("s(x_centro_km, y_centro_km, by = I_%s, bs = 'tp', k = %d)", names(k), as.integer(k))
}

alfa_respuesta <- function(d, y, fam, ctrl, etiqueta, tam_km, hipotesis = hip, ctrl_gen = FALSE) {
  h <- unname(hipotesis)
  if (!ctrl_gen) {
    m0 <- ajusta(fml(y, c(ctrl, h), character(), d), d, fam)
    m1 <- ajusta(fml(y, c(ctrl, h), character(), d, desv = h), d, fam)
    ms <- ajusta(fml(y, ctrl, h, d), d, fam)
  } else {                              # controles con pendiente propia de cada género
    m0 <- ajusta(fml(y, h, ctrl, d), d, fam)
    m1 <- ajusta(fml(y, h, ctrl, d, desv = h), d, fam)
    ms <- ajusta(fml(y, character(), c(ctrl, h), d), d, fam)
  }
  tab0 <- as.data.frame(summary(m0)$p.table)
  pruebas <- bind_rows(lapply(names(hipotesis), function(tipo) {
    p <- hipotesis[[tipo]]
    wi <- wald(m1, paste0("genero:", p))
    tibble(tam_km = tam_km, respuesta = etiqueta, tipo = tipo, predictor = p, n = nrow(d),
           generos = paste(levels(d$genero), collapse = ","),
           est_comun = tab0[p, 1], se_comun = tab0[p, 2], ic_inf = tab0[p, 1] - 1.96 * tab0[p, 2],
           ic_sup = tab0[p, 1] + 1.96 * tab0[p, 2], p_comun = tab0[p, 4], gl_int = wi[["gl"]], p_int = wi[["p"]])
  }))
  pruebas
}

alfa_tamano <- function(base, pres, tam_km, generos, hipotesis = hip, ctrl_r = ctrl_riq, ctrl_s = ctrl_ses,
                        y_ses = "ses_mpd", ctrl_gen = FALSE, solo = c("riqueza", "ses")) {
  h <- unname(hipotesis)
  base$bloque <- factor(paste(floor(base$x_centro_km / bloque_km), floor(base$y_centro_km / bloque_km)))
  base <- bind_rows(lapply(generos, function(g) {
    d <- base %>% filter(genero == g); if (!nrow(d)) return(NULL)
    pool <- pres %>% filter(genero_analisis == g, celda %in% d$celda) %>% summarise(n_distinct(especie)) %>% pull()
    d$S_fallo <- pool - d$S
    if (any(d$S_fallo < 0)) stop("Celdas con más especies que el grupo de especies: ", g)
    d
  }))
  out <- list()
  if ("riqueza" %in% solo) {
    vr <- c(h, ctrl_r)
    d0 <- base %>% filter(if_all(all_of(vr), ~ !is.na(.x)))
    ok <- names(which(table(d0$genero) >= n_min_modelo))
    dr <- preparar(d0 %>% filter(genero %in% ok), vr)
    out$riq <- alfa_respuesta(dr, "cbind(S, S_fallo)", binomial(), ctrl_r, "riqueza", tam_km, hipotesis, ctrl_gen)
  }
  if ("ses" %in% solo) {
    vs <- c(h, ctrl_s)
    d1 <- base %>% filter(!is.na(.data[[y_ses]]), if_all(all_of(vs), ~ !is.na(.x)))
    if (y_ses == "ses_mpd") d1 <- d1 %>% filter(filo_valida)
    ok <- names(which(table(d1$genero) >= n_min_modelo))
    ds <- preparar(d1 %>% filter(genero %in% ok), vs)
    out$ses <- alfa_respuesta(ds, y_ses, gaussian(), ctrl_s, y_ses, tam_km, hipotesis, ctrl_gen)
  }
  bind_rows(out)
}

# ---- variantes ------------------------------------------------------------------------------------------------------------------------
sust <- function(x, de, a) sub(paste0("^", de, "$"), a, x)
variantes <- function(base) {
  v <- list(
    principal            = list(),
    espacio_ninguno      = list(espacio = "ninguno"),
    espacio_rigido       = list(espacio = "rigido"),
    espacio_bloque       = list(espacio = "bloque"),
    controles_por_genero = list(ctrl_gen = TRUE),
    esfuerzo_genero      = list(ctrl_r = sust(ctrl_riq, "ctrl_esfuerzo_log", "esfuerzo_genero_log"),
                                ctrl_s = sust(ctrl_ses, "ctrl_esfuerzo_log", "esfuerzo_genero_log")),
    historia_anomalia    = list(hipotesis = c(presente = "topo_sd_elev", historica = "hist_anom_temp")),
    solo_presente        = list(hipotesis = c(presente = "topo_sd_elev")),
    solo_historica       = list(hipotesis = c(historica = "hist_sd_temp")),
    ses_sin_hermanas     = list(y_ses = "ses_mpd_sin_hermanas", solo = "ses"),
    ses_solo_adn         = list(y_ses = "ses_mpd_solo_adn", solo = "ses"),
    ses_pd               = list(y_ses = "ses_pd", solo = "ses"))
  rango <- intersect(candidatos_rango, names(base))
  if (length(rango)) v$topo_rango <- list(hipotesis = c(presente = rango[1], historica = "hist_sd_temp"))
  else cat("  (no hay columna de rango de elevación entre", paste(candidatos_rango, collapse = ", "), "; se omite topo_rango)\n")
  for (g in generos) v[[paste0("sin_", g)]] <- list(generos = setdiff(generos, g))
  v
}

correr <- function(nombre, opc, base, pres, tam_km) {
  espacio_modo <<- if (is.null(opc$espacio)) "principal" else opc$espacio
  on.exit(espacio_modo <<- "principal")
  cat(sprintf("\n-- %s --\n", nombre))
  r <- tryCatch(alfa_tamano(base, pres, tam_km,
                            generos   = if (is.null(opc$generos)) generos else opc$generos,
                            hipotesis = if (is.null(opc$hipotesis)) hip else opc$hipotesis,
                            ctrl_r    = if (is.null(opc$ctrl_r)) ctrl_riq else opc$ctrl_r,
                            ctrl_s    = if (is.null(opc$ctrl_s)) ctrl_ses else opc$ctrl_s,
                            y_ses     = if (is.null(opc$y_ses)) "ses_mpd" else opc$y_ses,
                            ctrl_gen  = isTRUE(opc$ctrl_gen),
                            solo      = if (is.null(opc$solo)) c("riqueza", "ses") else opc$solo),
                error = function(e) { cat("   falló:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  r %>% mutate(variante = nombre)
}

res_todo <- list()
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  fr <- file.path(dir_mod, sprintf("topo_rango_%dkm.csv", tam_km))      # rango de elevación por celda (script 05c)
  if (file.exists(fr)) base <- base %>% left_join(read_csv(fr, show_col_types = FALSE) %>% select(celda, topo_rango_elev), by = "celda")
  res <- list()
  vv <- variantes(base)                  # se arma una sola vez (avisa si falta la columna de rango)
  if (!is.null(solo_variantes)) vv <- vv[intersect(names(vv), solo_variantes)]
  for (nm in names(vv)) res[[nm]] <- correr(nm, vv[[nm]], base, pres, tam_km)
  # variantes de datos: archivos alternativos si existen
  for (nm in if (is.null(solo_variantes)) names(variantes_datos) else intersect(names(variantes_datos), solo_variantes)) {
    fb <- file.path(dir_mod, sprintf("base_modelos_%dkm%s.csv", tam_km, variantes_datos[[nm]]))
    fp <- file.path(dir_rej, sprintf("presencia_%dkm%s.csv", tam_km, variantes_datos[[nm]]))
    if (!file.exists(fb)) { cat(sprintf("\n-- %s: no existe %s; se omite --\n", nm, basename(fb))); next }
    bd <- read_csv(fb, show_col_types = FALSE)
    pd <- if (file.exists(fp)) read_csv(fp, show_col_types = FALSE) else { cat("   (sin archivo de presencia alterno; se usa el principal)\n"); pres }
    # si la base alterna no trae SES válido (p. ej. rango_estricto), solo riqueza
    hay_ses <- "ses_mpd" %in% names(bd) && any(!is.na(bd$ses_mpd)) &&
               (!"filo_valida" %in% names(bd) || any(bd$filo_valida))
    if (!hay_ses) cat("   (sin SES válido en esta base; solo se modela la riqueza)\n")
    res[[nm]] <- correr(nm, if (hay_ses) list() else list(solo = "riqueza"), bd, pd, tam_km)
  }
  res <- bind_rows(res)
  if (!nrow(res)) { cat("\n  (ninguna variante produjo resultados en", tam_km, "km)\n"); next }
  res <- res %>% select(variante, tam_km, respuesta, tipo, predictor, n, est_comun, ic_inf, ic_sup, p_comun, gl_int, p_int)
  write_csv(res, file.path(dir_mod, sprintf("alfa_robustez_%dkm%s.csv", tam_km, if (is.null(solo_variantes)) "" else "_parcial")))
  res_todo[[as.character(tam_km)]] <- res
}
todo <- bind_rows(res_todo)
write_csv(todo, file.path(dir_mod, if (is.null(solo_variantes)) "alfa_robustez_resumen.csv" else "alfa_robustez_resumen_parcial.csv"))
cat("\nResumen: pendiente común [IC] (p) | interacción p, por variante\n")
print(as.data.frame(todo %>% mutate(resumen = sprintf("%.3f [%.3f, %.3f] (p=%.3f) | int p=%.4f", est_comun, ic_inf, ic_sup, p_comun, p_int)) %>%
        select(tam_km, respuesta, tipo, predictor, variante, resumen) %>% arrange(tam_km, respuesta, tipo, variante)), row.names = FALSE)
