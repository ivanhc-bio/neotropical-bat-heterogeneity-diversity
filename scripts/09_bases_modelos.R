# 09_bases_modelos.R
# Construye las bases de datos para los modelos (una por tamaño de celda, con todos los géneros en formato largo):
#   - respuestas: riqueza observada del género (S), SES-MPD y SES-PD (media entre árboles) y su incertidumbre, y las
#     versiones de las dos sensibilidades filogenéticas (sin hermanas, solo ADN);
#   - predictores finales (tabla de 05), coordenadas de la celda (para el suavizado espacial) y esfuerzo del género;
#   - marcas de qué modelo corresponde a cada género y respuesta (completo o reducido);
#   - un VIF por género y base, para confirmar que la colinealidad sigue siendo aceptable con menos celdas.
#
# Entradas (script 04, 05, 06 y 08):
#   datos/rejilla/celdas_{tam}km.csv                      coordenadas de las celdas
#   datos/predictores/predictores_{tam}km.csv             predictores
#   datos/modelos/celdas_incluidas_{tam}km.csv            inclusión por género, S_<género> y n_reg_<género>
#   datos/modelos/ses_{etiqueta}_{criterio}_{tam}km.csv   SES (principal y sensibilidades)
# Salidas (datos/modelos):
#   base_modelos_{tam}km.csv, diccionario_base_modelos.csv, vif_bases_{tam}km.csv
#
# No probado: no tengo R en mi entorno. Corre por tamaño y revisa lo que imprime.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(tidyr)
library(tibble)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_pred <- file.path(dir_base, "datos", "predictores")
dir_mod  <- file.path(dir_base, "datos", "modelos")

tamanos_km <- c(50, 100, 200)
criterio   <- "principal"
generos    <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")

# Modelo completo (10 predictores) y modelo reducido (5 predictores, uno por bloque de hipótesis + controles mínimos)
pred_completo <- c("topo_sd_elev", "clima_bio12_sd", "clima_bio01_sd_res", "hist_sd_temp", "hist_anom_prec",
                   "clima_bio01_media", "clima_bio12_media", "hab_evi_media", "ctrl_esfuerzo_log", "ctrl_area_tierra_km2")
pred_reducido <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp", "clima_bio01_media", "ctrl_esfuerzo_log")

# Modelo filogenético: completo solo donde hay muestra suficiente (celdas con S >= 3); reducido en los demás
filo_completo <- c("Artibeus", "Carollia")

etiquetas_ses <- c(principal = "completo", sin_hermanas = "completo_sin_hermanas", solo_adn = "solo_adn_sin_hermanas")

# Esfuerzo de muestreo independiente de los géneros focales (script 05c). Sustituye a ctrl_esfuerzo_log de 05, que se
# calcula con los registros de los propios géneros focales y es circular con la riqueza.
# Opciones: "esf_mam_log" (Mammalia), "esf_quir_log" (Chiroptera sin géneros focales),
#           "esf_quir_ej_log" (Chiroptera sin géneros focales, solo ejemplares de museo), o NULL para no sustituirlo.
esfuerzo_independiente <- "esf_mam_log"

vif_desde_datos <- function(d) {
  m <- cor(d, use = "complete.obs")
  tibble(variable = colnames(m), vif = round(diag(solve(m)), 2)) %>% arrange(desc(vif))
}

leer_ses <- function(etiqueta, tam_km) {
  f <- file.path(dir_mod, sprintf("ses_%s_%s_%dkm.csv", etiqueta, criterio, tam_km))
  if (!file.exists(f)) stop("Falta ", f)
  read_csv(f, show_col_types = FALSE)
}

# ---- 1. Una base por tamaño de celda ------------------------------------------------------------
vifs <- list()
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  cel  <- read_csv(file.path(dir_rej, sprintf("celdas_%dkm.csv", tam_km)), show_col_types = FALSE) %>%
    select(celda, x_centro_km, y_centro_km, lon_centro, lat_centro)
  pred <- read_csv(file.path(dir_pred, sprintf("predictores_%dkm.csv", tam_km)), show_col_types = FALSE)
  inc  <- read_csv(file.path(dir_mod, sprintf("celdas_incluidas_%dkm.csv", tam_km)), show_col_types = FALSE)
  falta <- setdiff(pred_completo, names(pred))
  if (length(falta)) stop("Faltan predictores en la tabla de 05: ", paste(falta, collapse = ", "))
  # esfuerzo independiente: se conserva el focal como esfuerzo_focal_log y ctrl_esfuerzo_log pasa a ser el independiente
  pred <- pred %>% mutate(esfuerzo_focal_log = ctrl_esfuerzo_log)
  if (!is.null(esfuerzo_independiente)) {
    f_esf <- file.path(dir_pred, sprintf("esfuerzo_independiente_%dkm.csv", tam_km))
    if (!file.exists(f_esf)) stop("Falta ", f_esf, ": corre 05c_esfuerzo_gbif.R o pon esfuerzo_independiente <- NULL")
    esf <- read_csv(f_esf, show_col_types = FALSE) %>% select(celda, esf_indep = all_of(esfuerzo_independiente))
    pred <- pred %>% left_join(esf, by = "celda")
    # todas las celdas incluidas en algún criterio fueron consultadas en 05c: un NA aquí es una consulta que falló o una
    # corrida de prueba (n_prueba finito), no un esfuerzo de cero. Una celda consultada sin registros tiene 0 de verdad.
    usadas <- inc$celda[rowSums(inc[, grep("^incl_", names(inc))], na.rm = TRUE) > 0]
    n_na <- sum(is.na(pred$esf_indep[pred$celda %in% usadas]))
    cat(sprintf("  Esfuerzo independiente: %d de %d celdas incluidas sin valor\n", n_na, length(usadas)))
    if (n_na > 0) stop("Faltan valores de esfuerzo independiente en ", n_na, " celdas incluidas de ", tam_km,
                       " km. Vuelve a correr 05c_esfuerzo_gbif_api.R con n_prueba <- Inf (reintenta las celdas con NA).")
    pred <- pred %>% mutate(ctrl_esfuerzo_log = esf_indep) %>% select(-esf_indep)
  }
  # Predictores extra para la corrida puente y las sensibilidades (no entran a los modelos principales); se conservan
  # solo los que existan en la tabla de 05
  pred_extra <- intersect(c("hist_sd_prec_mm", "hist_sd_prec_log", "clima_bio01_sd", "topo_elev_media", "topo_tri_media",
                            "hist_anom_temp", "hist_vel_media_log10"), names(pred))
  pred_completo_sel <- unique(c(pred_completo, "esfuerzo_focal_log", pred_extra))

  ses <- lapply(names(etiquetas_ses), function(k) leer_ses(etiquetas_ses[[k]], tam_km) %>%
                  select(celda, genero, mpd_obs, pd_obs, ses_mpd, ses_pd, ses_mpd_de_arboles, ses_pd_de_arboles,
                         rango_mpd, rango_pd, S_ses = S))
  names(ses) <- names(etiquetas_ses)
  alt <- function(k) ses[[k]] %>% select(celda, genero, ses_mpd, ses_pd) %>%
    rename_with(~ paste0(.x, "_", k), c(ses_mpd, ses_pd))

  bases <- list()
  for (g in generos) {
    col_inc <- sprintf("incl_%s_%s", criterio, g)
    b <- inc %>%
      filter(.data[[col_inc]] %in% TRUE) %>%
      transmute(celda, genero = g, tam_km = tam_km, criterio = criterio,
                S = .data[[paste0("S_", g)]], n_reg_genero = .data[[paste0("n_reg_", g)]],
                esfuerzo_genero_log = log(n_reg_genero), log_S = log(pmax(S, 1))) %>%
      left_join(cel, by = "celda") %>%
      left_join(pred %>% select(celda, all_of(pred_completo_sel)), by = "celda") %>%
      left_join(ses$principal, by = c("celda", "genero")) %>%
      left_join(alt("sin_hermanas"), by = c("celda", "genero")) %>%
      left_join(alt("solo_adn"), by = c("celda", "genero")) %>%
      mutate(filo_valida = !is.na(ses_mpd),
             modelo_riqueza = "completo",
             modelo_filo = ifelse(g %in% filo_completo, "completo", "reducido"))
    bases[[g]] <- b

    # VIF del conjunto en cada base (riqueza: todas las celdas; filogenético: celdas con SES)
    v1 <- vif_desde_datos(b[, pred_completo]) %>% mutate(base = "riqueza", modelo = "completo")
    bf <- b %>% filter(filo_valida)
    vf <- vif_desde_datos(bf[, if (g %in% filo_completo) pred_completo else pred_reducido]) %>%
      mutate(base = "filogenetica", modelo = ifelse(g %in% filo_completo, "completo", "reducido"))
    vifs[[paste(tam_km, g)]] <- bind_rows(v1, vf) %>% mutate(genero = g, tam_km = tam_km)
    cat(sprintf("  %-13s riqueza n = %3d | filogenética n = %3d (%s) | VIF máx: riqueza %.1f, filogenética %.1f\n",
                g, nrow(b), nrow(bf), ifelse(g %in% filo_completo, "completo", "reducido"),
                max(v1$vif), max(vf$vif)))
  }
  base_tam <- bind_rows(bases)
  write_csv(base_tam, file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)))
  write_csv(bind_rows(vifs[startsWith(names(vifs), paste(tam_km, ""))]),
            file.path(dir_mod, sprintf("vif_bases_%dkm.csv", tam_km)))
  cat("  Guardado: base_modelos_", tam_km, "km.csv (", nrow(base_tam), " filas)\n", sep = "")
}

# ---- 2. Diccionario de variables ------------------------------------------------------------
dic <- tribble(
  ~variable, ~descripcion,
  "celda, genero, tam_km, criterio", "Identificador de la celda, género, tamaño de celda (km) y criterio de inclusión (06)",
  "S", "Riqueza observada del género en la celda (especies con registro)",
  "n_reg_genero, esfuerzo_genero_log", "Registros del género en la celda y su logaritmo (esfuerzo específico del género)",
  "log_S", "Logaritmo de S; covariable en los modelos filogenéticos de Carollia y Sturnira (SES correlaciona con S)",
  "x_centro_km, y_centro_km, lon_centro, lat_centro", "Centro de la celda (LAEA, km) y en grados; para el suavizado espacial",
  "mpd_obs, pd_obs", "MPD y PD de Faith observados, medios entre árboles (MPD desde S >= 2)",
  "ses_mpd, ses_pd", "SES medio entre árboles (S >= 3): positivo = sobredispersión, negativo = agrupamiento",
  "ses_mpd_de_arboles, ses_pd_de_arboles", "Desviación estándar del SES entre árboles (incertidumbre filogenética)",
  "rango_mpd, rango_pd", "Proporción de nulos menores que el observado, promediada entre árboles",
  "ses_mpd_sin_hermanas, ses_pd_sin_hermanas", "SES con árboles completados sin hermanas (especies injertadas al azar en el género)",
  "ses_mpd_solo_adn, ses_pd_solo_adn", "SES solo con especies con ADN (sin imputadas ni injertadas); NA si S < 3 tras quitarlas",
  "filo_valida", "TRUE si la celda tiene SES (S >= 3) y entra en los modelos filogenéticos",
  "modelo_riqueza, modelo_filo", "Modelo que corresponde: completo (10 predictores) o reducido (5: topo_sd_elev, clima_bio12_sd, hist_sd_temp, clima_bio01_media, ctrl_esfuerzo_log)",
  "topo_sd_elev ... ctrl_area_tierra_km2", "Predictores finales de 05 (sin estandarizar)")
write_csv(dic, file.path(dir_mod, "diccionario_base_modelos.csv"))
cat("\nDiccionario guardado. Listo.\n")
