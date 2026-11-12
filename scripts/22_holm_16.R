# 22_holm_16.R
# Junta las pruebas de alfa (script 20) y de beta (script 21) y aplica Holm sobre las 16 pruebas de cada tamaño de celda:
#   2 predictores (presente, histórica) x 4 respuestas (riqueza, SES-MPD, recambio taxonómico, recambio filogenético) x 2 pruebas
#   (pendiente común, interacción género x predictor).
# Entrada (datos/modelos): alfa_pruebas_{tam}km.csv (script 20), beta_robustez_{tam}km.csv (script 24)
# Salida: pruebas_16_holm.csv (formato largo, una fila por prueba) y pruebas_16_holm_ancho.csv (una fila por predictor y respuesta)
# Requiere: dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr); library(readr); library(tidyr)

# terra/raster enmascaran select(); forzamos el de dplyr
select <- dplyr::select
dir_mod <- file.path(RUTA_BASE, "datos", "modelos")
tamanos_km <- c(50, 100, 200)

# Inferencia de beta: "bloque_s1" (principal: errores robustos agrupados por bloques de 500 km, script 24) o
# "celda_s1" (errores agrupados por celda, script 21; versión del suplemento). Ambas salen de beta_robustez_{tam}km.csv.
# OJO: la inferencia PRIMARIA del manuscrito es por bloques de 500 km.
# Con "celda_s1" la Tabla 1 sale con errores agrupados por celda, que son optimistas.
beta_inferencia <- "bloque_s1"
sufijo <- if (beta_inferencia == "bloque_s1") "" else paste0("_", beta_inferencia)

leer <- function(tam) {
  a <- read_csv(file.path(dir_mod, sprintf("alfa_pruebas_%dkm.csv", tam)), show_col_types = FALSE) %>%
    mutate(faceta = "alfa", n = as.numeric(n), inferencia = "GAM con suavizado espacial") %>% rename(n_unidades = n)
  b <- read_csv(file.path(dir_mod, sprintf("beta_robustez_%dkm.csv", tam)), show_col_types = FALSE) %>%
    filter(variante == beta_inferencia) %>%
    mutate(faceta = "beta", inferencia = beta_inferencia) %>% rename(n_unidades = n_pares) %>% select(-variante, -n_grupos)
  bind_rows(a, b)
}
todo <- bind_rows(lapply(tamanos_km, leer)) %>%
  mutate(respuesta = factor(respuesta, c("riqueza", "ses_mpd", "beta_tax", "beta_filo")))

largo <- bind_rows(
  todo %>% transmute(tam_km, faceta, inferencia, respuesta, tipo, predictor, prueba = "pendiente_comun", gl = 1, estimado = est_comun,
                     ic_inf, ic_sup, p = p_comun, n_unidades),
  todo %>% transmute(tam_km, faceta, inferencia, respuesta, tipo, predictor, prueba = "interaccion_genero", gl = gl_int, estimado = NA_real_,
                     ic_inf = NA_real_, ic_sup = NA_real_, p = p_int, n_unidades)) %>%
  group_by(tam_km) %>%
  mutate(n_pruebas_familia = n(), p_holm = p.adjust(p, method = "holm")) %>%
  ungroup() %>% arrange(tam_km, respuesta, tipo, prueba)
stopifnot(all(largo$n_pruebas_familia == 16))
write_csv(largo, file.path(dir_mod, sprintf("pruebas_16_holm%s.csv", sufijo)))

ancho <- largo %>% select(tam_km, respuesta, tipo, predictor, prueba, estimado, ic_inf, ic_sup, p, p_holm) %>%
  pivot_wider(names_from = prueba, values_from = c(estimado, ic_inf, ic_sup, p, p_holm)) %>%
  select(tam_km, respuesta, tipo, predictor, estimado = estimado_pendiente_comun, ic_inf = ic_inf_pendiente_comun,
         ic_sup = ic_sup_pendiente_comun, p_comun = p_pendiente_comun, holm_comun = p_holm_pendiente_comun,
         p_interaccion = p_interaccion_genero, holm_interaccion = p_holm_interaccion_genero)
write_csv(ancho, file.path(dir_mod, sprintf("pruebas_16_holm_ancho%s.csv", sufijo)))
print(as.data.frame(ancho %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
