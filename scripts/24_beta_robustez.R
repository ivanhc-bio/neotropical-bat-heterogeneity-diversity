# 24_beta_robustez.R
# Robustez de las 8 pruebas de beta (script 21) frente a dos cosas:
#   1. Dependencia espacial: la covarianza robusta se agrupa por BLOQUES de 500 km en lugar de por celda. Cada par aporta a los
#      bloques de sus dos celdas. Con pocos bloques los errores crecen: es la comprobación conservadora de los p tan pequeños del 50 km.
#      (Sustituye a la permutación por bloques anunciada: da lo mismo que se quiere saber, el p con dependencia espacial, en un solo ajuste.)
#   2. Submuestreo de celdas: el script 21 usa hasta 400 celdas por género elegidas al azar; aquí se repite con otras dos semillas.
# Variantes: celda_s1 (= script 21), bloque_s1, celda_s2, bloque_s2, celda_s3, bloque_s3.
# Salidas (datos/modelos): beta_robustez_{tam}km.csv (formato largo) y beta_robustez_resumen.csv (todos los tamaños).
# Requiere correr antes el script 21 (comparte funciones). No probado: no hay R en mi entorno.
# Tarda: recalcula los pares por cada semilla (3 por tamaño). Prueba con tamanos_km <- 200, max_celdas <- 150, n_arboles_beta <- 3.

solo_funciones <- TRUE
source("21_beta_pruebas.R")     # funciones, rutas y parámetros (max_celdas, n_arboles_beta, ...)
select <- dplyr::select          # por si terra quedó cargado

tamanos_km <- c(100, 50, 200)
semillas   <- c(1, 2, 3)
bloque_km  <- 500

arboles <- readRDS(archivo_arboles); arboles <- arboles[seq_len(min(n_arboles_beta, length(arboles)))]
res_todo <- list()
for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  res <- list()
  for (s in semillas) {
    semilla <- s                                       # lo lee pares_tamano (set.seed(semilla))
    cat(sprintf("\n-- semilla %d --\n", s))
    P <- pares_tamano(base, pres, arboles, tam_km, generos)
    res[[paste0("celda_s", s)]]  <- pruebas_pares(P, tam_km)$pruebas %>% mutate(variante = paste0("celda_s", s))
    res[[paste0("bloque_s", s)]] <- pruebas_pares(a_bloques(P, base, bloque_km), tam_km)$pruebas %>% mutate(variante = paste0("bloque_s", s))
  }
  res <- bind_rows(res) %>% select(tam_km, variante, respuesta, tipo, predictor, n_pares, n_celdas, est_comun, ic_inf, ic_sup, p_comun, gl_int, p_int) %>%
    rename(n_grupos = n_celdas)                        # n_grupos = celdas (variante celda) o bloques (variante bloque)
  write_csv(res, file.path(dir_mod, sprintf("beta_robustez_%dkm.csv", tam_km)))
  res_todo[[as.character(tam_km)]] <- res
}
todo <- bind_rows(res_todo)
write_csv(todo, file.path(dir_mod, "beta_robustez_resumen.csv"))
cat("\nResumen: pendiente común (p) y p de la interacción por variante\n")
print(as.data.frame(todo %>% mutate(resumen = sprintf("%.3f (p=%.3f); int p=%.4f", est_comun, p_comun, p_int)) %>%
        select(tam_km, respuesta, tipo, variante, n_grupos, resumen) %>%
        arrange(tam_km, respuesta, tipo, variante)), row.names = FALSE)
