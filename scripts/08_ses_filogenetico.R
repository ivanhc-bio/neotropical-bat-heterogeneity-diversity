# 08_ses_filogenetico.R
# Diversidad filogenética alfa por género y celda: MPD y PD de Faith, con tamaño de efecto estandarizado (SES) frente
# a un modelo nulo, propagando la incertidumbre filogenética (cada árbol del conjunto de árboles se usa por separado).
#
# Qué hace:
#   - Para cada tamaño de celda y cada género, arma la matriz celdas x especies a partir de presencia_{tam}km.csv
#     (script 04), usando solo las celdas incluidas para ese género (script 06, columnas incl_<criterio>_<género>).
#   - Modelo nulo: se aleatoriza la matriz manteniendo la riqueza de cada celda y la frecuencia de cada especie
#     (vegan: "quasiswap"). El grupo de especies es el género, a escala continental (todas las celdas incluidas).
#     Opcional: grupos por región biogeográfica (archivo_regiones), aleatorizando dentro de cada región.
#   - Las matrices nulas no dependen del árbol: se generan UNA vez por género y se reutilizan con los 500 árboles.
#   - Para cada árbol se calcula MPD y PD (con raíz en el ancestro del género) observados y en los nulos, y el SES =
#     (observado - media nula) / DE nula. Se resume entre árboles: media del SES, DE entre árboles y rango (proporción
#     de nulos con valor menor que el observado, promediada entre árboles).
#   - SES solo en celdas con S >= s_min_ses (3); MPD observado desde S >= 2.
#
# Salidas (en datos/modelos): ses_{etiqueta}_{criterio}_{tam}km.csv, una fila por género y celda.
#
# Notas:
#   - La nula conserva las frecuencias observadas de cada especie, así que hereda en parte el sesgo de muestreo. Eso es
#     deliberado (compara la composición de cada celda con lo que cabría esperar dada la frecuencia observada de cada
#     especie), pero hay que decirlo en métodos.
#   - Con árboles de solo ADN (etiqueta "solo_adn") las especies sin ADN se quitan de la matriz antes de aleatorizar; las
#     celdas pueden quedar con menos especies.
#   - No probado: no tengo R en mi entorno. Prueba con tamanos_km <- 100, n_nulos <- 99, n_arboles_uso <- 5.
# Requiere: ape, vegan, dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(ape)
library(vegan)
library(dplyr)
library(readr)
library(tidyr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_arb  <- file.path(dir_base, "datos", "filogenia")

etiqueta         <- "completo"        # "completo", "completo_sin_hermanas" o "solo_adn_sin_hermanas" (nombre del .rds)
archivo_arboles  <- file.path(dir_arb, sprintf("arboles_%s_500.rds", etiqueta))
criterio         <- "principal"       # "principal", "estricto" o "laxo" (columnas incl_<criterio>_<género> del 06)
tamanos_km       <- c(200)            # primero 100; después 50 y 200
generos          <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
n_nulos          <- 499
s_min_ses        <- 3
n_arboles_uso    <- 500               # usa solo los primeros N árboles (para pruebas)
semilla          <- 2026
archivo_regiones <- NULL              # opcional: CSV con columnas celda, region (sensibilidad por región biogeográfica)

# ---- 1. Funciones -------------------------------------------------------------------
# Para cada fila de X (celdas x especies, 0/1): MPD y PD de Faith con raíz en el ancestro del grupo
# D: distancias cofenéticas (especies x especies); A: especies x ramas (1 si la especie desciende de la rama); len: ramas
metricas <- function(X, D, A, len) {
  S <- rowSums(X)
  mpd <- rowSums((X %*% D) * X) / (S * (S - 1))
  mpd[S < 2] <- NA
  pd <- as.numeric(((X %*% A) > 0) %*% len)
  list(mpd = as.numeric(mpd), pd = pd)
}

incidencia <- function(tr) {
  n <- Ntip(tr)
  pp <- prop.part(tr)                         # tips que descienden de cada nodo interno (n + 1, n + 2, ...)
  hijo <- tr$edge[, 2]
  A <- matrix(0, n, length(hijo))
  for (e in seq_along(hijo)) {
    h <- hijo[e]
    A[if (h <= n) h else pp[[h - n]], e] <- 1
  }
  list(A = A, len = tr$edge.length)
}

# media y desviación estándar por columna de una matriz (nulos x celdas)
media_col <- function(m) colMeans(m, na.rm = TRUE)
de_col <- function(m) { mu <- media_col(m); sqrt(pmax(colMeans((sweep(m, 2, mu))^2, na.rm = TRUE) * nrow(m) / (nrow(m) - 1), 0)) }

# ---- 2. Datos comunes -----------------------------------------------------------------
arboles <- readRDS(archivo_arboles)
if (is.finite(n_arboles_uso)) arboles <- arboles[seq_len(min(n_arboles_uso, length(arboles)))]
tips_arbol <- arboles[[1]]$tip.label
cat("Árboles usados:", length(arboles), "| especies en el árbol:", length(tips_arbol), "\n")

regiones <- if (!is.null(archivo_regiones)) read_csv(archivo_regiones, show_col_types = FALSE) else NULL

set.seed(semilla)
resultados <- list()

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  inc <- read_csv(file.path(dir_mod, sprintf("celdas_incluidas_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)

  for (g in generos) {
    col_inc <- sprintf("incl_%s_%s", criterio, g)
    if (!col_inc %in% names(inc)) stop("No existe la columna ", col_inc, " en celdas_incluidas_", tam_km, "km.csv")
    celdas_g <- inc$celda[inc[[col_inc]] %in% TRUE]

    p <- pres %>% filter(genero_analisis == g, celda %in% celdas_g)
    fuera <- setdiff(unique(p$especie), tips_arbol)
    if (length(fuera)) cat("  [", g, "] especies fuera del árbol (se omiten):", paste(fuera, collapse = ", "), "\n")
    p <- p %>% filter(especie %in% tips_arbol)

    # matriz celdas x especies (0/1); se quitan celdas sin especies tras el filtro
    X <- p %>% distinct(celda, especie) %>% mutate(v = 1L) %>%
      pivot_wider(names_from = especie, values_from = v, values_fill = 0L)
    cel <- X$celda; X <- as.matrix(X[, -1]); rownames(X) <- cel
    esp <- colnames(X)
    S <- rowSums(X)
    if (length(esp) < s_min_ses || sum(S >= s_min_ses) < 10) {
      cat("  [", g, "] pocas especies o celdas con S >=", s_min_ses, "; se omite\n"); next
    }
    idx <- which(S >= s_min_ses)
    cat(sprintf("  [%s] celdas: %d | especies: %d | celdas con S >= %d: %d\n", g, nrow(X), length(esp), s_min_ses, length(idx)))

    # ---- matrices nulas (una vez por género; no dependen del árbol)
    bloques <- if (is.null(regiones)) rep("todas", nrow(X)) else regiones$region[match(rownames(X), regiones$celda)]
    if (anyNA(bloques)) stop("Hay celdas sin región en archivo_regiones (", g, ")")
    sim <- array(0L, dim = c(nrow(X), ncol(X), n_nulos))
    for (b in unique(bloques)) {
      ib <- which(bloques == b)
      if (length(ib) < 2 || sum(colSums(X[ib, , drop = FALSE]) > 0) < 2) {
        cat("    bloque", b, "con pocas celdas o especies: se deja la matriz observada (sin aleatorizar)\n")
        for (k in seq_len(n_nulos)) sim[ib, , k] <- X[ib, ]
      } else {
        sim[ib, , ] <- simulate(nullmodel(X[ib, , drop = FALSE], "quasiswap"), nsim = n_nulos)
      }
    }
    Xn <- aperm(sim[idx, , , drop = FALSE], c(3, 1, 2))      # nulos x celdas x especies
    dim(Xn) <- c(n_nulos * length(idx), ncol(X)); storage.mode(Xn) <- "double"
    rm(sim)

    # ---- recorrer los árboles
    ni <- length(idx); nc <- nrow(X)
    acum <- list(mpd_obs = numeric(nc), pd_obs = numeric(nc),
                 s_mpd = numeric(ni), s_mpd2 = numeric(ni), r_mpd = numeric(ni),
                 s_pd = numeric(ni), s_pd2 = numeric(ni), r_pd = numeric(ni))
    t0 <- Sys.time()
    for (k in seq_along(arboles)) {
      tr <- keep.tip(arboles[[k]], esp)
      D <- cophenetic.phylo(tr)[esp, esp]
      ic <- incidencia(tr); A <- ic$A[match(esp, tr$tip.label), , drop = FALSE]
      obs <- metricas(X, D, A, ic$len)
      nul <- metricas(Xn, D, A, ic$len)
      for (m in c("mpd", "pd")) {
        nm <- matrix(nul[[m]], nrow = n_nulos)                  # nulos x celdas (solo las de S >= s_min_ses)
        mu <- media_col(nm); sdn <- de_col(nm)
        ses <- (obs[[m]][idx] - mu) / sdn; ses[!is.finite(ses)] <- NA
        rango <- colMeans(sweep(nm, 2, obs[[m]][idx], "<"), na.rm = TRUE)   # proporción de nulos MENORES al observado
        acum[[paste0("s_", m)]] <- acum[[paste0("s_", m)]] + ses
        acum[[paste0("s_", m, "2")]] <- acum[[paste0("s_", m, "2")]] + ses^2
        acum[[paste0("r_", m)]] <- acum[[paste0("r_", m)]] + rango
      }
      acum$mpd_obs <- acum$mpd_obs + ifelse(is.na(obs$mpd), 0, obs$mpd)
      acum$pd_obs <- acum$pd_obs + obs$pd
      if (k %% 50 == 0 || k == length(arboles)) cat("    árboles:", k, "/", length(arboles), "(",
                                                    round(difftime(Sys.time(), t0, units = "mins"), 1), "min )\n")
    }
    K <- length(arboles)
    ses_media <- function(s, s2) s / K
    ses_de <- function(s, s2) sqrt(pmax(s2 / K - (s / K)^2, 0) * K / max(K - 1, 1))
    out <- tibble(celda = rownames(X), genero = g, tam_km = tam_km, criterio = criterio, etiqueta_arboles = etiqueta,
                  S = S, mpd_obs = ifelse(S >= 2, acum$mpd_obs / K, NA_real_), pd_obs = acum$pd_obs / K,
                  ses_mpd = NA_real_, ses_mpd_de_arboles = NA_real_, rango_mpd = NA_real_,
                  ses_pd = NA_real_, ses_pd_de_arboles = NA_real_, rango_pd = NA_real_)
    out$ses_mpd[idx] <- ses_media(acum$s_mpd); out$ses_mpd_de_arboles[idx] <- ses_de(acum$s_mpd, acum$s_mpd2)
    out$rango_mpd[idx] <- acum$r_mpd / K
    out$ses_pd[idx] <- ses_media(acum$s_pd); out$ses_pd_de_arboles[idx] <- ses_de(acum$s_pd, acum$s_pd2)
    out$rango_pd[idx] <- acum$r_pd / K
    cat(sprintf("    SES-MPD medio: %.2f | SES-PD medio: %.2f | DE entre árboles (mediana): %.2f\n",
                mean(out$ses_mpd, na.rm = TRUE), mean(out$ses_pd, na.rm = TRUE), median(out$ses_mpd_de_arboles, na.rm = TRUE)))
    resultados[[paste(tam_km, g)]] <- out
  }
  res_tam <- bind_rows(resultados[startsWith(names(resultados), paste(tam_km, ""))])
  write_csv(res_tam, file.path(dir_mod, sprintf("ses_%s_%s_%dkm.csv", etiqueta, criterio, tam_km)))
  cat("Guardado:", file.path(dir_mod, sprintf("ses_%s_%s_%dkm.csv", etiqueta, criterio, tam_km)), "\n")
}
