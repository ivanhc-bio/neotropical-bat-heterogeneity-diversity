# 12_diversidad_beta.R
# Diversidad beta por celda y género: recambio (turnover) taxonómico y filogenético de Simpson.
#   - Para cada celda: media del βsim contra las otras celdas del MISMO género que están dentro de un radio
#     (radios_km; 500 km por defecto, para evitar el efecto de borde de compararse con todo el continente) y contra todas
#     las celdas (radio Inf, comparable con el borrador). Se necesitan al menos min_vecinos vecinas con valor.
#   - βsim taxonómico = min(b, c) / (min(b, c) + a), con a = especies compartidas, b y c = exclusivas de cada celda.
#   - βsim filogenético = igual, con a, b y c como longitud de rama compartida/exclusiva (PD, con raíz). Se calcula con
#     una matriz de incidencia especies x ramas: a = P diag(len) P', donde P marca las ramas presentes en cada celda.
#     Esto equivale a phylo.beta.pair pero permite cientos de árboles y cientos de nulos.
#   - Nulo filogenético: se permutan las etiquetas de las puntas (qué especie ocupa qué punta), manteniendo la composición
#     de cada celda y los vecinos. SES = (obs - media nulo) / DE nulo, por árbol; se promedia entre árboles.
#   - Con pocas especies en el género (Carollia, Platyrrhinus) hay pocas permutaciones distintas: el nulo es grueso.
#   El βsim taxonómico depende de la riqueza (más especies, más recambio posible): usa log_S como covariable al modelarlo.
#
# Salidas (datos/modelos): beta_{etiqueta}_{tam}km.csv, una fila por celda x género x radio.
#
# No probado: no tengo R en mi entorno. Prueba con n_arboles_beta <- 5 y n_nulos <- 49, y luego sube.
# Requiere: ape, dplyr, readr, tidyr

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(ape)
library(dplyr)
library(readr)
library(tidyr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_rej  <- file.path(dir_base, "datos", "rejilla")
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_arb  <- file.path(dir_base, "datos", "filogenia")

etiqueta        <- "completo"
archivo_arboles <- file.path(dir_arb, sprintf("arboles_%s_500.rds", etiqueta))
todas_las_escalas <- TRUE    # TRUE = corre 100, 50 y 200 km (en ese orden); FALSE = solo 100 km
tamanos_km <- if (todas_las_escalas) c(100, 50, 200) else c(100)
generos         <- c("Anoura", "Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
radios_km       <- c(500, Inf)       # Inf = todas las celdas del género
min_vecinos     <- 3                 # vecinas con valor para calcular la media de una celda
n_arboles_beta  <- 100               # árboles (de los 500) para el βsim filogenético
n_nulos         <- 499               # permutaciones de etiquetas por árbol
semilla         <- 1

# ---- 1. Funciones -------------------------------------------------------------------
incidencia <- function(tr) {          # matriz especies x ramas (1 si la especie desciende de la rama)
  n <- Ntip(tr)
  pp <- prop.part(tr)
  hijo <- tr$edge[, 2]
  A <- matrix(0, n, length(hijo))
  for (e in seq_along(hijo)) {
    h <- hijo[e]
    A[if (h <= n) h else pp[[h - n]], e] <- 1
  }
  list(A = A, len = tr$edge.length)
}

# βsim por pares a partir de la matriz de compartido (a) y el total por celda (S)
beta_simpson <- function(a, S) {
  b <- pmax(outer(S, rep(1, length(S))) - a, 0)
  c <- t(b)
  m <- pmin(b, c); den <- m + a
  out <- m / den; out[den <= 0] <- NA
  out
}

# media por fila contra las vecinas (Wb = matriz 0/1 de vecindad); NA si hay menos de min_v vecinas con valor
valor_celda <- function(M, Wb, min_v) {
  valido <- !is.na(M); M[!valido] <- 0
  nv <- rowSums(Wb * valido)
  out <- rowSums(Wb * M) / nv
  out[nv < min_v] <- NA
  out
}

beta_taxonomico <- function(X) { a <- X %*% t(X); beta_simpson(a, rowSums(X)) }

beta_filogenetico <- function(X, A, len) {
  P <- ((X %*% A) > 0) * 1
  a <- (P * rep(len, each = nrow(P))) %*% t(P)
  beta_simpson(a, as.vector(P %*% len))
}

# ---- 2. Datos comunes -----------------------------------------------------------------
arboles <- readRDS(archivo_arboles)
arboles <- arboles[seq_len(min(n_arboles_beta, length(arboles)))]
tips_arbol <- arboles[[1]]$tip.label
cat("Árboles usados:", length(arboles), "| nulos por árbol:", n_nulos, "\n")
set.seed(semilla)

for (tam_km in tamanos_km) {
  cat("\n=== Celdas de", tam_km, "km ===\n")
  base <- read_csv(file.path(dir_mod, sprintf("base_modelos_%dkm.csv", tam_km)), show_col_types = FALSE)
  pres <- read_csv(file.path(dir_rej, sprintf("presencia_%dkm.csv", tam_km)), show_col_types = FALSE)
  salida <- list()

  for (g in generos) {
    d <- base %>% filter(genero == g)
    p <- pres %>% filter(genero_analisis == g, celda %in% d$celda, especie %in% tips_arbol)
    X <- p %>% distinct(celda, especie) %>% mutate(v = 1) %>%
      pivot_wider(names_from = especie, values_from = v, values_fill = 0)
    cel <- X$celda; X <- as.matrix(X[, -1]); rownames(X) <- cel
    esp <- colnames(X)
    if (length(esp) < 3 || nrow(X) < 10) { cat("  [", g, "] pocas especies o celdas; se omite\n"); next }
    d <- d %>% filter(celda %in% cel) %>% arrange(match(celda, cel))
    stopifnot(identical(d$celda, cel))
    nc <- nrow(X)
    cat(sprintf("  [%s] celdas: %d | especies: %d\n", g, nc, length(esp)))

    dm <- as.matrix(dist(d[, c("x_centro_km", "y_centro_km")]))
    Wb <- lapply(radios_km, function(r) { w <- (dm <= r) * 1; diag(w) <- 0; w })
    names(Wb) <- as.character(radios_km)

    # taxonómico (una vez)
    Mt <- beta_taxonomico(X)
    tax <- lapply(Wb, function(w) valor_celda(Mt, w, min_vecinos))

    # filogenético: árboles x nulos
    obs_arb <- lapply(Wb, function(w) matrix(NA_real_, length(arboles), nc))
    ses_arb <- lapply(Wb, function(w) matrix(NA_real_, length(arboles), nc))
    nul_arb <- lapply(Wb, function(w) matrix(NA_real_, length(arboles), nc))
    for (t in seq_along(arboles)) {
      tr <- keep.tip(arboles[[t]], esp)
      inc <- incidencia(tr); Xt <- X[, tr$tip.label, drop = FALSE]
      obs <- lapply(Wb, function(w) valor_celda(beta_filogenetico(Xt, inc$A, inc$len), w, min_vecinos))
      s1 <- lapply(Wb, function(w) numeric(nc)); s2 <- s1; nn <- s1
      for (k in seq_len(n_nulos)) {
        Mk <- beta_filogenetico(Xt[, sample.int(ncol(Xt)), drop = FALSE], inc$A, inc$len)
        for (r in names(Wb)) {
          v <- valor_celda(Mk, Wb[[r]], min_vecinos); ok <- !is.na(v)
          s1[[r]][ok] <- s1[[r]][ok] + v[ok]; s2[[r]][ok] <- s2[[r]][ok] + v[ok]^2; nn[[r]][ok] <- nn[[r]][ok] + 1
        }
      }
      for (r in names(Wb)) {
        mu <- s1[[r]] / nn[[r]]
        de <- sqrt(pmax((s2[[r]] / nn[[r]] - mu^2) * nn[[r]] / (nn[[r]] - 1), 0))
        obs_arb[[r]][t, ] <- obs[[r]]; nul_arb[[r]][t, ] <- mu
        ses_arb[[r]][t, ] <- ifelse(de > 0, (obs[[r]] - mu) / de, NA)
      }
      if (t %% 10 == 0 || t == length(arboles)) cat(sprintf("    árbol %d de %d\n", t, length(arboles)))
    }

    for (r in names(Wb)) {
      salida[[paste(g, r)]] <- tibble(
        celda = cel, genero = g, tam_km = tam_km, radio_km = as.numeric(r), S = rowSums(X),
        n_vecinos = rowSums(Wb[[r]]), beta_tax = tax[[r]],
        beta_filo_obs = colMeans(obs_arb[[r]], na.rm = TRUE),
        beta_filo_nulo = colMeans(nul_arb[[r]], na.rm = TRUE),
        beta_filo_ses = colMeans(ses_arb[[r]], na.rm = TRUE),
        beta_filo_ses_de_arboles = apply(ses_arb[[r]], 2, sd, na.rm = TRUE)) %>%
        mutate(across(starts_with("beta_filo"), ~ ifelse(is.nan(.x), NA_real_, .x)))
    }
  }

  res <- bind_rows(salida)
  write_csv(res, file.path(dir_mod, sprintf("beta_%s_%dkm.csv", etiqueta, tam_km)))
  cat("\nGuardado:", file.path(dir_mod, sprintf("beta_%s_%dkm.csv", etiqueta, tam_km)), "\n")

  cat("\nResumen (mediana por género y radio; SES > 0 = más recambio filogenético que el azar):\n")
  print(as.data.frame(res %>% group_by(genero, radio_km) %>%
    summarise(celdas_con_valor = sum(!is.na(beta_tax)),
              r_tax_logS = round(suppressWarnings(cor(beta_tax, log(pmax(S, 1)), use = "complete.obs")), 2),
              mediana_beta_tax = round(median(beta_tax, na.rm = TRUE), 2),
              mediana_beta_filo = round(median(beta_filo_obs, na.rm = TRUE), 2),
              mediana_ses_filo = round(median(beta_filo_ses, na.rm = TRUE), 2),
              prop_ses_mayor_0 = round(mean(beta_filo_ses > 0, na.rm = TRUE), 2), .groups = "drop")),
    row.names = FALSE)
}
