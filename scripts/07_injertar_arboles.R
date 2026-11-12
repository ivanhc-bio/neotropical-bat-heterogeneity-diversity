# 07_injertar_arboles.R
# Prepara el conjunto de árboles filogenéticos para los análisis de diversidad filogenética:
#   1) lee una muestra de árboles de Upham et al. (2019) descargada de vertlife.org (p. ej. 500 árboles);
#   2) conserva solo las especies focales (lista MDD v2.5), renombrando las puntas con el nombre de MDD v2.5
#      según la tabla equivalencias_especies.csv;
#   3) INJERTA, en cada árbol por separado, las especies focales que no están en el árbol, dentro del clado de su
#      género, de modo que cada árbol del conjunto lleva una colocación distinta (así la incertidumbre de la
#      colocación se propaga a los resultados);
#   4) guarda los árboles y una tabla de control con el origen de cada especie.
#
# Modos:
#   modo = "completo"  : árboles "completados" de Upham (con especies imputadas por Upham) + injertos propios.
#   modo = "solo_adn"  : NO se injerta nada; se poda a las especies con ADN (con_adn = TRUE en la tabla). Se puede usar
#                        con los MISMOS árboles completados (se podan las imputadas), sin descargar el conjunto
#                        "DNA-only" de vertlife. Las especies sin ADN quedan fuera de ese análisis de sensibilidad
#                        (hay que recalcular las celdas sin ellas).
#
# Colocación de las especies injertadas:
#   - Con 'hermana' en la tabla: una especie (rama terminal) o un clado de varias separadas por ';' (rama que lleva a su
#     ancestro común), en un punto al azar (uniforme) de esa rama.
#   - Sin 'hermana': se elige una rama al azar dentro del clado del género (probabilidad proporcional a su longitud) y
#     un punto al azar sobre ella. La rama que antecede al género no se usa, así el género no deja de ser monofilético.
#   En ambos casos el árbol sigue siendo ultramétrico (la nueva punta llega a la misma altura que las demás). Si no,
#   el script se detiene.
#
# No probado con árboles reales: no tengo R en mi entorno. Probar primero con n_arboles <- 5.
# Requiere: ape, phytools, dplyr, readr (install.packages(c("ape", "phytools", "dplyr", "readr")))

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(ape)
library(phytools)
library(dplyr)
library(readr)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_arb  <- file.path(dir_base, "datos", "filogenia")
dir.create(dir_arb, showWarnings = FALSE, recursive = TRUE)

archivo_arboles <- file.path(dir_arb, "arboles_vertlife_completados.nex")   # ajusta al nombre del archivo descargado
archivo_equiv   <- file.path(dir_arb, "equivalencias_especies.csv")
modo            <- "completo"        # "completo" o "solo_adn"
usar_hermanas   <- TRUE              # FALSE = todas las injertadas al azar dentro del género (sensibilidad)
n_arboles       <- 500               # prueba con 5 antes de la corrida final
semilla         <- 2026
etiqueta        <- paste0(modo, if (!usar_hermanas) "_sin_hermanas" else "")

prefijos_genero <- list(Anoura = "Anoura", Artibeus = c("Artibeus", "Dermanura"), Carollia = "Carollia",
                        Platyrrhinus = "Platyrrhinus", Sturnira = "Sturnira")

# ---- 1. Funciones -------------------------------------------------------------------
leer_arboles <- function(f) {
  x <- if (grepl("\\.(nex|nexus)$", f, ignore.case = TRUE)) read.nexus(f) else read.tree(f)
  if (inherits(x, "phylo")) x <- structure(list(x), class = "multiPhylo")
  x
}

# "Anoura_geoffroyi_PHYLLOSTOMIDAE_CHIROPTERA" -> "Anoura_geoffroyi"
base_punta <- function(x) sub("^([^_]+_[^_]+).*$", "\\1", x)

# altura de cada nodo sobre las puntas (válido en árboles ultramétricos)
altura_nodo <- function(tr) { d <- node.depth.edgelength(tr); max(d) - d }

# Las longitudes de rama de vertlife vienen redondeadas, así que las puntas pueden diferir en la última cifra y
# is.ultrametric() da FALSE. Si la desviación relativa es diminuta se corrige extendiendo las ramas terminales; si es
# grande, el script se detiene (algo estaría mal con el archivo).
hacer_ultrametrico <- function(tr, tol = 1e-4) {
  d <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]
  desv <- (max(d) - min(d)) / max(d)
  if (desv > tol) stop("El árbol no es ultramétrico (desviación relativa ", signif(desv, 3), ")")
  if (desv > 0) tr <- force.ultrametric(tr, method = "extend", message = FALSE)
  tr
}

# aristas del clado del género (sin la rama que lo antecede)
aristas_clado <- function(tr, especies_grupo) {
  presentes <- intersect(especies_grupo, tr$tip.label)
  if (length(presentes) == 0) stop("Ninguna especie del grupo está en el árbol: ", paste(especies_grupo, collapse = ", "))
  if (length(presentes) == 1) return(which(tr$edge[, 2] == match(presentes, tr$tip.label)))
  desc <- getDescendants(tr, getMRCA(tr, presentes))
  which(tr$edge[, 2] %in% desc)
}

# 'hermana' puede ser una especie ("Sturnira lilium") o un clado de varias separadas por ";" ("Sturnira luisi;Sturnira paulsoni").
# Con una sola especie presente en el árbol se usa su rama terminal; con varias, la rama que lleva a su ancestro común (MRCA),
# como en 04b_graft_missing_species.R del borrador. Las especies de la lista que no estén en el árbol se ignoran (se avisa una vez).
avisadas <- character()
injertar_una <- function(tr, nuevo, especies_grupo, hermana = NA_character_) {
  H <- altura_nodo(tr)
  sis <- if (is.na(hermana)) character() else trimws(strsplit(hermana, ";")[[1]])
  sis_presentes <- intersect(sis, tr$tip.label)
  if (length(sis) > length(sis_presentes) && !(nuevo %in% avisadas)) {
    message("Aviso: para ", nuevo, " faltan en el árbol las hermanas: ", paste(setdiff(sis, sis_presentes), collapse = ", "))
    avisadas <<- c(avisadas, nuevo)
  }
  if (length(sis_presentes) >= 1) {
    nodo <- if (length(sis_presentes) == 1) match(sis_presentes, tr$tip.label) else getMRCA(tr, sis_presentes)
    ar <- which(tr$edge[, 2] == nodo)
  } else {
    ars <- aristas_clado(tr, especies_grupo)
    ar <- ars[sample.int(length(ars), 1, prob = tr$edge.length[ars])]
    nodo <- tr$edge[ar, 2]
  }
  largo <- tr$edge.length[ar]
  pos <- runif(1, 0, largo)                       # distancia desde el nodo hacia la raíz (uniforme; el borrador usaba el punto medio)
  tr2 <- bind.tip(tr, nuevo, edge.length = H[nodo] + pos, where = nodo, position = pos)
  d <- node.depth.edgelength(tr2)
  if (abs(d[match(nuevo, tr2$tip.label)] - max(d)) > 1e-6 * max(d))
    stop("El injerto de ", nuevo, " rompió la ultrametría; revisa la versión de phytools (bind.tip).")
  tr2
}

# ---- 2. Tabla de equivalencias ------------------------------------------------------------
equiv <- read_csv(archivo_equiv, show_col_types = FALSE) %>%
  mutate(across(c(punta_arbol, hermana), ~ ifelse(is.na(.x) | .x == "", NA_character_, .x)))
stopifnot(all(c("especie_mdd", "genero_grupo", "punta_arbol", "estado", "hermana") %in% names(equiv)))
if (any(equiv$estado == "pendiente")) {
  print(as.data.frame(equiv %>% filter(estado == "pendiente") %>% select(especie_mdd, nota)), row.names = FALSE)
  stop("Hay especies con estado 'pendiente' en la tabla; resuélvelas (en_arbol o injertar) antes de seguir.")
}
if (!all(equiv$estado %in% c("en_arbol", "injertar"))) stop("Valores de 'estado' permitidos: en_arbol, injertar")

arboles <- leer_arboles(archivo_arboles)
cat("Árboles leídos:", length(arboles), "\n")
puntas1 <- base_punta(arboles[[1]]$tip.label)

# todas las puntas marcadas 'en_arbol' deben existir; si no, la tabla tiene un nombre mal escrito
falta <- setdiff(equiv$punta_arbol[equiv$estado == "en_arbol"], puntas1)
if (length(falta)) {
  if (modo == "completo") stop("Puntas de la tabla que no están en el árbol: ", paste(falta, collapse = ", "))
  cat("Modo solo_adn: especies marcadas en_arbol pero sin ADN (quedan fuera):", paste(falta, collapse = ", "), "\n")
}
# puntas de los géneros focales que la tabla no usa (posibles sinónimos o nombres distintos)
gen_arbol <- sub("_.*$", "", puntas1)
sin_uso <- setdiff(puntas1[gen_arbol %in% unlist(prefijos_genero)], equiv$punta_arbol)
if (length(sin_uso)) cat("Puntas de géneros focales en el árbol que la tabla NO usa:", paste(sin_uso, collapse = ", "), "\n")

# ---- 3. Injertar árbol por árbol ---------------------------------------------------------
set.seed(semilla)
idx <- if (length(arboles) > n_arboles) sample(length(arboles), n_arboles) else seq_along(arboles)
en_arbol <- equiv %>% filter(estado == "en_arbol", punta_arbol %in% puntas1)
# solo_adn: se conservan únicamente las especies con ADN (columna con_adn de la tabla); sirve tanto con árboles
# completados (se podan las imputadas por Upham) como con árboles solo ADN descargados de vertlife
if (modo == "solo_adn") en_arbol <- en_arbol %>% filter(con_adn %in% TRUE)
inj <- equiv %>% filter(estado == "injertar")
if (modo == "solo_adn") inj <- inj[0, ]
if (!usar_hermanas) inj$hermana <- NA_character_
inj <- inj[order(is.na(inj$hermana)), ]            # primero las que tienen hermana
grupo_de <- split(equiv$especie_mdd, equiv$genero_grupo)

res <- vector("list", length(idx))
for (k in seq_along(idx)) {
  tr <- arboles[[idx[k]]]
  tr$tip.label <- base_punta(tr$tip.label)
  tr <- keep.tip(tr, en_arbol$punta_arbol)
  tr$tip.label <- en_arbol$especie_mdd[match(tr$tip.label, en_arbol$punta_arbol)]
  tr <- hacer_ultrametrico(tr)
  for (i in seq_len(nrow(inj)))
    tr <- injertar_una(tr, inj$especie_mdd[i], grupo_de[[inj$genero_grupo[i]]], inj$hermana[i])
  res[[k]] <- tr
  if (k %% 50 == 0 || k == length(idx)) cat("  árboles listos:", k, "/", length(idx), "\n")
}
res <- structure(res, class = "multiPhylo")

# ---- 4. Controles y salidas ------------------------------------------------------------
tr1 <- res[[1]]
qc <- equiv %>%
  mutate(origen = case_when(
    !especie_mdd %in% tr1$tip.label ~ "ausente (solo_adn)",
    estado == "en_arbol" & con_adn %in% c(FALSE, "FALSE") ~ "en el árbol de Upham (imputada por Upham, sin ADN)",
    estado == "en_arbol" ~ "en el árbol de Upham (con ADN)",
    !is.na(hermana) ~ "injertada junto a su hermana",
    TRUE ~ "injertada al azar en el género")) %>%
  select(especie_mdd, genero_grupo, origen, hermana, con_adn)
cat("\nOrigen de las especies focales por género:\n")
print(as.data.frame(qc %>% count(genero_grupo, origen) %>% tidyr::pivot_wider(names_from = origen, values_from = n, values_fill = 0)), row.names = FALSE)
cat("\nUltramétricos:", all(vapply(res, is.ultrametric, logical(1), tol = 1e-5 * max(node.depth.edgelength(tr1)))), "\n")
cat("Especies por árbol:", unique(vapply(res, Ntip, numeric(1))), "\n")

write_csv(qc, file.path(dir_arb, paste0("qc_arboles_", etiqueta, ".csv")))
saveRDS(res, file.path(dir_arb, paste0("arboles_", etiqueta, "_", length(res), ".rds")))
write.tree(res, file.path(dir_arb, paste0("arboles_", etiqueta, "_", length(res), ".tre")))

# un árbol de ejemplo con las especies injertadas en rojo
col <- ifelse(tr1$tip.label %in% qc$especie_mdd[grepl("injertada", qc$origen)], "red", "black")
pdf(file.path(dir_arb, paste0("ejemplo_arbol_", etiqueta, ".pdf")), width = 8, height = 14)
plot(ladderize(tr1), tip.color = col[match(ladderize(tr1)$tip.label, tr1$tip.label)], cex = 0.6, no.margin = TRUE)
dev.off()
cat("Listo. Revisa el PDF de ejemplo: en rojo, las especies injertadas.\n")
