# 30_figuras_tablas.R
# Genera las figuras y tablas principales (y las de apoyo del suplemento) a partir de los CSV de los scripts 20-25.
#   Tabla 1   -> tabla1_16_pruebas.csv            (16 pruebas por tamaño de celda: estimador, IC, p de Holm)   [texto principal]
#   Tabla 2   -> tabla2_robustez_alfa_100km.csv   (robustez de H1, H2 y H3 a 100 km)                           [texto principal]
#   Tabla S2  -> tablaS2_robustez_alfa_{50,200}km.csv                                                          [suplemento]
#   Figura 1  -> fig1_area_estudio.png/pdf        (predictores y número de géneros por celda, 100 km)         [texto principal]
#   Figura 2  -> fig2_riqueza_topografia.png/pdf  (riqueza contra topografía por género + pendientes por género)[texto principal]
#   Figura 3  -> fig3_pendientes_comunes.png/pdf  (pendiente común de las 8 combinaciones en los 3 tamaños)     [texto principal]
#   Figura 4  -> fig4_beta_por_genero.png/pdf     (pendientes por género de recambio)                         [texto principal]
#   Figura S1 -> figS1_alfa_por_genero.png/pdf    (pendientes por género de SES-MPD y de historia)            [suplemento]
# Entrada (datos/modelos): pruebas_16_holm.csv (script 22), alfa_robustez_resumen.csv (25), alfa_pendientes_genero_{tam}km.csv (20),
#   beta_pendientes_genero_{tam}km.csv (21), base_modelos_{tam}km.csv.
# Salida: dir_fig = <raiz del repositorio>/figuras (se crea). Textos de las figuras en inglés (van al manuscrito).
# Requiere: dplyr, readr, tidyr, ggplot2, patchwork. Opcional (costa en la Figura 1): sf, rnaturalearth. No probado: no hay R en mi entorno.
# Los IC de las pendientes por género de beta (Figura 4) vienen del script 21, que ahora agrupa por bloques de 500 km,
# igual que la inferencia principal.

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork)

# terra/raster enmascaran select(); forzamos el de dplyr
select <- dplyr::select

dir_base <- RUTA_BASE
dir_mod  <- file.path(dir_base, "datos", "modelos")
dir_fig  <- file.path(dir_base, "figuras"); dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)
tamanos_km <- c(50, 100, 200)
generos_frug <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")
pal_gen <- c(Artibeus = "#E69F00", Carollia = "#0072B2", Platyrrhinus = "#009E73", Sturnira = "#CC79A7", Anoura = "#999999")
pal_tam <- c("50" = "#56B4E9", "100" = "#0072B2", "200" = "#D55E00")
# nombres de género en cursivas: etiquetas de ejes (etq_gen) y de paneles (lab_gen)
etq_gen <- function(x) lapply(x, function(g) if (g %in% names(pal_gen)) bquote(italic(.(g))) else g)
lab_gen <- as_labeller(function(x) paste0("italic(", x, ")"), default = label_parsed)
tema <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
guarda <- function(p, nombre, w, h) {
  ggsave(file.path(dir_fig, paste0(nombre, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(dir_fig, paste0(nombre, ".pdf")), p, width = w, height = h)
}
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<.001", sub("^0", "", sprintf("%.3f", p))))
fmt_est <- function(e, lo, hi) ifelse(is.na(e), "", sprintf("%.3f [%.3f, %.3f]", e, lo, hi))

et_resp <- c(riqueza = "Richness", ses_mpd = "SES-MPD", beta_tax = "Taxonomic turnover", beta_filo = "Phylogenetic turnover")
et_tipo <- c(presente = "Present (topography)", historica = "Historical (thermal instability)")

# ---- Tabla 1 y Figura 3 ------------------------------------------------------------------------------------------------
largo <- read_csv(file.path(dir_mod, "pruebas_16_holm.csv"), show_col_types = FALSE) %>%
  mutate(respuesta = factor(respuesta, names(et_resp), et_resp), tipo = factor(tipo, names(et_tipo), et_tipo))
ancho <- read_csv(file.path(dir_mod, "pruebas_16_holm_ancho.csv"), show_col_types = FALSE) %>%
  mutate(respuesta = factor(respuesta, names(et_resp), et_resp), tipo = factor(tipo, names(et_tipo), et_tipo))
t1 <- ancho %>% arrange(respuesta, tipo, tam_km) %>%
  transmute(Response = respuesta, Predictor = tipo, `Cell size (km)` = tam_km,
            `Common slope [95% CI]` = fmt_est(estimado, ic_inf, ic_sup), `p` = fmt_p(p_comun), `Holm p` = fmt_p(holm_comun),
            `Genus interaction p` = fmt_p(p_interaccion), `Genus interaction Holm p` = fmt_p(holm_interaccion))
write_csv(t1, file.path(dir_fig, "tabla1_16_pruebas.csv"))

f3 <- largo %>% filter(prueba == "pendiente_comun") %>%
  mutate(sig = ifelse(p_holm < 0.05, "Significant after Holm (p < .05)", "Not significant after Holm"), tam = factor(tam_km, c(200, 100, 50)))
p3 <- ggplot(f3, aes(estimado, tam, colour = tipo, shape = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ic_inf, xmax = ic_sup), position = position_dodge(width = 0.6), size = 0.5) +
  scale_shape_manual(values = c("Significant after Holm (p < .05)" = 16, "Not significant after Holm" = 1), name = NULL) +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = "Predictor") +
  facet_wrap(~ respuesta, scales = "free_x", nrow = 1) +
  labs(x = "Slope shared by the four genera (per SD of the predictor; turnover on the logit scale)", y = "Cell size (km)") +
  tema + theme(legend.position = "bottom")
guarda(p3, "fig3_pendientes_comunes", 11, 4)

# ---- Tabla 2 y S2: robustez de alfa --------------------------------------------------------------------------------------
rob <- read_csv(file.path(dir_mod, "alfa_robustez_resumen.csv"), show_col_types = FALSE)
# Las variantes corridas por separado (script 25 con solo_variantes) quedan en el archivo _parcial.
ruta_parcial <- file.path(dir_mod, "alfa_robustez_resumen_parcial.csv")
if (file.exists(ruta_parcial)) {
  extra <- read_csv(ruta_parcial, show_col_types = FALSE)
  rob <- bind_rows(rob, extra) %>%
    distinct(tam_km, variante, respuesta, tipo, predictor, .keep_all = TRUE)
  cat("Variantes añadidas desde _parcial:", paste(unique(extra$variante), collapse = ", "), "\n")
}
et_var <- c(principal = "Main model", espacio_ninguno = "No spatial term", espacio_rigido = "Rigid spatial smooth (k = 7)",
            espacio_bloque = "Random effect of 500 km blocks", controles_por_genero = "Genus-specific slopes for controls",
            esfuerzo_genero = "Sampling effort of the genus", historia_anomalia = "Temperature anomaly instead of instability",
            solo_presente = "Present predictor only", solo_historica = "Historical predictor only",
            ses_sin_hermanas = "SES-MPD without grafted sister species", ses_solo_adn = "SES-MPD with DNA species only",
            ses_pd = "SES of Faith's PD", topo_rango = "Elevation range instead of SD of elevation",
            rango_estricto = "Records beyond the range polygon removed",
            sin_Artibeus = "Without Artibeus", sin_Carollia = "Without Carollia", sin_Platyrrhinus = "Without Platyrrhinus",
            sin_Sturnira = "Without Sturnira")
tabla_rob <- function(tam) {
  r <- rob %>% filter(tam_km == tam)
  h1 <- r %>% filter(respuesta == "riqueza", tipo == "presente") %>% select(variante, est1 = est_comun, lo1 = ic_inf, hi1 = ic_sup, p1 = p_comun, pint = p_int)
  h3 <- r %>% filter(grepl("^ses", respuesta), tipo == "historica") %>% select(variante, est3 = est_comun, lo3 = ic_inf, hi3 = ic_sup, p3 = p_comun)
  tibble(variante = names(et_var)) %>% left_join(h1, by = "variante") %>% left_join(h3, by = "variante") %>%
    filter(variante %in% r$variante) %>%
    transmute(Variant = unname(et_var[variante]),
              `H1: richness ~ topography [95% CI]` = fmt_est(est1, lo1, hi1), `H1 p` = fmt_p(p1),
              `H2: genus interaction p` = fmt_p(pint),
              `H3: SES-MPD ~ historical instability [95% CI]` = fmt_est(est3, lo3, hi3), `H3 p` = fmt_p(p3))
}
write_csv(tabla_rob(100), file.path(dir_fig, "tabla2_robustez_alfa_100km.csv"))
for (tam in c(50, 200)) write_csv(tabla_rob(tam), file.path(dir_fig, sprintf("tablaS2_robustez_alfa_%dkm.csv", tam)))
cat("Nota: los p de las tablas de robustez son SIN corregir; solo el modelo principal tiene p de Holm (Tabla 1).\n")

# ---- Figura 2: riqueza contra topografía por género ------------------------------------------------------------------------
b100 <- read_csv(file.path(dir_mod, "base_modelos_100km.csv"), show_col_types = FALSE) %>%
  filter(genero %in% generos_frug) %>% group_by(genero) %>%
  mutate(topo_z = (topo_sd_elev - mean(topo_sd_elev, na.rm = TRUE)) / sd(topo_sd_elev, na.rm = TRUE)) %>% ungroup()
p2a <- ggplot(b100, aes(topo_z, S, colour = genero)) +
  geom_point(alpha = 0.35, size = 1) + geom_smooth(method = "glm", method.args = list(family = quasipoisson), formula = y ~ x, se = TRUE, linewidth = 0.8) +
  scale_colour_manual(values = pal_gen, guide = "none") + facet_wrap(~ genero, scales = "free_y", nrow = 1, labeller = labeller(genero = lab_gen)) +
  labs(x = "Topographic heterogeneity (SD of elevation, standardized within genus)", y = "Species richness in the cell") + tema

pend_alfa <- bind_rows(lapply(tamanos_km, function(t) read_csv(file.path(dir_mod, sprintf("alfa_pendientes_genero_%dkm.csv", t)), show_col_types = FALSE)))
comun <- largo %>% filter(prueba == "pendiente_comun", faceta == "alfa") %>%
  transmute(tam_km, respuesta = as.character(respuesta), predictor, genero = "All genera (common)", estimate = estimado, ic_inf, ic_sup)
f2b <- bind_rows(pend_alfa %>% mutate(respuesta = unname(et_resp[respuesta])), comun) %>%
  filter(respuesta == "Richness", predictor == "topo_sd_elev") %>%
  mutate(genero = factor(genero, rev(c(generos_frug, "All genera (common)"))), tam = factor(tam_km, c(50, 100, 200)))
p2b <- ggplot(f2b, aes(estimate, genero, colour = tam)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ic_inf, xmax = ic_sup), position = position_dodge(width = 0.6), size = 0.45) +
  scale_colour_manual(values = pal_tam, name = "Cell size (km)") +
  scale_y_discrete(labels = etq_gen) +
  labs(x = "Slope of richness on topographic heterogeneity (per SD)", y = NULL) + tema + theme(legend.position = "bottom")
guarda(p2a / p2b + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "a"), "fig2_riqueza_topografia", 10, 8)

# ---- Figura 4 (y S4): pendientes por género ----------------------------------------------------------------------------------
forest_genero <- function(d, titulo_x) {
  ggplot(d %>% mutate(genero = factor(genero, rev(generos_frug)), tam = factor(tam_km, c(50, 100, 200))),
         aes(estimate, genero, colour = tam)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = ic_inf, xmax = ic_sup), position = position_dodge(width = 0.6), size = 0.4) +
    scale_colour_manual(values = pal_tam, name = "Cell size (km)") + scale_y_discrete(labels = etq_gen) + labs(x = titulo_x, y = NULL) + tema + theme(legend.position = "bottom")
}
pend_beta <- bind_rows(lapply(tamanos_km, function(t)
  read_csv(file.path(dir_mod, sprintf("beta_pendientes_genero_%dkm.csv", t)),
           show_col_types = FALSE))) %>%
  mutate(
    resp = factor(respuesta, c("beta_tax", "beta_filo"),
                  c("Taxonomic turnover", "Phylogenetic turnover")),
    pred = factor(
      ifelse(grepl("topo", predictor),
             "Topographic differentiation (absolute difference in the SD of elevation)",
             "Historical differentiation (absolute difference in thermal instability)"),
      levels = c("Topographic differentiation (absolute difference in the SD of elevation)",
                 "Historical differentiation (absolute difference in thermal instability)")))
p4 <- forest_genero(pend_beta, "Slope of turnover on environmental differentiation (logit scale, per SD)") + facet_grid(resp ~ pred)
guarda(p4, "fig4_beta_por_genero", 9, 6)

pend_a2 <- pend_alfa %>% filter(respuesta == "ses_mpd" | predictor == "hist_sd_temp") %>%
  mutate(resp = unname(et_resp[respuesta]), pred = ifelse(predictor == "topo_sd_elev", "topographic heterogeneity", "historical thermal instability"),
         panel = paste0(resp, " ~ ", pred))          # 3 paneles (richness ~ topografia ya esta en la Figura 2b)
pS4 <- forest_genero(pend_a2, "Slope (per SD of the predictor)") + facet_wrap(~ panel, scales = "free_x", ncol = 2)
guarda(pS4, "figS1_alfa_por_genero", 9, 6)

# ---- Figura 1: área de estudio ----------------------------------------------------------------------------------------------
bm <- read_csv(file.path(dir_mod, "base_modelos_100km.csv"), show_col_types = FALSE)
celdas <- bm %>% group_by(celda) %>% summarise(x = first(x_centro_km), y = first(y_centro_km), topo = first(topo_sd_elev), hist = first(hist_sd_temp),
                                              n_gen = n_distinct(genero[genero %in% generos_frug]), .groups = "drop")
costa <- tryCatch({
  library(sf); library(rnaturalearth)
  w <- ne_countries(scale = "medium", returnclass = "sf")
  st_transform(st_crop(w, xmin = -125, xmax = -30, ymin = -60, ymax = 35), "+proj=laea +lat_0=0 +lon_0=-75 +x_0=0 +y_0=0 +datum=WGS84 +units=km")
}, error = function(e) NULL)
mapa <- function(fill, nombre, escala) {
  p <- ggplot()
  if (!is.null(costa)) p <- p + geom_sf(data = costa, fill = NA, colour = "grey60", linewidth = 0.2)
  p + geom_tile(data = celdas, aes(x, y, fill = .data[[fill]]), width = 100, height = 100) + escala +
    labs(fill = nombre, x = NULL, y = NULL) + (if (!is.null(costa)) coord_sf(expand = FALSE) else coord_equal()) + tema +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position = "top", barwidth = unit(4.5, "cm"), barheight = unit(0.4, "cm")))
}
p1 <- mapa("topo", "SD of elevation (m)", scale_fill_viridis_c(option = "magma", direction = -1)) +
  mapa("hist", "Historical thermal instability (SD of temperature)", scale_fill_viridis_c(option = "viridis")) +
  mapa("n_gen", "Frugivorous genera analysed", scale_fill_viridis_c(option = "cividis", breaks = 1:4)) +
  plot_annotation(tag_levels = "a")
guarda(p1, "fig1_area_estudio", 12, 5)

cat("\nListo. Archivos en:", dir_fig, "\n")
