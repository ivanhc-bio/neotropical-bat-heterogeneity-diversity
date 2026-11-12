# 19_figura_tamanos_efecto.R
# Figura única de tamaños de efecto con intervalos de confianza al 95% en los tres tamaños de celda (50, 100 y 200 km),
# a partir de robustez_espacial_{tam}km.csv (script 15, tratamiento espacial "actual").
#
#   Panel A: pendiente común de los predictores de hipótesis (riqueza: topo_sd_elev, clima_bio12_sd, hist_sd_temp;
#            SES-MPD: los mismos tres). Las pruebas confirmatorias (H1 y H3) van en negro y el resto, en gris (exploratorio).
#   Panel B: pendiente de topo_sd_elev por género sobre la riqueza (base de H2), con los cuatro frugívoros.
# Las pendientes son por DE dentro del género y del tamaño de celda: se comparan en signo e incertidumbre, no en magnitud
# entre tamaños (las celdas de 50 km están anidadas en las de 100 km y no son independientes).
#
# Entrada: datos/modelos/robustez_espacial_{tam}km.csv para los tamaños que existan (ejecuta antes 15 con todas_las_escalas = TRUE)
# Salida (datos/modelos): figura_tamanos_efecto.png (300 dpi), figura_tamanos_efecto.pdf y tabla_tamanos_efecto.csv
# No probado: no tengo R en mi entorno. Requiere: dplyr, readr, ggplot2

source("config.R")   # rutas del proyecto; ver README. Sustituye la ruta absoluta original.

library(dplyr)
library(readr)
library(ggplot2)

# ---- 0. Rutas y parámetros ------------------------------------------------------------
dir_base <- RUTA_BASE
dir_mod  <- file.path(dir_base, "datos", "modelos")

tamanos_km  <- c(50, 100, 200)
trat_figura <- "actual"
hip         <- c("topo_sd_elev", "clima_bio12_sd", "hist_sd_temp")
confirm     <- tibble(respuesta = c("riqueza", "ses_mpd"), termino = c("topo_sd_elev", "hist_sd_temp"),
                      hipotesis = c("H1", "H3"))
frugivoros  <- c("Artibeus", "Carollia", "Platyrrhinus", "Sturnira")

etiq_pred <- c(topo_sd_elev = "Elevation SD (topographic)", clima_bio12_sd = "Precipitation SD (climatic)",
               hist_sd_temp = "Temperature instability (historical)")
etiq_resp <- c(riqueza = "Species richness", ses_mpd = "SES-MPD")
col_tam   <- c("50" = "#E69F00", "100" = "#000000", "200" = "#0072B2")   # paleta apta para daltonismo

# ---- 1. Lectura -----------------------------------------------------------------------
leer <- function(tam) {
  f <- file.path(dir_mod, sprintf("robustez_espacial_%dkm.csv", tam))
  if (!file.exists(f)) { message("No encuentro ", f, "; se omite ", tam, " km"); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>% mutate(tam_km = tam)
}
d <- bind_rows(lapply(tamanos_km, leer))
if (!nrow(d)) stop("No hay archivos robustez_espacial_*km.csv en ", dir_mod)
d <- d %>% filter(tratamiento == trat_figura) %>%
  mutate(ic_inf = estimate - 1.96 * se, ic_sup = estimate + 1.96 * se, tam = factor(tam_km, levels = tamanos_km))

# ---- 2. Paneles ----------------------------------------------------------------------------
a <- d %>% filter(modelo == "comun", termino %in% hip) %>%
  left_join(confirm, by = c("respuesta", "termino")) %>%
  mutate(confirmatoria = !is.na(hipotesis),
         termino = factor(termino, levels = rev(hip), labels = rev(etiq_pred[hip])),
         respuesta = factor(respuesta, levels = names(etiq_resp), labels = etiq_resp))
b <- d %>% filter(modelo == "por_genero", respuesta == "riqueza", termino == "topo_sd_elev", genero %in% frugivoros) %>%
  mutate(genero = factor(genero, levels = rev(frugivoros)))

tema <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", axis.title.y = element_blank())
pos <- position_dodge(width = 0.6, reverse = TRUE)   # 50 km arriba y 200 km abajo, como en la leyenda

pa <- ggplot(a, aes(estimate, termino, colour = tam, shape = confirmatoria, group = tam)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = ic_inf, xmax = ic_sup), height = 0, position = pos, linewidth = 0.5) +
  geom_point(position = pos, size = 2.2) +
  facet_wrap(~ respuesta, scales = "free_x") +
  scale_colour_manual(values = col_tam, name = "Cell size (km)") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Confirmatory (H1, H3)", `FALSE` = "Exploratory"),
                     name = NULL) +
  labs(x = "Common slope per SD within genus (95% CI)", title = "A  Common effects of heterogeneity predictors") + tema

pb <- ggplot(b, aes(estimate, genero, colour = tam, group = tam)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = ic_inf, xmax = ic_sup), height = 0, position = pos, linewidth = 0.5) +
  geom_point(position = pos, size = 2.2) +
  scale_colour_manual(values = col_tam, name = "Cell size (km)") +
  labs(x = "Slope of elevation SD on richness per SD (95% CI)", title = "B  Genus-specific effect of topographic heterogeneity (H2)") +
  tema + theme(axis.text.y = element_text(face = "italic"))

# ---- 3. Guardar -----------------------------------------------------------------------------
tabla <- bind_rows(
  a %>% transmute(panel = "A", respuesta = as.character(respuesta), termino = as.character(termino), genero = "todos",
                  tam_km, estimate, se, ic_inf, ic_sup, p, confirmatoria),
  b %>% transmute(panel = "B", respuesta = "Species richness", termino = "Elevation SD (topographic)", genero = as.character(genero),
                  tam_km, estimate, se, ic_inf, ic_sup, p, confirmatoria = NA)) %>% arrange(panel, respuesta, termino, genero, tam_km)
write_csv(tabla, file.path(dir_mod, "tabla_tamanos_efecto.csv"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1.1, 1), guides = "collect") +
    patchwork::plot_annotation(caption = paste("Slopes are per SD within genus and cell size: compare sign and uncertainty, not magnitude across sizes or panels (axes differ).",
                                               "50-km cells are nested in 100-km cells. Intervals are unadjusted 95% CIs.", sep = "\n")) &
    theme(legend.position = "bottom")
} else {
  message("Instala patchwork (install.packages('patchwork')) para unir los paneles; se guardan por separado.")
  fig <- NULL
  ggsave(file.path(dir_mod, "figura_tamanos_efecto_A.png"), pa, width = 8, height = 4.5, dpi = 300)
  ggsave(file.path(dir_mod, "figura_tamanos_efecto_B.png"), pb, width = 8, height = 4, dpi = 300)
}
if (!is.null(fig)) {
  ggsave(file.path(dir_mod, "figura_tamanos_efecto.png"), fig, width = 8, height = 9, dpi = 300)
  ggsave(file.path(dir_mod, "figura_tamanos_efecto.pdf"), fig, width = 8, height = 9)
  print(fig)
}
cat("\nGuardado en", dir_mod, "\n")
print(as.data.frame(tabla %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
