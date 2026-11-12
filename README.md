# Environmental heterogeneity and the diversity of Neotropical frugivorous bats

Data and code for:

> Hernández-Chávez, I., & León-Paniagua, L. (2026). Environmental
> heterogeneity and the diversity of Neotropical frugivorous bats: a shared
> response in richness, genus-specific responses in turnover.
> *Journal of Biogeography*. [DOI del artículo]

Archived release: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## What this is

Taxonomic and phylogenetic alpha and beta diversity of four frugivorous
phyllostomid genera (*Artibeus*, including *Dermanura*; *Carollia*;
*Platyrrhinus*; *Sturnira*), with the nectarivorous *Anoura* as a secondary
contrast, modelled against present topographic heterogeneity and historical
climatic instability in equal-area cells of 50, 100 and 200 km across the
Neotropics.

The repository holds the curation pipeline, the modelling scripts, the derived
tables the models read, and the code that builds every figure and table in the
paper. It does **not** hold the raw climate and elevation rasters, which are
large and publicly archived; see **Data sources** for how to obtain them.

## Repository layout

```
datos/
  curacion/     record curation, one .rds per step (see Supplementary Table S1)
  rejilla/      grids, cell assignments, presence tables
  predictores/  per-cell environmental predictors
  modelos/      model input tables and model outputs
  filogenia/    pruned trees and null communities
  ambiente/     DEM and climate layers (NOT distributed; see below)
scripts/        analysis scripts, run in numerical order
config.R        project paths (inside scripts/)
figuras/        figures and tables as produced by script 30
tablas/         summary tables
```

## Reproducing the analysis

Run the scripts in `scripts/` in numerical order. Each writes its outputs to `datos/`
and reads only what earlier scripts wrote.

| Script | What it does |
|---|---|
| `config.R` | project paths; every script loads it with `source("config.R")` |
| `00_recorte_trace_remoto.R` / `.py` | downloads only the Neotropical window of CHELSA-TraCE21k (bio01, bio12) from the remote GeoTIFFs |
| `01_descarga_gbif.R` | GBIF download of the Phyllostomidae in a single request |
| `02_taxonomia_referencia.R` | reference species list from the Mammal Diversity Database v2.5 |
| `03_curacion_coordenadas.R` | coordinate curation of the specimens matched to MDD |
| `03b_auditoria_paises.R` | country audit against the MDD country lists |
| `04_rejilla_celdas.R` | nested equal-area grids; assignment of specimens to cells |
| `04b_filtro_rango.R` | distance of each record to its species' range polygon |
| `05_predictores_celdas.R` | per-cell environmental predictors at the three cell sizes |
| `05b_historia_climatica.R` | historical instability, anomaly and velocity from the TraCE21k series |
| `05c_esfuerzo_gbif.R`, `05c_esfuerzo_gbif_api.R` | sampling-effort layer, independent of the focal genera (file and API versions) |
| `05c_topo_rango.R` | elevation range per cell, the alternative topographic measure |
| `06_celdas_incluidas.R` | cell-inclusion criteria |
| `07_injertar_arboles.R` | prunes the 500 trees and grafts the missing species |
| `08_ses_filogenetico.R` | MPD, PD and their standardized effect sizes |
| `09_bases_modelos.R` | builds `base_modelos_{size}km.csv` |
| `10_modelos.R` | per-genus GAMs |
| `11_modelo_conjunto.R` | joint GAM with the genus interaction |
| `12_diversidad_beta.R` | Simpson turnover, taxonomic and phylogenetic |
| `13_modelos_beta.R` | GAMs of beta diversity |
| `15_robustez_espacial.R` | robustness to the treatment of space |
| `16_sensibilidad_datos.R` | robustness to data decisions |
| `17_gdm_beta.R` | generalized dissimilarity models (exploratory) |
| `18_h2_especificacion.R`, `18b_h2_wald_controles_por_genero.R` | diagnosis of the H2 interaction |
| `19_figura_tamanos_efecto.R` | effect-size figure |
| `20_alfa_pruebas.R` | **confirmatory alpha tests** (8 per cell size) |
| `20_topo_rugosidad_fina.R` | fine roughness measures (Supplementary Table S7) |
| `21_beta_pruebas.R` | **confirmatory beta tests** (8 per cell size) |
| `22_holm_16.R` | Holm correction over the 16 tests of each cell size |
| `23_anoura.R` | *Anoura* contrasts (exploratory) |
| `24_beta_robustez.R` | beta robustness checks |
| `25_alfa_robustez.R` | alpha robustness checks |
| `26_rango_estricto.R`, `26a_estructura.R` | strict-range variant (richness only) |
| `27_diagnosticos.R` | concurvity, Moran's I, VIF (Supplementary Table S6) |
| `28_precipitacion.R` | exploratory precipitation analyses (alpha) |
| `29_choco.R` | regional check |
| `30_figuras_tablas.R` | every figure and table in the paper |
| `30a_satura.R`, `30b_sin_saturadas.R` | checks on the truncated precipitation pixels |
| `31_beta_precipitacion.R` | beta models with differentiation in precipitation variability |
| `31_tabla_s1.R`, `31b_resumen.R` | the counts of Supplementary Table S1 |
| `32_rarefaccion_carollia.R` | rarefaction check of the *Carollia* contrast |
| `33_metodos_numeros.R` | the counts reported in Methods |
| `34_carollia.R` | the two *Carollia* diagnostics of Section 4 |

Numbers 14, 14b and 14c are absent: they were method-development runs that are
not part of any result reported here.

Scripts are run from the `scripts/` directory. `config.R` resolves every path
from `CAP4_BASE` or, failing that, from the location of the `datos/` folder, so
no script contains an absolute path. The large rasters are not in the
repository; point `CAP4_TRACE`, `CAP4_CHELSA` and `CAP4_TOPO` at them, or place
them under `datos/ambiente/` as `config.R` expects.

### Software

Analyses were run in R with **mgcv** (GAMs), **dplyr**, **readr**, **sf**,
**ape**, **vegan** (null communities), **gdm**, **ggplot2** and **patchwork**.
Run `sessionInfo()` and paste the output into `sessionInfo.txt` before the
release, so the archived version records exact package versions.

## Data sources

| Source | What | Reference |
|---|---|---|
| GBIF | 343,553 Phyllostomidae occurrence records from 680 datasets, downloaded 30 October 2026 via **rgbif** | https://doi.org/10.15468/dl.gfvm4t |
| GBIF | 3,921,476 Mammalia records used as the sampling-effort layer, downloaded 6 November 2026 | https://doi.org/10.15468/dl.2u7xpj |
| EarthEnv-DEM (GMTED2010) | elevation, ~1 km | Amatulli et al. (2018), https://doi.org/10.1038/sdata.2018.40 |
| CHELSA v2.1 | present climate, 1981–2010 | Brun et al. (2022), https://doi.org/10.5194/essd-14-5573-2022 |
| CHELSA-TraCE21k | transient climate, 45 layers at 500-year steps over 22,000 years | Karger et al. (2023), https://doi.org/10.5194/cp-19-439-2023 |
| Mammal Diversity Database v2.5 | reference taxonomy | https://doi.org/10.5281/zenodo.17033774 |
| VertLife / Upham et al. | 500 mammal trees | https://doi.org/10.1371/journal.pbio.3000494 |
| Marsh et al. (2022) | expert range polygons | https://doi.org/10.1111/jbi.14330 |

## Variable dictionary

### `datos/rejilla/presencia_{50,100,200}km.csv`

One row per genus, cell and species.

| Column | Type | Description |
|---|---|---|
| `celda` | chr | cell identifier, `c<size>_<col>_<row>` in the LAEA grid |
| `genero_analisis` | chr | genus as analysed (*Dermanura* folded into *Artibeus*) |
| `especie` | chr | species name in the MDD v2.5 reference taxonomy |
| `n_reg` | int | preserved-specimen records of that species in that cell |

A `_rango_estricto` variant of each file holds the strict-range sensitivity
analysis, recomputed without records beyond the 100 km polygon margin.

### `datos/modelos/base_modelos_{50,100,200}km.csv`

One row per genus and cell. At 100 km: 1,402 rows, 490 cells, 5 genera.

**Identifiers and geography**

| Column | Type | Range (100 km) | Description |
|---|---|---|---|
| `celda` | chr | 490 values | cell identifier |
| `genero` | chr | 5 values | genus |
| `tam_km` | int | 100 | cell size |
| `criterio` | chr | `principal` | cell-inclusion criterion |
| `x_centro_km`, `y_centro_km` | int | −3550…3950, −3050…3550 | cell centre, LAEA km (lat_0 = 0, lon_0 = −75) |
| `lon_centro`, `lat_centro` | dbl | −110.8…−38.2, −27.6…31.0 | cell centre, degrees |

**Responses**

| Column | Type | Range | Description |
|---|---|---|---|
| `S` | int | 1…10 | species of the genus recorded in the cell |
| `log_S` | dbl | 0…2.30 | natural log of `S` |
| `S_ses` | int | 1…10 | species of the genus in the cell that are present in the phylogeny |
| `mpd_obs` | dbl | 2.32…17.5 | observed mean pairwise phylogenetic distance (NA when S < 2) |
| `pd_obs` | dbl | 5.15…42.4 | observed phylogenetic diversity |
| `ses_mpd`, `ses_pd` | dbl | −4.67…2.55 | standardized effect sizes, mean across the 500 trees |
| `ses_mpd_de_arboles`, `ses_pd_de_arboles` | dbl | 0.04…1.72 | **[VERIFICAR]** between-tree standard deviation of the SES |
| `rango_mpd`, `rango_pd` | dbl | 0.0002…0.999 | **[VERIFICAR]** rank of the observed value within the null distribution |
| `ses_mpd_sin_hermanas`, `ses_pd_sin_hermanas` | dbl | −4.66…2.32 | SES excluding species grafted next to a sister |
| `ses_mpd_solo_adn`, `ses_pd_solo_adn` | dbl | −4.90…2.43 | SES using only species with DNA sequences |
| `filo_valida` | lgl | T/F | whether the cell enters the phylogenetic models (S_ses ≥ 3) |

**Hypothesis predictors**

| Column | Type | Range | Description |
|---|---|---|---|
| `topo_sd_elev` | dbl | 3.3…1587 | **topographic heterogeneity**: SD of elevation within the cell, m |
| `hist_sd_temp` | dbl | 0.58…3.08 | **historical thermal instability**: SD of mean annual temperature across the 45 TraCE21k layers, °C |

**Controls and alternatives**

| Column | Type | Range | Description |
|---|---|---|---|
| `ctrl_esfuerzo_log` | dbl | 1.23…4.89 | sampling effort: log10 of Mammalia records in the cell, 1960–2025, plus one |
| `ctrl_area_tierra_km2` | dbl | 5048…10000 | land area in the cell |
| `clima_bio01_media`, `clima_bio12_media` | dbl | 5.8…28.1, 362…6518 | mean annual temperature (°C) and precipitation (mm) |
| `topo_elev_media` | dbl | 2.1…4071 | mean elevation, a control in the beta models |
| `topo_tri_media` | dbl | 0.46…229 | mean terrain ruggedness index |
| `clima_bio01_sd` | dbl | 0.07…7.07 | within-cell SD of mean annual temperature (r = 0.985 with `topo_sd_elev`) |
| `clima_bio12_sd` | dbl | 31…1852 | within-cell SD of annual precipitation |
| `hist_anom_temp` | dbl | 0.54…8.25 | temperature anomaly since the Last Glacial Maximum |
| `hist_sd_prec_mm` | dbl | 18…2786 | SD of annual precipitation across time steps, mm |
| `hist_sd_prec_log` | dbl | 0.02…0.57 | SD across time steps of log(precipitation + 1), negatives clamped to zero |
| `hist_anom_prec` | dbl | 0.006…0.88 | **[VERIFICAR]** precipitation anomaly since the LGM |
| `hist_vel_media_log10` | dbl | −0.38…2.01 | **[VERIFICAR]** mean climate velocity, log10 |
| `clima_bio01_sd_res` | dbl | −1.53…1.44 | **[VERIFICAR]** residual of `clima_bio01_sd` on `topo_sd_elev` |
| `hab_evi_media` | dbl | 0.12…0.59 | **[VERIFICAR]** mean EVI in the cell |
| `n_reg_genero` | int | 5…2023 | records of the genus in the cell |
| `esfuerzo_genero_log` | dbl | 1.61…7.61 | natural log of `n_reg_genero` |
| `esfuerzo_focal_log` | dbl | 1.04…3.43 | **[VERIFICAR]** log10 of records of the focal genera in the cell |
| `modelo_riqueza`, `modelo_filo` | chr | `completo`/`reducido` | which model variant applies to that row |

> **[POR COMPLETAR]** Las siete filas marcadas `[VERIFICAR]` las deduje del
> nombre y del rango, no de tu código. Confírmalas o corrígelas: una definición
> equivocada en el diccionario es peor que ninguna.

## Licence

Code: [MIT / GPL-3.0 — elige una] (`LICENSE`).

Data: the derived occurrence tables descend from a GBIF download released under
**CC BY-NC 4.0**, so they are redistributed here under the same terms. Note that
this is more restrictive than the CC BY 4.0 that journals usually assume for
supplementary data.

## Contact

Iván Hernández-Chávez — https://github.com/ivanhc-bio
Posgrado en Ciencias Biológicas, UNAM, and Laboratorio de Mastozoología
Evolutiva y Colecciones Científicas, Facultad de Ciencias, UNAM.

Corresponding author of the article: Livia León-Paniagua (llp@ciencias.unam.mx).

## Funding

Secretaría de Ciencia, Humanidades, Tecnología e Innovación (SECIHTI),
grant CF-2023-G-222.
