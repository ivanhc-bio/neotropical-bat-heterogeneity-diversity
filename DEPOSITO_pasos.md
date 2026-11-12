# Pasos para depositar en GitHub + Zenodo

Esto no va al repositorio: es la lista para ti.

## 1. Antes de subir nada

- [ ] Elige la licencia del **código** (MIT si quieres máxima reutilización,
      GPL-3.0 si quieres que las obras derivadas sigan siendo abiertas) y crea
      el archivo `LICENSE`. Ajusta el campo `license` de `CITATION.cff`.
- [ ] **Revisa el asunto de la licencia de los datos.** Tu descarga de GBIF es
      CC BY-NC 4.0. Eso restringe el uso comercial de todo lo derivado de ella,
      y es más restrictivo de lo que las revistas suelen asumir. Dos salidas:
      (a) depositar las tablas derivadas también como CC BY-NC 4.0, que es lo
      que dejé puesto, o (b) rehacer la descarga filtrando a conjuntos de datos
      con licencia CC BY o CC0, lo que reduciría los registros. La (a) es la
      razonable; solo hay que declararla, y ya está declarada en el manuscrito.
- [ ] Completa las siete filas `[VERIFICAR]` del diccionario del README.
- [ ] Completa la tabla de scripts del README con los que faltan (01–04, 06–09,
      12–16, 19).
- [ ] Corre `sessionInfo()` y guarda la salida en `sessionInfo.txt`.
- [ ] Añade tu ORCID y los de tus coautores a `CITATION.cff` y `.zenodo.json`.
- [ ] Revisa que ningún script lleve rutas absolutas tuyas (`D:/capitulo4_diversity`).
      Lo más limpio es un `config.R` con una sola línea `base <- "."` que los
      demás lean con `source()`.

## 2. Qué subir y qué no

El `.gitignore` que incluí deja fuera los rásteres y las descargas crudas de
GBIF, que son grandes y se reconstruyen desde los DOI. Sube:

- todos los scripts de `R/`
- `datos/rejilla/`, `datos/modelos/`, `datos/predictores/` (las tablas, no las capas)
- `datos/curacion/` si cabe; si pesa mucho, deja solo el último paso
- `figuras/` y `tablas/`

GitHub rechaza archivos de más de 100 MB. Si alguna tabla se pasa, súbela solo
a Zenodo, que admite hasta 50 GB por depósito, y dilo en el README.

## 3. Crear el repositorio

```bash
cd /ruta/a/capitulo4_diversity
git init
git add .
git commit -m "Data and code for the Neotropical frugivorous bats paper"
git branch -M main
git remote add origin https://github.com/ivanhc-bio/NOMBRE-DEL-REPO.git
git push -u origin main
```

Nombre sugerido: `heterogeneidad-murcielagos-neotropicales` o
`neotropical-bat-heterogeneity`. Que sea descriptivo: va a vivir en la cita.

## 4. Enlazar Zenodo

1. Entra a https://zenodo.org con tu cuenta de GitHub.
2. Ve a **Settings → GitHub** y activa el interruptor del repositorio.
   Importante: Zenodo solo archiva *releases* creadas **después** de activarlo.
3. En GitHub, **Releases → Create a new release**, etiqueta `v1.0.0`,
   título "Version accompanying the Journal of Biogeography submission".
4. Zenodo genera el DOI en unos minutos. Vas a recibir **dos**: uno de la
   versión (`v1.0.0`) y uno de *concept* que siempre apunta a la más reciente.
   **Usa el de concept en el manuscrito**, para que la cita no caduque si
   publicas una versión corregida tras la revisión.
5. Edita el depósito en Zenodo si hace falta: `.zenodo.json` rellena los
   metadatos automáticamente, pero conviene revisarlos a ojo.

## 5. Cerrar el manuscrito

Sustituye en la sección *Data and code availability*:

- `10.5281/zenodo.XXXXXXX` por el DOI de concept
- `https://github.com/ivanhc-bio/REPO` por la URL real

Y en el README, el badge y la cita del artículo.

## 6. Después de la revisión

Si cambias análisis durante la revisión, haz una `v1.1.0` y una release nueva.
El DOI de concept seguirá apuntando a la última, así que el manuscrito no
necesita cambios.
