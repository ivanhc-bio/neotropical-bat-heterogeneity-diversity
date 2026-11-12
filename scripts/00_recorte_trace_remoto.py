#!/usr/bin/env python3
# 00_recorte_trace_remoto.py
# Versión en Python de 00_recorte_trace_remoto.R: recorta la ventana de interés de CHELSA-TraCE21k (bio01 y bio12)
# desde los GeoTIFF remotos y la agrega a ~5 km (media 5x5, ignorando los píxeles de relleno), pero:
#   - descarga varias capas A LA VEZ (hilos),
#   - lee por franjas, sin crear archivos temporales enormes,
#   - escribe la capa final con el MISMO nombre y carpeta que el script en R, así 05b no cambia.
# Es reanudable: salta lo que ya existe (escribe a .part y renombra al terminar).
#
# Requiere: pip install rasterio numpy
#
# Uso (desde una terminal):
#   python 00_recorte_trace_remoto.py --prueba                 # solo bio01 ids 20 y 0, mide tiempos
#   python 00_recorte_trace_remoto.py --verificar              # compara modo "bloques" y modo "rapido" en una ventana chica
#   python 00_recorte_trace_remoto.py --hilos 6                # todo, 6 capas en paralelo
#   python 00_recorte_trace_remoto.py --solo bio12 --hilos 6   # solo precipitación (lo que hay que rehacer)
#   python 00_recorte_trace_remoto.py --metodo rapido ...      # usa vistas reducidas (overviews) si el archivo las tiene
#
# Modos de lectura:
#   bloques (por defecto): lee la resolución completa (~1 km) por franjas y promedia 5x5 ignorando el relleno.
#       Equivale exactamente a lo que hacía el script en R (recortar, descartar relleno, agregar con na.rm).
#   rapido: pide a GDAL directamente la malla de ~5 km; si el archivo trae overviews lee muy poco. Solo úsalo después
#       de correr --verificar y ver que la diferencia con "bloques" es despreciable (el relleno podría contaminar
#       las vistas reducidas).

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window, from_bounds

# ---- Parámetros ------------------------------------------------------------------------------------------------
BASE_URL = "https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim"
DIR_SALIDA = os.environ.get("CAP4_TRACE", "../datos/ambiente/chelsa_trace_5km")
VENTANA = (-118.0, -38.0, -34.0, 33.0)          # xmin, ymin, xmax, ymax (grados)
FACTOR = 5                                       # 30 arcsec (~1 km) -> ~5 km
TIME_IDS = list(range(-200, 21, 5))              # 45 capas, una cada 500 años (1 unidad = 100 años; 20 = 1900-1990)
VARS = {
    "bio01": None,                               # sin relleno que descartar
    "bio12": 20000,                              # solo descarta rellenos tipo 65535; la lluvia real llega a ~9600 mm
}
FRANJA_FILAS = 400                               # filas de 1 km por franja; debe ser múltiplo de FACTOR
INTENTOS = 3


def fmt_id(i):
    return f"-{abs(i):03d}" if i < 0 else f"{i:04d}"


def ruta_origen(var, nom, base):
    ruta = f"{base}/{var}/{nom}"
    return "/vsicurl/" + ruta if base.startswith("http") else ruta


def entorno_gdal():
    return rasterio.Env(
        GDAL_HTTP_MAX_RETRY="6",
        GDAL_HTTP_RETRY_DELAY="5",
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_TIMEOUT="120",
        GDAL_HTTP_MULTIPLEX="YES",
        VSI_CACHE="TRUE",
        GDAL_CACHEMAX=256,
    )


def ventana_pixeles(src, ventana, factor):
    """Ventana entera, alineada a la malla del origen y con ancho y alto múltiplos de `factor`."""
    xmin, ymin, xmax, ymax = ventana
    w = from_bounds(xmin, ymin, xmax, ymax, src.transform)
    col0, row0 = int(round(w.col_off)), int(round(w.row_off))
    ncol = int(round(w.width)) // factor * factor
    nrow = int(round(w.height)) // factor * factor
    if ncol <= 0 or nrow <= 0:
        raise ValueError("la ventana pedida no cae dentro del raster")
    return Window(col0, row0, ncol, nrow)


def leer_valores(src, win, relleno, out_shape=None, resampling=Resampling.nearest):
    """Lee una ventana, aplica escala/desplazamiento y pone NaN en nodata y en valores > relleno."""
    a = src.read(1, window=win, masked=True, out_shape=out_shape, resampling=resampling).astype("float32")
    a = a.filled(np.nan) if np.ma.isMaskedArray(a) else a
    esc, off = src.scales[0], src.offsets[0]
    if esc != 1 or off != 0:
        a = a * np.float32(esc) + np.float32(off)
    if relleno is not None:
        a[a > relleno] = np.nan
    return a


def promedio_bloques(a, f):
    """Media f x f ignorando NaN (como aggregate(..., fun='mean', na.rm=TRUE))."""
    n, m = a.shape
    b = a.reshape(n // f, f, m // f, f)
    cuenta = np.sum(~np.isnan(b), axis=(1, 3))
    suma = np.nansum(b, axis=(1, 3), dtype="float64")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(cuenta > 0, suma / cuenta, np.nan)
    return out.astype("float32")


def agregar_capa(src, win, relleno, metodo):
    nrow, ncol = int(win.height) // FACTOR, int(win.width) // FACTOR
    if metodo == "rapido":
        # GDAL elige el overview adecuado; el relleno se descarta después (aproximación)
        return leer_valores(src, win, relleno, out_shape=(nrow, ncol), resampling=Resampling.average)
    out = np.full((nrow, ncol), np.nan, dtype="float32")
    for r0 in range(0, int(win.height), FRANJA_FILAS):
        h = min(FRANJA_FILAS, int(win.height) - r0)
        sub = Window(win.col_off, win.row_off + r0, win.width, h)
        out[r0 // FACTOR:(r0 + h) // FACTOR, :] = promedio_bloques(leer_valores(src, sub, relleno), FACTOR)
    return out


def procesar_capa(var, ident, relleno, base, dir_salida, metodo):
    nom = f"CHELSA_TraCE21k_{var}_{fmt_id(ident)}_V.1.0.tif"
    dir_v = os.path.join(dir_salida, var)
    os.makedirs(dir_v, exist_ok=True)
    salida = os.path.join(dir_v, nom)
    if os.path.exists(salida):
        return nom, "ya existe", 0.0
    t0 = time.time()
    with entorno_gdal(), rasterio.open(ruta_origen(var, nom, base)) as src:
        win = ventana_pixeles(src, VENTANA, FACTOR)
        datos = agregar_capa(src, win, relleno, metodo)
        tr0 = src.window_transform(win)
        perfil = dict(driver="GTiff", height=datos.shape[0], width=datos.shape[1], count=1, dtype="float32",
                      crs=src.crs, nodata=np.nan, compress="deflate", tiled=True, blockxsize=256, blockysize=256,
                      transform=rasterio.Affine(tr0.a * FACTOR, tr0.b, tr0.c, tr0.d, tr0.e * FACTOR, tr0.f))
    tmp = salida + ".part.tif"
    with rasterio.open(tmp, "w", **perfil) as dst:
        dst.write(datos, 1)
    os.replace(tmp, salida)
    return nom, "ok", time.time() - t0


def con_reintentos(*args):
    for k in range(1, INTENTOS + 1):
        try:
            return procesar_capa(*args)
        except Exception as e:  # noqa: BLE001
            print(f"  intento {k} falló para {args[0]} {args[1]}: {e}", file=sys.stderr)
            time.sleep(10)
    raise RuntimeError(f"no se pudo procesar {args[0]} {args[1]}")


def verificar(base, metodo_rapido=True):
    """Compara 'bloques' y 'rapido' en una ventana chica del Chocó con bio12 (id 20)."""
    global VENTANA
    VENTANA = (-79.0, 3.0, -76.0, 8.0)
    nom = f"CHELSA_TraCE21k_bio12_{fmt_id(20)}_V.1.0.tif"
    with entorno_gdal(), rasterio.open(ruta_origen("bio12", nom, base)) as src:
        print("overviews:", src.overviews(1), "| nodata:", src.nodata, "| escala/desplazamiento:", src.scales[0], src.offsets[0],
              "| dtype:", src.dtypes[0], "| tamaño:", src.width, "x", src.height)
        win = ventana_pixeles(src, VENTANA, FACTOR)
        t0 = time.time(); a = agregar_capa(src, win, VARS["bio12"], "bloques"); t_b = time.time() - t0
        t0 = time.time(); b = agregar_capa(src, win, VARS["bio12"], "rapido"); t_r = time.time() - t0
    ok = ~np.isnan(a) & ~np.isnan(b)
    print(f"bloques: {t_b:.1f} s | rapido: {t_r:.1f} s")
    print(f"NA bloques: {np.isnan(a).mean():.4f} | NA rapido: {np.isnan(b).mean():.4f}")
    print(f"max {np.nanmax(a):.0f} (bloques) vs {np.nanmax(b):.0f} (rapido)")
    print(f"diferencia absoluta media: {np.abs(a[ok] - b[ok]).mean():.2f} mm; máxima: {np.abs(a[ok] - b[ok]).max():.1f} mm")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prueba", action="store_true", help="solo bio01 ids 20 y 0")
    ap.add_argument("--verificar", action="store_true", help="compara los dos métodos en una ventana chica")
    ap.add_argument("--hilos", type=int, default=4)
    ap.add_argument("--solo", choices=list(VARS), help="procesa solo una variable")
    ap.add_argument("--metodo", choices=["bloques", "rapido"], default="bloques")
    ap.add_argument("--base", default=BASE_URL, help="URL base (o carpeta local para pruebas)")
    ap.add_argument("--dir-salida", default=DIR_SALIDA)
    args = ap.parse_args()

    if args.verificar:
        verificar(args.base)
        return
    variables = {args.solo: VARS[args.solo]} if args.solo else dict(VARS)
    ids = TIME_IDS
    if args.prueba:
        ids = [20, 0]
        variables = {"bio01": None}
    tareas = [(v, i, variables[v], args.base, args.dir_salida, args.metodo) for v in variables for i in ids]
    print(f"{len(tareas)} capas, {args.hilos} en paralelo, método: {args.metodo}")
    t0, fallos = time.time(), []
    with ThreadPoolExecutor(max_workers=args.hilos) as ex:
        futs = {ex.submit(con_reintentos, *t): t for t in tareas}
        for k, f in enumerate(as_completed(futs), 1):
            t = futs[f]
            try:
                nom, estado, seg = f.result()
                print(f"[{k}/{len(tareas)}] {nom}: {estado} ({seg/60:.1f} min)")
            except Exception as e:  # noqa: BLE001
                fallos.append(f"{t[0]} {t[1]}")
                print(f"[{k}/{len(tareas)}] FALLÓ {t[0]} {t[1]}: {e}", file=sys.stderr)
    print(f"\nTerminado en {(time.time() - t0)/60:.1f} min. Capas fallidas:", "; ".join(fallos) if fallos else "ninguna")


if __name__ == "__main__":
    main()
