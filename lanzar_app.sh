#!/usr/bin/env bash
# Lanza SARA con el entorno de herramientas activo y en un puerto fijo, y abre
# la interfaz en una ventana sin barra de direcciones ni pestanas.
#
# El PATH del entorno conda es lo que permite que la comprobación de entorno de
# la app encuentre las ocho herramientas del pipeline.
#
# Para cambiar como se abre la interfaz:
#   SARA_UI=app        ventana de aplicación, sin adornos (por defecto)
#   SARA_UI=navegador  pestana normal en el navegador de siempre
#   SARA_UI=ninguno    no abre nada; se entra a mano por la URL
#   SARA_UI_SIZE=1280,800   tamaño inicial de la ventana
export PATH="/Users/usuario/miniforge3/envs/rnaseq_ecoli/bin:$PATH"
export SARA_HOST=127.0.0.1
export SARA_PORT=3838
cd "$(dirname "${BASH_SOURCE[0]}")"
exec Rscript app.R
