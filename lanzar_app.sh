#!/usr/bin/env bash
# Lanza la aplicacion con el entorno de herramientas activo y en un puerto fijo.
# El PATH del entorno conda es lo que permite que la comprobacion de entorno de
# la app encuentre las ocho herramientas del pipeline.
export PATH="/Users/usuario/miniforge3/envs/rnaseq_ecoli/bin:$PATH"
export SHINY_HOST=127.0.0.1
export SHINY_PORT=3838
cd "$(dirname "${BASH_SOURCE[0]}")"
exec Rscript app.R
