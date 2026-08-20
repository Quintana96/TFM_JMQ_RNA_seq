#!/usr/bin/env bash
# Verificacion de la imagen.
#
# Una imagen no es valida por construir sin error. Es valida si dentro de ella
# se cumple lo mismo que se cumplia en el entorno de desarrollo. Esto comprueba
# las tres cosas, en orden de coste creciente:
#
#   1. Que las herramientas del pipeline estan y responden.
#   2. Que los paquetes de R que la aplicacion da por presentes se cargan.
#   3. Que la bateria de tests pasa entera DENTRO del contenedor.
#
# El paso 3 es el que importa: es la misma bateria que valida el codigo en
# desarrollo, ejecutada sobre el entorno empaquetado.
set -Eeuo pipefail

APP_DIR="${APP_DIR:-/app}"
cd "$APP_DIR"

fallos=0
titulo() { printf '\n=== %s ===\n' "$1"; }

titulo "1. Herramientas del pipeline"
for t in fastqc fastp multiqc bowtie2 samtools featureCounts salmon kallisto; do
  if command -v "$t" >/dev/null 2>&1; then
    printf '  OK      %-14s %s\n' "$t" "$(command -v "$t")"
  else
    printf '  FALTA   %s\n' "$t"; fallos=$((fallos + 1))
  fi
done

titulo "2. Paquetes de R"
Rscript -e '
  req <- c("shiny","bslib","shinyjs","shinyFiles","DT","ggplot2","plotly","pheatmap",
           "processx","withr","DESeq2","edgeR","limma","tximport","clusterProfiler",
           "apeglm","fgsea")
  opc <- c("sva","IHW","qvalue","ashr","rtracklayer","dearseq","fishpond",
           "RNASeqPower","ReactomePA","reactome.db","org.Hs.eg.db","org.EcK12.eg.db")
  falta <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  sin_opc <- opc[!vapply(opc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(sin_opc)) cat("  opcionales ausentes:", paste(sin_opc, collapse=", "), "\n")
  if (length(falta)) { cat("  FALTAN requeridos:", paste(falta, collapse=", "), "\n"); quit(status = 1) }
  cat("  OK     ", length(req), "paquetes requeridos disponibles\n")
' || fallos=$((fallos + 1))

titulo "3. Bateria de tests"
Rscript tests/testthat.R || fallos=$((fallos + 1))

titulo "Resultado"
if [ "$fallos" -eq 0 ]; then
  printf 'Imagen verificada: los tres bloques han pasado.\n'
else
  printf 'Verificacion fallida: %d bloque(s) con errores.\n' "$fallos" >&2
fi
exit "$fallos"
