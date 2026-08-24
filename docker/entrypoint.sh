#!/usr/bin/env bash
# Punto de entrada del contenedor.
#
# Comprueba el entorno ANTES de arrancar la aplicación. El motivo no es
# ceremonial: si una herramienta del pipeline no está en el PATH, la aplicación
# arranca perfectamente y el fallo aparece veinte minutos después, a mitad de una
# ejecución, en forma de código de salida. Comprobarlo al arrancar convierte ese
# fallo tardio en un mensaje inmediato que nombra lo que falta.
set -Eeuo pipefail

APP_DIR="${APP_DIR:-/app}"
cd "$APP_DIR"

herramientas=(fastqc fastp multiqc bowtie2 samtools featureCounts salmon kallisto)

comprobar_entorno() {
  local faltan=()
  for t in "${herramientas[@]}"; do
    command -v "$t" >/dev/null 2>&1 || faltan+=("$t")
  done
  if [ ${#faltan[@]} -gt 0 ]; then
    printf 'AVISO: herramientas del pipeline ausentes en el PATH: %s\n' "${faltan[*]}" >&2
    printf 'El análisis a partir de una matriz de conteos seguira funcionando;\n' >&2
    printf 'la ejecución del pipeline desde FASTQ, no.\n' >&2
  else
    printf 'Entorno completo: %d herramientas del pipeline disponibles.\n' "${#herramientas[@]}"
  fi
}

case "${1:-app}" in
  app)
    comprobar_entorno
    printf 'SARA en http://%s:%s\n' \
      "${SARA_HOST:-0.0.0.0}" "${SARA_PORT:-3838}"
    exec Rscript app.R
    ;;
  verify)
    exec "$APP_DIR/docker/verify.sh"
    ;;
  versions)
    cat /opt/env/versions_r.txt /opt/env/versions_cli.txt
    ;;
  shell)
    exec /bin/bash
    ;;
  *)
    # Cualquier otra cosa se ejecuta tal cual: permite `docker run ... Rscript -e ...`
    exec "$@"
    ;;
esac
