# syntax=docker/dockerfile:1
#
# SARA (Shiny App for RNA-seq Analysis) — imagen de ejecución autocontenida.
#
# El problema que resuelve: la aplicación depende a la vez de R con una docena
# de paquetes de CRAN y Bioconductor, y de ocho herramientas de línea de comandos
# externas. `requirements.sh` instala todo eso, pero instala LAS VERSIONES
# VIGENTES EL DÍA DE LA INSTALACIÓN, no aquellas con las que la aplicación fue
# validada. La consecuencia práctica es que los informes pueden AUDITAR con que
# versiones se ejecuto un análisis, pero nadie puede RECONSTRUIR ese entorno.
# Esta imagen cierra esa brecha.
#
# Construcción en dos etapas, con un motivo concreto para cada una:
#   1. `tools`  — resuelve el entorno conda con las herramientas del pipeline.
#      Se separa porque el resolutor de dependencias es lo más lento y lo más
#      voluminoso del proceso, y porque así ni micromamba ni sus cachés acaban
#      dentro de la imagen final.
#   2. `runtime` — parte de la imagen oficial de Bioconductor, que ya trae R y
#      las dependencias de sistema de los paquetes de Bioconductor compiladas.
#      Compilar DESeq2 y sus dependencias desde cero en una imagen genérica de R
#      añade decenas de minutos por construcción.
#
# Uso:
#   docker build -t sara:1.0 .
#   docker run --rm -p 3838:3838 -v "$PWD/data:/data:ro" -v "$PWD/outputs:/app/outputs" sara:1.0

# ── Etapa 1: herramientas del pipeline ──────────────────────────────────────
FROM mambaorg/micromamba:1.5.8 AS tools

USER root
COPY docker/environment.yml /tmp/environment.yml
# `-p /opt/env` en lugar de un entorno con nombre: la ruta es la misma en las dos
# etapas, que es la condición para que los binarios copiados sigan encontrando
# sus bibliotecas.
RUN micromamba create -y -p /opt/env -f /tmp/environment.yml \
    && micromamba clean --all --yes \
    && find /opt/env -follow -type f -name '*.a' -delete \
    && find /opt/env -follow -type f -name '*.pyc' -delete

# ── Etapa 2: entorno de ejecución ───────────────────────────────────────────
# La release de Bioconductor fija también la versión de R y la compatibilidad
# cruzada entre paquetes. Es la versión con la que la aplicación se válido.
FROM bioconductor/bioconductor_docker:RELEASE_3_22 AS runtime

LABEL org.opencontainers.image.title="SARA" \
      org.opencontainers.image.description="Aplicación Shiny para el análisis integral de datos RNA-seq, desde FASTQ hasta enriquecimiento funcional." \
      org.opencontainers.image.licenses="MIT"

# Instantanea fechada de CRAN. Sin esto, shiny y bslib llegarian en la versión
# publicada el día de la construcción y dos imagenes con el mismo Dockerfile
# tendrían software distinto.
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-08-20
ARG BUILD_CPUS=4

COPY --from=tools /opt/env /opt/env
# El entorno conda va DELANTE en el PATH para que multiqc use su propio Python.
# No se toca LD_LIBRARY_PATH a propósito: hacerlo pondría las bibliotecas del
# entorno conda por delante de las que R ya tiene enlazadas y es una via
# conocida de romper paquetes compilados.
ENV PATH=/opt/env/bin:${PATH}

COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN CRAN_SNAPSHOT="${CRAN_SNAPSHOT}" BUILD_CPUS="${BUILD_CPUS}" \
    Rscript /tmp/install_r_packages.R && rm /tmp/install_r_packages.R

# Registro de las versiones de las herramientas CLI realmente instaladas. Junto
# con versions_r.txt es la procedencia del entorno: lo que la memoria del TFM
# puede citar como "esta imagen contiene exactamente esto".
RUN { for t in fastqc fastp multiqc bowtie2 samtools featureCounts salmon kallisto; do \
        printf '%s\t%s\n' "$t" "$(command -v "$t" >/dev/null 2>&1 && ("$t" --version 2>&1 | head -1) || echo AUSENTE)"; \
      done; } > /opt/env/versions_cli.txt \
    && cat /opt/env/versions_cli.txt

# El código de la aplicación va al FINAL, y en su propia capa: es lo único que
# cambia entre construcciones habituales, así que todo lo anterior se reutiliza
# desde la cache y una modificación en un modulo de R no reinstala el entorno.
WORKDIR /app
COPY --chown=1000:1000 . /app

# Usuario sin privilegios. La escritura se limita a outputs/, que es lo único
# que la aplicación necesita modificar.
RUN useradd --uid 1000 --create-home --shell /bin/bash app 2>/dev/null || true \
    && mkdir -p /app/outputs /data \
    && chown -R 1000:1000 /app/outputs \
    && chmod +x /app/workflow.sh /app/docker/entrypoint.sh /app/docker/verify.sh
USER 1000

# Escucha en todas las interfaces: con el 127.0.0.1 por defecto el servidor
# arranca pero es inalcanzable desde fuera del contenedor, que es el fallo más
# comun al contenerizar una aplicación Shiny.
ENV SARA_HOST=0.0.0.0 \
    SARA_PORT=3838
EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=5s --start-period=40s --retries=3 \
  CMD Rscript -e "con <- try(url(paste0('http://127.0.0.1:', Sys.getenv('SARA_PORT', '3838')), open = 'rb'), silent = TRUE); quit(status = if (inherits(con, 'try-error')) 1 else {close(con); 0})"

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["app"]
