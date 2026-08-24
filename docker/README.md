# Despliegue en contenedor

## Por qué

La aplicación depende a la vez de R con una docena de paquetes de CRAN y
Bioconductor, y de ocho herramientas de línea de comandos externas.
`requirements.sh` instala todo eso correctamente, pero instala **las versiones
vigentes el día de la instalación**, no aquellas con las que la aplicación fue
validada.

La consecuencia es concreta y está declarada como limitación en `VALIDACION.md`:
los informes registran las versiones *a posteriori*, lo que permite **auditar**
un análisis pero no **reconstruir** el entorno en el que se calculó. La imagen
cierra esa brecha: fija R, Bioconductor, los paquetes de CRAN y las herramientas
del pipeline, y deja constancia de todos ellos dentro de la propia imagen.

## Qué fija cada cosa

| Artefacto | Qué fija | Mecanismo |
|---|---|---|
| `bioconductor/bioconductor_docker:RELEASE_3_22` | R 4.5 y todos los paquetes de Bioconductor | Una release de Bioconductor es un conjunto probado entre sí, no versiones sueltas |
| `ARG CRAN_SNAPSHOT` | shiny, bslib, plotly, DT… | Instantánea fechada del repositorio de Posit |
| `docker/environment.yml` | fastqc, fastp, multiqc, bowtie2, samtools, subread, salmon, kallisto | Versiones exactas de bioconda |
| `/opt/env/versions_r.txt` y `/opt/env/versions_cli.txt` | Lo realmente instalado | Se escriben durante la construcción |

Los dos últimos ficheros son el registro autoritativo: si el resolutor de conda
sustituye una versión, o si un paquete opcional no se instala, queda constancia
dentro de la imagen en lugar de descubrirse en ejecución.

## Construir

```bash
docker build -t sara:1.0 .
```

La construcción tarda del orden de decenas de minutos en frío, dominada por la
resolución del entorno conda y la instalación de los paquetes de R. Con la caché
de capas caliente, un cambio en el código de la aplicación reconstruye solo la
última capa: el `COPY` del código va deliberadamente al final del `Dockerfile`.

Para fijar otra instantánea de CRAN:

```bash
docker build --build-arg CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-08-20 -t sara:1.0 .
```

## Ejecutar

```bash
docker compose up
```

La aplicación queda en `http://localhost:3838`. O sin compose:

```bash
docker run --rm -p 3838:3838 -v "$PWD/data:/data:ro" -v "$PWD/outputs:/app/outputs" sara:1.0
```

Dos decisiones de montaje que conviene no cambiar sin motivo:

- **`data/` en solo lectura.** La aplicación no tiene ningún motivo para
  modificar los FASTQ, el genoma o la anotación. Montarlos así elimina la
  posibilidad de que un fallo del pipeline los toque.
- **`outputs/` fuera de la imagen.** Los resultados sobreviven a `docker rm` y a
  reconstruir la imagen. Los datos nunca viven dentro de la imagen: `data/` son
  187 MB, más que todo el software junto.

## Verificar

Una imagen no es válida por construir sin error, sino por comportarse dentro
como lo hacía el entorno de desarrollo:

```bash
docker run --rm sara:1.0 verify
```

Comprueba, en orden de coste creciente, que las ocho herramientas del pipeline
responden, que los paquetes de R que la aplicación da por presentes se cargan, y
que **la batería completa de tests pasa dentro del contenedor**. Devuelve código
de salida distinto de cero si algo falla, de modo que sirve tal cual en
integración continua.

Para consultar el entorno empaquetado:

```bash
docker run --rm sara:1.0 versions
```

## Recursos y paralelismo

`workflow.sh` acepta `--THREADS`. Ese valor debe ajustarse al límite del
**contenedor**, no a los núcleos de la máquina anfitriona: dentro del contenedor
`nproc` sigue viendo los del host, así que pedir más hilos que la cuota asignada
ralentiza el análisis en lugar de acelerarlo. Los límites se declaran en
`compose.yaml` (4 CPU y 16 GB por defecto); la indexación del genoma y el
alineamiento son los pasos que más memoria piden.

## Limitaciones conocidas de la imagen

- **Sin `renv.lock`.** Las versiones de R se fijan por release de Bioconductor y
  por instantánea de CRAN, que es suficiente para reconstruir la imagen, pero no
  produce un fichero de bloqueo utilizable fuera de ella.
- **Sin navegador para `webshot2`.** La descarga de figuras entrega el HTML
  interactivo; la conversión a PNG requiere Chrome/Chromium, que multiplicaría el
  tamaño de la imagen. La aplicación ya contempla esta ausencia y lo explica en
  el fichero que acompaña a la descarga.
- **Arquitectura.** La imagen base y los paquetes de bioconda son `linux/amd64`.
  En equipos Apple Silicon corre bajo emulación, con la penalización de
  rendimiento correspondiente; para uso intensivo conviene construirla en la
  arquitectura de destino.
- **El pipeline sigue sin ser reanudable.** Contenerizar no cambia eso: un fallo
  en la penúltima muestra obliga a repetir las anteriores.
