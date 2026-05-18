# TFM_JMQ_RNA_seq

Aplicación Shiny para ejecutar, monitorizar y revisar un workflow de RNA-seq paired-end desde una interfaz gráfica. La app está pensada como punto de entrada operativo para análisis de calidad, trimming, alineamiento o pseudoalineamiento, cuantificación y revisión de resultados sin tener que encadenar manualmente todos los comandos desde terminal.

## Objetivo de la aplicación

`TFM_JMQ_RNA_seq` combina una interfaz Shiny con el script `workflow.sh` para que un usuario pueda:

- Configurar un análisis de RNA-seq a partir de archivos FASTQ paired-end.
- Elegir entre alineamiento clásico con Bowtie2 o pseudoalineamiento con Salmon/Kallisto.
- Ejecutar control de calidad inicial y posterior al trimming.
- Generar cuantificaciones y matrices/resúmenes de conteos.
- Monitorizar el proceso en tiempo real desde la app.
- Revisar ejecuciones guardadas en `outputs/`.
- Consultar tablas, informes MultiQC, logs, gráficos interactivos y alertas interpretativas.
- Descargar resultados clave directamente desde la interfaz.

## Requisitos

### Requisitos de R

- R >= 4.2
- Paquetes CRAN:
  - `shiny`
  - `bslib`
  - `shinyjs`
  - `DT`
  - `ggplot2`
  - `plotly`
  - `pheatmap`
  - `ggrepel`
  - `RColorBrewer`
  - `processx`
  - `shinyFiles`
- Paquetes Bioconductor recomendados:
  - `DESeq2`
  - `tximport`
  - `clusterProfiler`
  - `enrichplot`

### Herramientas de línea de comandos

El workflow usa herramientas externas que deben estar disponibles en el `PATH`:

- `bash`
- `Rscript`
- `fastqc`
- `fastp`
- `multiqc`
- `bowtie2`
- `samtools`
- `featureCounts`
- `salmon`
- `kallisto`

Puedes comprobar e instalar los paquetes R requeridos y verificar las herramientas CLI ejecutando:

```bash
chmod +x requirements.sh workflow.sh
./requirements.sh
```

> Nota: `requirements.sh` instala paquetes de R cuando faltan, pero las herramientas externas deben instalarse con el gestor adecuado para tu sistema, por ejemplo con `conda`, `mamba`, `brew` o `apt`.

## Ejecución

Desde la raíz del repositorio:

```bash
Rscript -e 'shiny::runApp()'
```

También puedes ejecutar directamente:

```bash
Rscript app.R
```

La app crea y consulta por defecto una carpeta `outputs/` dentro del repositorio. Cada ejecución se guarda en un subdirectorio propio con fecha, hora, tipo de análisis y herramienta usada.

## Datos de entrada esperados

### FASTQ

La app detecta muestras paired-end mediante archivos R1 con estos sufijos:

- `_1.fastq.gz`
- `_R1.fastq.gz`
- `_1.fastq`
- `_R1.fastq`

Para cada R1 debe existir un R2 equivalente:

- `_2.fastq.gz`
- `_R2.fastq.gz`
- `_2.fastq`
- `_R2.fastq`

Recomendación: usa nombres de muestra simples, con letras, números, punto, guion o guion bajo. Evita espacios y caracteres especiales.

### Referencia

Según el tipo de análisis, la aplicación solicita:

- Genoma o transcriptoma de referencia en formato FASTA.
- Archivo de anotación GFF/GTF para alineamiento con Bowtie2/featureCounts.
- Anotación opcional para pseudoalineamiento, según el uso previsto.

### Matriz externa de conteos

El modo “Análisis a partir de matriz de conteos” permite cargar una matriz TSV/CSV con genes o transcritos en filas y muestras en columnas. La primera columna debe contener el identificador del gen/transcrito.

## Usabilidad de la interfaz

La interfaz está organizada como un flujo de trabajo progresivo.

### 1 · Configuración

En esta pestaña se prepara el análisis antes de ejecutarlo.

Funciones principales:

- Selección del modo de inicio:
  - Ejecutar workflow completo.
  - Cargar una matriz de conteos externa.
- Selección del tipo de análisis:
  - Alineamiento clásico con Bowtie2 + featureCounts.
  - Pseudoalineamiento con Salmon o Kallisto.
- Selección del directorio de FASTQ mediante selector visual.
- Carga de archivos FASTA y GFF/GTF desde la interfaz.
- Vista resumida de entrada, salida y script de workflow.
- Checklist automático de configuración.
- Detección preliminar de muestras.
- Aviso de pares R2 ausentes o nombres de muestra potencialmente problemáticos.

Diseño de uso recomendado:

1. Selecciona el modo de inicio.
2. Elige alineamiento o pseudoalineamiento.
3. Selecciona el directorio con FASTQ.
4. Carga los archivos de referencia necesarios.
5. Revisa el checklist y las muestras detectadas.
6. Pulsa “Continuar al procesamiento”.

### 2 · Procesamiento

Esta pestaña ejecuta y monitoriza el workflow.

Funciones principales:

- Validación de configuración antes de lanzar el proceso.
- Creación automática de una carpeta de salida única.
- Construcción del comando que se enviará a `workflow.sh`.
- Ejecución integrada mediante `processx` cuando está disponible.
- Registro del log en tiempo real.
- Barra/estado de progreso por checkpoints.
- Seguimiento del procesamiento por muestra.
- Detección del código de salida del proceso.
- Carga automática de la matriz de conteos cuando el workflow termina correctamente.
- Redirección automática a Resultados al finalizar.

El log completo se guarda como `workflow_live.log` dentro de la carpeta de salida.

### 3 · Resultados

Esta pestaña centraliza la revisión de ejecuciones terminadas o guardadas.

Funciones principales:

- Selector de ejecuciones guardadas en `outputs/`.
- Opción para seleccionar manualmente una carpeta de resultados.
- Resumen visual del estado de la ejecución:
  - completado
  - error
  - incompleto
  - sin log
- Métricas rápidas:
  - número de muestras
  - número de genes/transcritos detectados
  - mapeo medio
  - tamaño y número de archivos generados
- Apertura directa del informe `multiqc_report.html` cuando existe.
- Tablas interactivas con filtrado, paginación y desplazamiento horizontal.
- Descarga de tablas y artefactos principales.
- Visualización de las últimas líneas del log.

Subsecciones disponibles:

- **Resumen**: interpretación rápida y estadísticas principales de MultiQC.
- **Calidad**: tablas FastQC, trimming y métricas de cuantificación/alineamiento.
- **Conteos**: resumen por muestra y genes/transcritos más abundantes.
- **Informes y archivos**: archivos principales, rutas y listado completo.
- **Log**: últimas líneas del workflow.

### Guía integrada

La app incluye una pestaña “Guía” con explicación interna del propósito de la aplicación, flujo de uso, características principales, buenas prácticas y estructura de resultados. Esta pestaña permite usar la app como punto de partida oficial sin depender únicamente de documentación externa.

## Características técnicas principales

### Control de calidad

- FastQC sobre FASTQ originales.
- fastp para trimming, filtrado y detección de adaptadores paired-end.
- FastQC sobre lecturas trimadas.
- MultiQC para consolidar métricas de calidad, trimming y alineamiento/cuantificación.

### Alineamiento clásico

Cuando se selecciona alineamiento:

1. Se construye o reutiliza un índice Bowtie2.
2. Se alinean lecturas paired-end contra la referencia.
3. Se ordenan los BAM con `samtools sort`.
4. Se indexan los BAM con `samtools index`.
5. Se cuantifican lecturas por gen con `featureCounts`.
6. Se genera una matriz `04_counts/count_matrix.tsv`.

### Pseudoalineamiento

Cuando se selecciona pseudoalineamiento:

- Salmon:
  - Construcción/reutilización de índice Salmon.
  - Ejecución de `salmon quant` por muestra.
  - Lectura de archivos `quant.sf`.
- Kallisto:
  - Construcción/reutilización de índice Kallisto.
  - Ejecución de `kallisto quant` por muestra.
  - Lectura de archivos `abundance.tsv`.

La app intenta cargar los conteos con `tximport` cuando está disponible y, si no, usa una lectura directa de los archivos de cuantificación.

### Resultados adicionales de QC

La interfaz añade una capa interpretativa adicional para facilitar la revisión de calidad.

Para alineamiento:

- resumen de lecturas totales
- lecturas únicas
- lecturas multimapeadas
- lecturas no alineadas
- tasa de mapeo
- tasa de multimapeo
- lecturas asignadas/no asignadas
- alertas automáticas por umbrales configurables

Para pseudoalineamiento:

- fragmentos procesados
- lecturas pseudoalineadas
- pseudoalignment rate
- transcritos detectados
- fracción de TPM cercano a cero
- distribución de TPM por muestra
- relación TPM vs NumReads
- alertas automáticas por umbrales configurables

Las alertas son orientativas y deben interpretarse como cribado inicial de calidad, no como diagnóstico definitivo.

## Estructura de salida

Cada ejecución genera una carpeta dentro de `outputs/`. La estructura esperada es:

```text
outputs/
└── <fecha>_<hora>_<tipo>_<herramienta>/
    ├── 01_quality/
    ├── 02_trimmed_reads/
    ├── 03_alignments/
    │   └── bowtie2|salmon|kallisto/
    ├── 04_counts/
    ├── indices/
    ├── multiqc_data/
    ├── multiqc_report.html
    └── workflow_live.log
```

Archivos especialmente relevantes:

- `multiqc_report.html`: informe global de calidad y ejecución.
- `workflow_live.log`: log completo del análisis.
- `04_counts/count_matrix.tsv`: matriz de conteos cuando está disponible.
- `02_trimmed_reads/*_fastp.html`: informes individuales de trimming.
- `01_quality/*_fastqc.html`: informes FastQC.
- `03_alignments/<tool>/`: BAM o cuantificaciones por muestra.

## Descargas desde la app

La interfaz permite descargar:

- Estadísticas principales MultiQC.
- Estado FastQC.
- Métricas de trimming y alineamiento/cuantificación.
- Librerías por muestra.
- Genes/transcritos más abundantes.
- Informes principales.
- Listado completo de archivos generados.
- Log de ejecución.
- Tablas de QC adicional.
- Gráficos interactivos exportados como HTML.

## Ejecución directa del workflow por terminal

Aunque la app es la vía recomendada, el workflow puede ejecutarse manualmente:

```bash
./workflow.sh \
  --INPUT /ruta/a/fastq \
  --OUTPUT /ruta/a/salida \
  --GENOME_FILE /ruta/referencia.fasta \
  --ANNOTATION_FILE /ruta/anotacion.gff \
  --ALIGNMENT_TYPE bowtie2
```

Valores permitidos para `--ALIGNMENT_TYPE`:

- `bowtie2`
- `salmon`
- `kallisto`

## Buenas prácticas

- Verifica los requerimientos antes de ejecutar análisis largos.
- Usa rutas sin espacios si tu entorno de terminal es problemático.
- Mantén juntas las parejas R1/R2 de cada muestra.
- Usa nombres de muestra consistentes.
- Revisa el checklist antes de iniciar el procesamiento.
- Conserva la carpeta completa de salida si necesitas auditar el análisis.
- Consulta el log si una ejecución aparece como incompleta o con error.
- No elimines `multiqc_data/` si quieres que las tablas de resumen se carguen correctamente.

## Limitaciones conocidas

- El workflow está orientado a datos paired-end.
- Algunas métricas dependen de que MultiQC genere los archivos esperados en `multiqc_data/`.
- Las alertas automáticas usan umbrales generales; pueden requerir ajuste según organismo, profundidad, protocolo o diseño experimental.
- La ejecución de análisis grandes puede consumir mucho tiempo, RAM y espacio en disco.
- La app no reemplaza la interpretación estadística/biológica experta; organiza resultados y facilita la revisión técnica.

## Comprobaciones rápidas para desarrollo

Comprobar sintaxis del script Bash:

```bash
bash -n workflow.sh
```

Lanzar la app desde la raíz del repositorio:

```bash
Rscript app.R
```
