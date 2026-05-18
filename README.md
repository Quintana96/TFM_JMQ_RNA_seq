# TFM_JMQ_RNA_seq

Aplicación Shiny para ejecutar y monitorizar un workflow de RNA-seq basado en `workflow.sh`.

## Requisitos

- R >= 4.2
- Paquetes R:
  - shiny
  - bslib
  - shinyjs
  - DT
  - ggplot2
  - plotly
  - pheatmap
  - ggrepel
  - RColorBrewer
  - processx
  - DESeq2
  - tximport
  - clusterProfiler
  - enrichplot
- Herramientas de línea de comandos usadas por `workflow.sh`:
  - bowtie2
  - samtools
  - fastqc
  - fastp
  - multiqc
  - featureCounts
  - salmon
  - kallisto

Puedes comprobar e instalar los requerimientos ejecutando:

```bash
./requirements.sh
```

## Ejecutar la app

Desde la raíz del repositorio:

```bash
Rscript -e "shiny::runApp('.')"
```

También puedes ejecutar directamente:

```bash
Rscript app.R
```

## Flujo de trabajo en la app

1. Rellena:
   - Directorio de FASTQ de entrada.
   - Directorio de salida.
   - Archivo de genoma (FASTA).
   - Archivo de anotación (GFF/GTF).
   - Tipo de alineamiento (`bowtie2`, `salmon` o `kallisto`).
   - Tipo de lectura (`paired-end` o `single-end`).
2. Revisa el comando en **Command preview**.
3. Pulsa **Run workflow** para lanzar `workflow.sh`.
4. Consulta la salida en **Run log**.
5. Pulsa **Refresh output status** para listar archivos generados.

`workflow.sh` acepta `--READ_TYPE pe|se`. En modo `pe` espera pares con sufijos
`_1/_2` o `_R1/_R2`; en modo `se` procesa archivos `.fastq` o `.fastq.gz`
individuales y omite archivos con sufijo de R2. Para kallisto single-end se
pueden ajustar `--FRAGMENT_LENGTH` y `--FRAGMENT_SD` desde la App.
