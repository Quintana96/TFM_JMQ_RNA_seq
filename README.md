# TFM_JMQ_RNA_seq

Aplicación Shiny para ejecutar y monitorizar un workflow de RNA-seq basado en `workflow.sh`.

## Requisitos

- R >= 4.2
- Paquetes R: `shiny`, `bslib`
- Herramientas de línea de comandos usadas por `workflow.sh` (bowtie2/salmon/kallisto, fastqc, fastp, samtools, featureCounts, multiqc)

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
2. Revisa el comando en **Command preview**.
3. Pulsa **Run workflow** para lanzar `workflow.sh`.
4. Consulta la salida en **Run log**.
5. Pulsa **Refresh output status** para listar archivos generados.
