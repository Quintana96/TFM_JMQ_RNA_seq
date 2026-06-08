# RNA-seq Workflow Runner (Nombre provisional)

Una aplicación web interactiva desarrollada con R Shiny para el análisis integral de datos de RNA-seq. Integra un pipeline bioinformático completo —desde los archivos FASTQ crudos hasta el análisis de expresión diferencial— en una interfaz gráfica accesible.

> Desarrollada como Trabajo de Fin de Máster (TFM) en la Universidad Europea de Madrid.

---

## Tabla de contenidos

1. [Descripción general](#descripción-general)
2. [Características principales](#características-principales)
3. [Requisitos del sistema](#requisitos-del-sistema)
4. [Instalación](#instalación)
5. [Cómo ejecutar](#cómo-ejecutar)
6. [Guía de uso](#guía-de-uso)
7. [Estructura del proyecto](#estructura-del-proyecto)
8. [Pipeline bioinformático](#pipeline-bioinformático)
9. [Arquitectura de la aplicación](#arquitectura-de-la-aplicación)
10. [Notas técnicas](#notas-técnicas)

---

## Descripción general

RNA-seq Workflow Runner es una aplicación Shiny que actúa como interfaz gráfica para un pipeline de análisis de RNA-seq. Permite configurar y ejecutar el análisis completo desde la interfaz web, monitorizar el progreso en tiempo real, explorar los resultados de calidad y cuantificación, y realizar análisis de expresión diferencial con múltiples motores estadísticos, todo sin necesidad de escribir código.

La aplicación puede operar en dos modos:

- **Workflow completo**: parte de archivos FASTQ y ejecuta todo el pipeline (QC → trimming → alineamiento → cuantificación).
- **Análisis desde matriz de conteos**: carga directamente una matriz TSV/CSV preexistente y salta al análisis de expresión diferencial.

---

## Características principales

### Tab 1 — Configuración
- Selección del modo de inicio: workflow completo o carga de matriz de conteos.
- Elección del tipo de análisis: **alineamiento clásico** (Bowtie2 + featureCounts) o **pseudoalineamiento** (Salmon / Kallisto).
- Soporte para lecturas **paired-end** y **single-end**.
- Selección interactiva de directorios para archivos FASTQ de entrada, genoma de referencia y anotación.
- Detección automática de muestras y previsualización del número de archivos FASTQ detectados.
- Checklist de validación que verifica todos los parámetros antes de permitir continuar.

### Tab 2 — Procesamiento
- Ejecución del pipeline `workflow.sh` directamente desde la interfaz.
- Log en tiempo real con estado del proceso (running / success / error).
- Indicador de progreso y tiempo transcurrido.
- Botón de cancelación del proceso en curso.
- Guardado automático del log completo en disco (`workflow_live.log`).

### Tab 3 — Resultados
- Exploración de todas las ejecuciones previas guardadas en el directorio `outputs/`.
- Resumen de métricas de calidad con alertas automáticas basadas en umbrales configurables:
  - Tasa de alineamiento / pseudoalineamiento.
  - Tasa de asignación de features (featureCounts).
  - Multimapping rate.
  - Distribución de reads y genes detectados.
- Tablas interactivas de conteos (genes × muestras) con opción de descarga en CSV.
- Tabla de artefactos de salida con acceso a los archivos generados.
- Resumen MultiQC integrado.

### Tab 4 — Expresión diferencial
- Tres motores estadísticos seleccionables: **DESeq2**, **edgeR** y **limma-voom**.
- Fuente de datos flexible: ejecución actual, ejecución guardada o matriz subida manualmente.
- Editor de metadatos (samplesheet) con soporte para carga de CSV/TSV o edición directa.
- Filtros interactivos: log2 fold-change mínimo, p-valor ajustado máximo, counts mínimos por gen.
- Visualizaciones:
  - Volcano plot interactivo (Plotly).
  - Heatmap de genes diferencialmente expresados (pheatmap).
  - PCA de muestras.
  - MA plot.
- Análisis de enriquecimiento funcional GO (requiere OrgDb de Bioconductor).
- Descarga de resultados en CSV.

---

## Requisitos del sistema

### Herramientas de sistema (CLI)

| Herramienta | Uso |
|---|---|
| `R` ≥ 4.1 | Entorno de ejecución de la aplicación |
| `bash` | Ejecución del pipeline |
| `fastqc` | Control de calidad de reads crudos |
| `fastp` | Trimming y filtrado de reads |
| `multiqc` | Reporte agregado de QC |
| `bowtie2` | Alineamiento clásico (modo alignment) |
| `samtools` | Procesamiento de BAM (modo alignment) |
| `featureCounts` (subread) | Cuantificación (modo alignment) |
| `salmon` | Pseudoalineamiento (modo pseudo) |
| `kallisto` | Pseudoalineamiento alternativo (modo pseudo) |

> Se requieren al menos las herramientas del modo de análisis que se vaya a usar. No es necesario tener instalados todos los alineadores a la vez.

### Paquetes R requeridos

**CRAN:**
```
shiny, bslib, shinyjs, DT, ggplot2, plotly, pheatmap, ggrepel, RColorBrewer, processx
```

**Bioconductor:**
```
DESeq2, edgeR, limma, tximport, clusterProfiler, enrichplot
```

**Bioconductor (opcional):**
```
org.EcK12.eg.db  # Para enriquecimiento GO en E. coli
```

---

## Instalación

### 1. Instalar herramientas de sistema

Con [conda/bioconda](https://bioconda.github.io/) (recomendado):

```bash
conda create -n rnaseq_app -c conda-forge -c bioconda \
  bowtie2 samtools fastqc fastp subread multiqc salmon kallisto
conda activate rnaseq_app
```

### 2. Instalar dependencias R

Ejecuta el script incluido, que comprueba e instala automáticamente todo lo necesario:

```bash
bash requirements.sh
```

El script verifica la disponibilidad de los programas CLI, instala los paquetes CRAN y Bioconductor que falten, y confirma que `workflow.sh` está presente.

---

## Cómo ejecutar

Desde el directorio raíz del proyecto (donde se encuentra `app.R`):

```r
# Opción 1: desde R o RStudio
shiny::runApp(".")
```

```bash
# Opción 2: desde la terminal
Rscript app.R
```

La aplicación se abrirá en el navegador por defecto. Por defecto acepta uploads de hasta **10 GB** para soportar archivos FASTQ grandes.

---

## Guía de uso

### Flujo típico: workflow completo

```
Tab 1 → Tab 2 → Tab 3 → Tab 4
Configurar  Ejecutar   Revisar   Analizar DEG
```

**1. Configurar** (`Tab 1 · Configuración`)

1. Selecciona **"Ejecutar workflow completo"** en el modo de inicio.
2. Elige el **tipo de análisis**: Alineamiento (Bowtie2) o Pseudoalineamiento (Salmon/Kallisto).
3. Indica si tus datos son **paired-end** o **single-end**.
4. Selecciona el **directorio de FASTQ** de entrada — la app detecta automáticamente las muestras.
5. Selecciona el **archivo de genoma** (FASTA) y el **archivo de anotación** (GTF/GFF).
6. Elige el **directorio de salida** donde se guardarán los resultados.
7. Verifica que el checklist de validación esté completo y haz clic en **"Continuar al procesamiento"**.

**2. Procesar** (`Tab 2 · Procesamiento`)

1. Revisa el resumen de la configuración y el número de muestras detectadas.
2. Haz clic en **"Ejecutar pipeline"**.
3. Monitoriza el log en tiempo real; el panel indica el estado (running / success / error).
4. Al finalizar con éxito, la app navega automáticamente a la pestaña de Resultados.

**3. Revisar resultados** (`Tab 3 · Resultados`)

1. Selecciona una ejecución del selector desplegable (se listan todas las carpetas en `outputs/`).
2. Explora las métricas de QC y las alertas automáticas.
3. Revisa la tabla de conteos y descárgala en CSV si es necesario.
4. Consulta el resumen MultiQC y los artefactos generados.

**4. Expresión diferencial** (`Tab 4 · Expresión diferencial`)

1. Selecciona la fuente de datos (ejecución actual, guardada o matriz subida).
2. Sube o edita el **samplesheet** con los metadatos (columna `sample` + columna de condición).
3. Configura la comparación: condición de referencia vs condición de tratamiento.
4. Selecciona el motor estadístico: DESeq2, edgeR o limma-voom.
5. Ajusta los filtros (log2FC, p-adj, counts mínimos) y lanza el análisis.
6. Explora los gráficos (volcano, heatmap, PCA) y descarga los resultados.

### Flujo alternativo: análisis desde matriz de conteos

Si ya tienes una matriz de conteos preexistente (TSV/CSV, genes como filas, muestras como columnas):

1. En `Tab 1`, selecciona **"Análisis a partir de matriz de conteos"**.
2. Sube el archivo y haz clic en **"Análisis a partir de matriz de conteos"**.
3. La app pasa directamente al `Tab 4 · Expresión diferencial`.

---

## Estructura del proyecto

```
rnaseq_app/
├── app.R                    # Entrypoint: carga archivos y arranca shinyApp()
├── global.R                 # Librerías, constantes, tema Bootstrap 5 y CSS
├── ui.R                     # Navbar de 4 pestañas (page_navbar)
├── server.R                 # Server principal: instancia state y delega en módulos
├── workflow.sh              # Pipeline Bash (QC → trimming → índice → alineamiento → conteos)
├── requirements.sh          # Comprueba e instala todas las dependencias
├── BRCA_exp_matrix-1.tsv    # Matriz de ejemplo (datos TCGA-BRCA)
├── clinical_info_TCGA-BRCA.tsv # Metadatos clínicos de ejemplo
└── R/                       # Módulos auto-sourced por Shiny
    ├── state.R              # create_app_state(): reactivos compartidos entre tabs
    │
    ├── utils_format.R       # Formateo: bytes, tiempo, porcentajes, logs
    ├── utils_fastq.R        # Detección de muestras FASTQ, validación de pares R1/R2
    ├── utils_io.R           # Lectura de archivos, directorios de resultados, tail de logs
    ├── utils_counts.R       # Carga de matrices de conteos (tximport, TSV, workflow)
    ├── utils_multiqc.R      # Parseo de estadísticas MultiQC
    ├── utils_status.R       # Estado de ejecución desde log, badges de status
    ├── utils_tables.R       # Tablas DT: FASTQ, alineamiento, conteos, artefactos
    ├── utils_qc.R           # Métricas QC, alertas, gráficos Plotly de calidad
    ├── utils_deg.R          # Análisis DEG: DESeq2 / edgeR / limma-voom
    ├── utils_enrich.R       # Enriquecimiento funcional GO (clusterProfiler)
    │
    ├── ui_helpers.R         # Helpers UI: botones de descarga, headers con CSV/Plotly
    ├── ui_tab_config.R      # UI del Tab 1 (6 cards en grid 3×2)
    ├── ui_tab_processing.R  # UI del Tab 2 (activa y versión bloqueada)
    ├── ui_tab_results.R     # UI del Tab 3 (activa y versión vacía)
    ├── ui_tab_deg.R         # UI del Tab 4 (configuración + visualizaciones DEG)
    │
    ├── server_tab_config.R      # Lógica Tab 1: validación, detección de muestras
    ├── server_tab_processing.R  # Lógica Tab 2: ejecución del pipeline, log en vivo
    ├── server_tab_results.R     # Lógica Tab 3: carga y renderizado de resultados
    └── server_tab_deg.R         # Lógica Tab 4: DEG, gráficos, enriquecimiento
```

---

## Pipeline bioinformático

El script `workflow.sh` implementa el siguiente flujo de análisis:

```
FASTQ (crudos)
      │
      ▼
  FastQC          ← Control de calidad inicial
      │
      ▼
  fastp           ← Trimming de adaptadores y filtrado por calidad
      │
      ▼
  FastQC          ← Control de calidad post-trimming
      │
  ┌───┴───────────────────┐
  │                       │
  ▼                       ▼
Bowtie2              Salmon / Kallisto
(alineamiento)       (pseudoalineamiento)
  │                       │
  ▼                       ▼
samtools sort/index   cuantificación directa
  │
  ▼
featureCounts
  │
  └───────────────────────┤
                          ▼
                       MultiQC    ← Reporte agregado de QC
                          │
                          ▼
                  Matriz de conteos
```

**Parámetros del pipeline:**

| Flag | Descripción | Valores |
|---|---|---|
| `--INPUT` | Directorio con archivos FASTQ | ruta absoluta |
| `--OUTPUT` | Directorio de salida | ruta absoluta |
| `--GENOME_FILE` | Genoma de referencia (FASTA) | ruta al archivo |
| `--ANNOTATION_FILE` | Anotación génica (GTF/GFF) | ruta al archivo |
| `--ALIGNMENT_TYPE` | Motor de alineamiento | `bowtie2` (defecto), `salmon`, `kallisto` |
| `--READ_TYPE` | Tipo de lectura | `pe` (defecto), `se` |
| `--FRAGMENT_LENGTH` | Longitud media del fragmento (SE + kallisto) | número (defecto: 200) |
| `--FRAGMENT_SD` | Desviación estándar del fragmento (SE + kallisto) | número (defecto: 20) |

**Estructura de salida generada:**

```
outputs/<nombre_ejecucion>/
├── 01_quality/          # FastQC antes y después de trimming
├── 02_trimmed_reads/    # FASTQ filtrados por fastp
├── 03_alignments/       # BAM ordenados e indexados (modo bowtie2)
│   └── bowtie2/
├── 04_counts/           # Matrices de conteos finales
├── indices/             # Índices de Bowtie2 / Salmon / Kallisto
└── workflow_live.log    # Log completo de la ejecución
```

---

## Arquitectura de la aplicación

La aplicación sigue una arquitectura modular multi-archivo estándar de Shiny:

- **No se usan `moduleServer()` / `NS()`**: las cuatro pestañas comparten un objeto `state` con todos los `reactiveVal` / `reactiveValues` comunes. Esto evita renombrar los IDs originales y simplifica la comunicación entre tabs.
- **`state` se crea en `create_app_state(session)`** (`R/state.R`) y expone rutas del sistema de archivos, reactivos de configuración validada, resultados cargados y el slot `state$shared` con los reactivos derivados del Tab 1 que necesitan los demás tabs.
- **`R/` se sourcea automáticamente**: Shiny carga todos los archivos `*.R` del directorio antes de `ui.R` y `server.R`. No hay llamadas `source()` manuales (excepto en el modo `Rscript app.R`).

**Optimizaciones de rendimiento:**
- `debounce` de 600 ms sobre el input del directorio FASTQ para evitar reactivaciones en cascada mientras el usuario escribe.
- `bindCache` en los reactivos pesados de Tab 3 (clave `dir|timestamp`), con fallback automático para versiones de Shiny < 1.6.
- Lectura tail-only del log (`read_tail_text`) para no cargar el archivo entero en cada actualización.
- Detección de muestras centralizada en un único reactivo `samples_eff()` reutilizado en validación, checklist y observer del botón de procesamiento.

**Umbrales de alerta QC** (configurables en `global.R`):

| Métrica | Warning | Error |
|---|---|---|
| Tasa de alineamiento | < 70 % | < 50 % |
| Tasa de asignación (featureCounts) | < 60 % | — |
| Multimapping rate | > 20 % | — |
| Tasa de pseudoalineamiento | < 70 % | — |

---

## Notas técnicas

- **Shiny >= 1.6** recomendado para aprovechar `bindCache`. La app funciona con versiones anteriores sin caché (sin pérdida de funcionalidad, solo de rendimiento).
- **Tamaño máximo de upload**: 10 GB por defecto (`options(shiny.maxRequestSize)`). Ajustable en `global.R`.
- **`org.EcK12.eg.db`** es opcional. Si no está instalado, el análisis de enriquecimiento GO no estará disponible, pero el resto de la app funciona con normalidad.
- El pipeline está optimizado por defecto para **8 threads** (`THREADS=8` en `workflow.sh`). Ajustar según los recursos disponibles.
- Los datos de ejemplo incluidos (`BRCA_exp_matrix-1.tsv`, `clinical_info_TCGA-BRCA.tsv`) corresponden al dataset TCGA-BRCA y pueden usarse para probar el flujo de análisis desde matriz de conteos.
