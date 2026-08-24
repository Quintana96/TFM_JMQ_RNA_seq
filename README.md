<p align="center">
  <img src="inst/logo/sara_logo.svg" alt="SARA — Shiny App for RNA-seq Analysis" width="440">
</p>

<p align="center">
  <strong>Del FASTQ al gen diferencial, sin escribir código.</strong><br>
  Interfaz gráfica en R Shiny sobre un pipeline completo de RNA-seq.
</p>

<p align="center">
  <img alt="R" src="https://img.shields.io/badge/R-4.5.2-2E7D5B">
  <img alt="Bioconductor" src="https://img.shields.io/badge/Bioconductor-3.22-2E7D5B">
  <img alt="Shiny" src="https://img.shields.io/badge/Shiny-1.11-2A6F87">
  <img alt="Docker" src="https://img.shields.io/badge/Docker-disponible-2A6F87">
  <img alt="Tests" src="https://img.shields.io/badge/tests-569%20pasan-2E7D5B">
</p>

---

**SARA** (*Shiny App for RNA-seq Analysis*) integra un pipeline bioinformático completo —desde los FASTQ crudos hasta el análisis funcional— en una interfaz gráfica. Ocho herramientas de línea de comandos y una veintena de paquetes de Bioconductor, gobernados desde el navegador, con la trazabilidad que exige un resultado publicable: versiones de cada herramienta, sumas de verificación de las entradas, semillas fijadas y un script de R equivalente que reproduce el análisis fuera de la aplicación.

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
10. [Despliegue con Docker](#despliegue-con-docker)
11. [Validación](#validación)
12. [Notas técnicas](#notas-técnicas)

---

## Descripción general

SARA es una aplicación Shiny que actúa como interfaz gráfica para un pipeline de análisis de RNA-seq. Permite configurar y ejecutar el análisis completo desde la interfaz web, monitorizar el progreso en tiempo real, explorar los resultados de calidad y cuantificación, y realizar análisis de expresión diferencial con múltiples motores estadísticos, todo sin necesidad de escribir código.

La aplicación puede operar en dos modos:

- **Workflow completo**: parte de archivos FASTQ y ejecuta todo el pipeline (QC → trimming → alineamiento → cuantificación).
- **Análisis desde matriz de conteos**: carga directamente una matriz TSV/CSV preexistente y salta al análisis de expresión diferencial.

---

## Características principales

### Tab 0 — Inicio
- Portada con los cuatro pasos del análisis y el estado de cada uno (pendiente, listo, en ejecución, completado).
- Resumen de la sesión: fuente de los conteos, muestras detectadas, matriz en memoria y número de ejecuciones guardadas.
- Cada tarjeta es un acceso directo al paso correspondiente.

### Tab 1 — Configuración
- Selección del modo de inicio: workflow completo o carga de matriz de conteos.
- Elección del tipo de análisis: **alineamiento clásico** (Bowtie2 + featureCounts) o **pseudoalineamiento** (Salmon / Kallisto).
- Soporte para lecturas **paired-end** y **single-end**.
- Selección del directorio de FASTQ con el **diálogo del sistema** —Finder, Explorador o el del escritorio—, o pegando la ruta directamente. No hay raíces ni árbol web que aprender: se navega a cualquier carpeta del equipo, dentro o fuera de la del proyecto. Genoma y anotación usan el diálogo nativo del navegador.
- Cuando no hay escritorio (contenedor o servidor remoto) se cae al selector web de shinyFiles, que sigue funcionando con el servidor en otra máquina.
- Detección automática de muestras y previsualización del número de archivos FASTQ detectados.
- Checklist de validación que verifica todos los parámetros antes de permitir continuar, con la lista de errores concretos junto al botón de continuar.
- En modo "matriz de conteos" se ocultan los campos del pipeline, que en ese modo no aplican.
- Opciones avanzadas del pipeline y calculadora de potencia en un acordeón plegado.

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
- **Coste de la ejecución**: duración y memoria máxima, por paso y en total. Es lo que permite responder a «¿esto cabe en mi portátil o necesito un servidor?» antes de lanzar un conjunto de datos grande, y no después. La memoria se mide muestreando el árbol de procesos completo una vez por segundo, no con `/usr/bin/time`, cuyas opciones y unidades difieren entre BSD y GNU y harían que el contenedor midiese distinto que el equipo de desarrollo.
- Tabla de artefactos de salida con acceso a los archivos generados.
- Resumen MultiQC integrado.
- Control de calidad específico del método (alineamiento o pseudoalineamiento) como pestaña más del mismo navegador de pestañas; solo se muestra el que corresponde a la herramienta de la ejecución.

### Tab 4 — Expresión diferencial
- Tres motores estadísticos seleccionables: **DESeq2**, **edgeR** y **limma-voom**.
- Fuente de datos flexible: ejecución actual, ejecución guardada o matriz subida manualmente.
- Editor de metadatos (samplesheet) con soporte para carga de CSV/TSV o edición directa.
- **Controles en vivo**: el FDR objetivo y el umbral |log2FC| del test se recalculan sobre el ajuste ya hecho, sin relanzar el análisis. No es un recorte de la lista ya calculada: se vuelve a extraer la tabla del modelo, de modo que el filtrado independiente y la hipótesis nula se recalculan como corresponde. Cuesta un 4 % de lo que costó ajustar (0,20 s frente a 5,05 s sobre 20.000 genes × 8 muestras).
- Aviso explícito cuando se cambia un parámetro que sí exige reajustar (motor, diseño, batch, variables sustitutas, prefiltrado o encogido).
- Filtros de visualización independientes: log2 fold-change mínimo, counts mínimos por gen.
- Visualizaciones:
  - Volcano plot interactivo (Plotly).
  - Heatmap de genes diferencialmente expresados (pheatmap).
  - PCA de muestras.
  - MA plot.
- Distribución de la expresión por muestra (densidad y diagrama de cajas), en escala normalizada por composición o sin normalizar, coloreada por condición.
- Análisis de enriquecimiento funcional sobre **GO**, **KEGG**, **Reactome** y conjuntos de genes propios en formato GMT, con ORA y GSEA (GO y Reactome requieren un OrgDb de Bioconductor).
- Seis representaciones del resultado funcional, además del dotplot y el running score: barras, **red gen-concepto** coloreada por log2FC, **mapa de términos** por genes compartidos, **upset** de solapamientos, **ridge** de la distribución del estadístico y el **diagrama oficial de la ruta KEGG** con los cambios pintados encima. Todas descargables en PNG a 300 ppp.
- **Traducción de identificadores con la anotación.** La matriz de conteos trae los identificadores con los que se contaron los genes —locus tags, si el pipeline se lanzó con `featureCounts -g locus_tag`— y ninguna base de anotación funcional los conoce. SARA deduce el atributo de origen comparando la lista contra el GTF y traduce la lista *y el universo* antes del test. Sin ese paso, el enriquecimiento mapea el 0 % y devuelve «sin términos», que no se distingue de la ausencia de señal.
- Descarga de resultados en CSV.
- Los parámetros viven en una barra lateral plegable (acordeón de seis secciones) y los resultados ocupan el área principal, agrupados en seis pestañas: Genes, Muestras, Diagnósticos, Robustez, Enriquecimiento y Reproducibilidad.

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
DESeq2, edgeR, limma, tximport, clusterProfiler, enrichplot, apeglm, fgsea
```

**Bioconductor (opcional):**
```
org.EcK12.eg.db   # Anotación GO/KEGG de E. coli
org.Hs.eg.db      # Anotación humana
ReactomePA        # Enriquecimiento sobre rutas de Reactome
reactome.db       # Conjuntos de Reactome (necesario para el running score)
sva, IHW, qvalue, ashr, rtracklayer, RNASeqPower, ggupset, ggridges, pathview
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


### La interfaz se abre sola

Al arrancar, SARA abre la interfaz en una **ventana de aplicación**: sin barra de direcciones, sin pestañas y sin marcadores. Por dentro sigue siendo el navegador —Chrome, Edge, Brave, Chromium o Vivaldi, el primero que encuentre— en su modo `--app`, pero visualmente no se distingue de una ventana propia.

No es una ventana nativa de verdad. Shiny es una aplicación web y la única forma de tener una ventana nativa real es envolverla en Electron, lo que obliga a empaquetar R entero y convierte el arranque en un proceso de compilación. Para lo que aporta, no compensa.

| Variable | Efecto |
|---|---|
| `SARA_UI=app` | Ventana de aplicación, sin adornos **(por defecto)** |
| `SARA_UI=navegador` | Pestaña normal en el navegador de siempre |
| `SARA_UI=ninguno` | No abre nada; se entra por la URL |
| `SARA_UI_SIZE=1280,800` | Tamaño inicial de la ventana |

Dentro de un contenedor, o en una máquina sin entorno gráfico, no se intenta abrir nada: se detecta y se omite, porque ahí no hay ventana posible y el único efecto sería un error en el log de arranque.

Si no hay ningún navegador compatible, se avisa por consola y se abre en el navegador por defecto: mejor una ventana con barra de direcciones que ninguna.

---

## Guía de uso

### Flujo típico: workflow completo

```
Tab 0 (Inicio) → Tab 1 → Tab 2 → Tab 3 → Tab 4
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
├── Dockerfile               # Imagen de ejecucion autocontenida (dos etapas)
├── compose.yaml             # Despliegue local: puertos, volumenes y limites
├── docker/                  # Todo lo relativo al contenedor
│   ├── environment.yml      #   Herramientas del pipeline, con version fijada
│   ├── install_r_packages.R #   Paquetes de R + registro de versiones
│   ├── entrypoint.sh        #   Comprueba el entorno y arranca la app
│   ├── verify.sh            #   Verifica la imagen (tests dentro del contenedor)
│   └── README.md            #   Construccion, montajes y limitaciones
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
    ├── utils_enrich.R       # Enriquecimiento GO / KEGG / Reactome / GMT (ORA y GSEA)
    ├── utils_diag.R         # Diagnosticos post-ajuste y distribucion de la expresion
    │
    ├── ui_helpers.R         # Helpers UI: botones de descarga, headers con CSV/Plotly
    ├── ui_tab_home.R        # UI del Tab 0 (portada: pasos y estado de la sesión)
    ├── ui_tab_config.R      # UI del Tab 1 (configuración + estado, en dos columnas)
    ├── ui_tab_processing.R  # UI del Tab 2 (activa y versión bloqueada)
    ├── ui_tab_results.R     # UI del Tab 3 (activa y versión vacía)
    ├── ui_tab_deg.R         # UI del Tab 4 (configuración + visualizaciones DEG)
    │
    ├── server_tab_home.R        # Lógica Tab 0: estado de cada paso y navegación
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

**Separación entre ajustar y extraer:**

Los motores devuelven, además de la tabla, el objeto ajustado (`dds`, `glmQLFit`
o `lmFit` según el motor). `deg_reextract()` vuelve a extraer la tabla de ese
objeto con otro FDR o umbral, en lugar de repetir el ajuste. Es lo que permite
que esos dos controles sean reactivos sin caer en el filtro *post hoc* que
invalidaría la FDR declarada. `lfcShrink()` —más caro que el propio ajuste— se
calcula una vez y se reutiliza, porque depende del coeficiente y no del nivel de
significación. La batería de tests comprueba que reextraer da un resultado
**idéntico** a reajustar.

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

## Despliegue con Docker

La instalación nativa (`requirements.sh` + conda) instala **las versiones vigentes
el día de la instalación**, no aquellas con las que la aplicación fue validada:
permite auditar un análisis a posteriori, pero no reconstruir su entorno. La
imagen de Docker cierra esa brecha fijando R, la release de Bioconductor, una
instantánea fechada de CRAN y las versiones exactas de las ocho herramientas del
pipeline.

```bash
docker compose up
```

La aplicación queda en `http://localhost:3838`, con `data/` montado en solo
lectura y `outputs/` fuera de la imagen. Para verificar que la imagen se comporta
como el entorno de desarrollo —herramientas en el PATH, paquetes de R cargables y
**la batería de tests completa dentro del contenedor**:

```bash
docker run --rm sara:1.0 verify
```

Detalles de construcción, montajes, límites de recursos y limitaciones conocidas
en [`docker/README.md`](docker/README.md).

---

## Validación

SARA se ha contrastado contra los resultados publicados de **GSE273773**
(Salazar-Alemán & Turner, 2025, *Sci Rep* 15:1389), reprocesando los FASTQ
originales por las dos rutas del pipeline.

Primero se comprueba que el criterio aplicado es el de los autores: con
`padj < 0,05` y `|log2FC| > 1` sobre *su propia tabla* salen **581 genes
inducidos y 791 reprimidos**, exactamente las cifras del artículo. A partir de
ahí, cualquier diferencia es de datos o de método.

| | Alineamiento (bowtie2) | Pseudoalineamiento (salmon) |
|---|---|---|
| Genes diferenciales publicados recuperados | **96,8 %** | 95,2 % |
| Índice de Jaccard de las listas | 0,912 | 0,888 |
| Coincidencia en el sentido del cambio | **99,9 %** | 99,9 % |
| log2FC · correlación de Pearson | 0,974 | 0,976 |
| log2FC · diferencia absoluta mediana | **0,042** | 0,073 |

Los diez genes más significativos del artículo se recuperan los diez. Las
discrepancias se concentran en el borde de los umbrales y en genes poco
expresados: los que solo detecta SARA tienen un `|log2FC|` mediano de 0,97,
justo por debajo del corte de 1. Hay **una sola discrepancia de sentido en
1.320 genes**, y ambas rutas coinciden entre sí en ella, lo que apunta a la
asignación de lecturas y no a la estadística.

Las dos rutas concuerdan más entre ellas (Pearson 0,9925, Jaccard 0,937) que
cualquiera de las dos con la publicada, que es lo esperable: comparten recorte,
anotación y modelo, y solo difieren en cómo asignan las lecturas.

### Un caso que justifica ofrecer los tres enfoques de enriquecimiento

El artículo destaca homeostasis del hierro, metabolismo del sulfato,
biosíntesis de cisteína y reparación del ADN. A nivel de gen se reproducen sin
discusión (`cysD` +4,11; `cysK` +3,86; `fecA` +2,43; `recA` +1,42;
`sodA` −2,29), pero **el ORA sobre la lista completa no los saca**.

La causa es que esa lista junta 587 genes al alza con 790 a la baja, y la
mezcla cancela cualquier señal direccional:

| Término | Lista completa | Solo al alza |
|---|---|---|
| Homeostasis del hierro | 20/52, esperados 19,3 → p = 0,475 | 18/52, esperados 8,2 → **p = 0,0006** |
| Biosíntesis de cisteína | 7/11 → p = 0,068 | 6/11 → **p = 0,0035** |
| Asimilación de sulfato | 7/9 → p = 0,016 | 6/9 → **p = 0,0008** |

Se recupera de dos maneras, ambas disponibles en la aplicación: marcando el
**ORA direccional**, que parte la lista por sentido, o usando **GSEA**, que
recorre el ranking con signo sin umbralizar (homeostasis del hierro NES +2,22,
p ajustado 1,6e-05).

> Descartadas por medición otras dos explicaciones: no es el tamaño de la lista
> —recortándola de 1.377 a 100 genes el término sigue sin salir— ni la
> traducción de identificadores —el 97,9 % de los diferenciales tiene símbolo
> frente al 92,1 % de los no diferenciales, así que recorta lista y fondo por
> igual—.

---

## Notas técnicas

- **Shiny >= 1.6** recomendado para aprovechar `bindCache`. La app funciona con versiones anteriores sin caché (sin pérdida de funcionalidad, solo de rendimiento).
- **Tamaño máximo de upload**: 10 GB por defecto (`options(shiny.maxRequestSize)`). Ajustable en `global.R`.
- **`org.EcK12.eg.db`** es opcional. Si no está instalado, el análisis de enriquecimiento GO no estará disponible, pero el resto de la app funciona con normalidad.
- **Reactome** requiere `ReactomePA` **y** `reactome.db`; si falta alguno, la colección no aparece en el selector en lugar de aparecer y fallar. Reactome cubre un catálogo cerrado de eucariotas y **no incluye procariotas**: con datos de *E. coli* hay que usar KEGG o un fichero GMT propio. Trabaja en `ENTREZID`, así que los identificadores se traducen internamente y la tasa de mapeo se muestra junto al resultado.
- El pipeline está optimizado por defecto para **8 threads** (`THREADS=8` en `workflow.sh`). Ajustar según los recursos disponibles.
- Los datos de ejemplo incluidos (`BRCA_exp_matrix-1.tsv`, `clinical_info_TCGA-BRCA.tsv`) corresponden al dataset TCGA-BRCA y pueden usarse para probar el flujo de análisis desde matriz de conteos.
