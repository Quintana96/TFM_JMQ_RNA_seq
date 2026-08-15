#' global.R
#' Librerias, opciones globales, constantes y tema de la aplicacion.
#' Se sourcea automaticamente por Shiny antes de app.R / ui.R / server.R.

library(shiny)
library(bslib)
library(shinyjs)
library(DT)
library(shinyFiles)

# Permitir uploads grandes (~10 GB) y cache de sass en tempdir
options(shiny.maxRequestSize = 10 * 1024^3)
options(sass.cache = file.path(tempdir(), "rnaseq_shiny_sass_cache"))

#' Detector simple de paquete instalado
pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)

HAS_PROCESSX <- pkg_ok("processx")
HAS_TXIMPORT <- pkg_ok("tximport")

# Paquetes para el modulo DEG (Tab 4)
HAS_DESEQ2          <- pkg_ok("DESeq2")
HAS_EDGER           <- pkg_ok("edgeR")
HAS_LIMMA           <- pkg_ok("limma")
HAS_CLUSTERPROFILER <- pkg_ok("clusterProfiler")
HAS_PHEATMAP        <- pkg_ok("pheatmap")
HAS_ORGECDB         <- pkg_ok("org.EcK12.eg.db")
# Encogido de log2FC (lfcShrink). apeglm es la opcion recomendada por la
# vinieta de DESeq2; ashr sirve de alternativa. Sin ninguno de los dos se cae a
# type = "normal".
HAS_APEGLM          <- pkg_ok("apeglm")
HAS_ASHR            <- pkg_ok("ashr")
# Fase 2: GSEA, ponderacion de hipotesis y segunda estimacion de pi0. La interfaz
# oculta las opciones cuyo paquete falta, en lugar de ofrecerlas y fallar.
HAS_FGSEA           <- pkg_ok("fgsea")
HAS_IHW             <- pkg_ok("IHW")
HAS_QVALUE          <- pkg_ok("qvalue")
# Fase 3: variacion no deseada y lectura de la anotacion. rtracklayer es lo que
# permite construir el mapa transcrito-gen que necesita tximport.
HAS_SVA             <- pkg_ok("sva")
HAS_RTRACKLAYER     <- pkg_ok("rtracklayer")
# Fase 4: motores robustos, Swish y calculo de potencia.
HAS_DEARSEQ         <- pkg_ok("dearseq")
HAS_FISHPOND        <- pkg_ok("fishpond")
HAS_RNASEQPOWER     <- pkg_ok("RNASeqPower")

#' Operador nulo-coalescente
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


#' @section Determinismo
#' Semilla unica de todo lo estocastico del analisis: el numero de variables
#' sustitutas de `sva::num.sv` (metodo "be", por permutacion), las permutaciones
#' de `fishpond::swish`, el remuestreo del panel de replicabilidad y las
#' miniaturas simuladas de los diagnosticos.
#'
#' Se centraliza para que el informe pueda declararla: una semilla fija que no
#' se registra no sirve de nada, porque quien reproduzca el analisis no sabe
#' cual era. Todos los usos van por `withr::with_seed()`, que restaura el estado
#' del RNG al salir y evita que fijar la semilla en un sitio altere la
#' aleatoriedad de cualquier otro.
ANALYSIS_SEED <- 1L

#' Replicas inferenciales que permuta Swish. Es un parametro del resultado, no
#' un detalle interno, asi que viaja al informe junto a la semilla.
SWISH_NPERMS <- 30L


#' @section Constantes FASTQ
FASTQ_R1_SUFFIXES <- c("_1.fastq.gz", "_R1.fastq.gz", "_1.fastq", "_R1.fastq")
FASTQ_R2_SUFFIXES <- c("_2.fastq.gz", "_R2.fastq.gz", "_2.fastq", "_R2.fastq")
FASTQ_R1_PATTERN  <- "(_1\\.fastq\\.gz$|_R1\\.fastq\\.gz$|_1\\.fastq$|_R1\\.fastq$)"
FASTQ_R2_PATTERN  <- "(_2\\.fastq\\.gz$|_R2\\.fastq\\.gz$|_2\\.fastq$|_R2\\.fastq$)"
FASTQ_ANY_PATTERN <- "(\\.fastq\\.gz$|\\.fastq$)"


#' @section Umbrales para alertas de QC adicional
qc_thresholds <- list(
  mapping_rate_warning            = 0.70,
  mapping_rate_error              = 0.50,
  assigned_rate_warning           = 0.60,
  multimapping_rate_warning       = 0.20,
  pseudoalignment_rate_warning    = 0.70,
  low_reads_fraction_warning      = 0.50,
  low_detected_fraction_warning   = 0.50,
  near_zero_tpm_fraction_warning  = 0.80,
  tpm_distribution_shift_warning  = 2
)


#' @section Tema y CSS
app_theme <- bs_theme(
  version = 5,
  primary = "#BEE8C8",
  secondary = "#D8F1DD",
  success = "#A8DDB8",
  info = "#A8DADC",
  warning = "#F6D58A",
  danger = "#F4A6A6",
  base_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Arial", "sans-serif"),
  heading_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Arial", "sans-serif"),
  bg = "#F6FBF7",
  fg = "#20332A"
)

app_css <- HTML("
  :root {
    --pastel-green: #D8F1DD;
    --pastel-green-strong: #A8DDB8;
    --pastel-green-soft: #F1FAF3;
    --pastel-mint: #EAF8EE;
    --pastel-blue: #D7EEF1;
    --pastel-lavender: #E7DDF5;
    --pastel-peach: #F9DCC4;
    --pastel-rose: #FADDE1;
    --pastel-yellow: #FFF1C9;
    --pastel-coral: #FAD1CF;
    --ink: #20332A;
    --muted: #60756A;
    --line: #5F6F66;
    --thin-border: 0.75px solid #5F6F66;
  }

  body,
  .bslib-page-navbar {
    background:
      linear-gradient(180deg, #FAFEFB 0%, #F1FAF3 48%, #FBFEFC 100%);
    color: var(--ink);
  }

  .navbar {
    background: linear-gradient(90deg, #D7EEF1 0%, #E7DDF5 56%, #F9DCC4 100%) !important;
    border: var(--thin-border);
    border-width: 0 0 1px 0;
    box-shadow: 0 8px 24px rgba(93, 137, 108, 0.12);
  }

  .navbar-brand,
  .navbar .nav-link {
    color: var(--ink) !important;
    font-weight: 600;
    letter-spacing: 0;
  }

  .navbar .nav-link.active {
    background: rgba(255, 255, 255, 0.62);
    border: var(--thin-border);
    border-radius: 8px;
  }

  .card {
    border: var(--thin-border);
    border-radius: 8px;
    background: rgba(255, 255, 255, 0.92);
    box-shadow: 0 10px 24px rgba(94, 128, 105, 0.10);
  }

  .card-header {
    background: linear-gradient(90deg, #E7DDF5 0%, #F8FCF9 100%);
    border-bottom: var(--thin-border);
    color: var(--ink);
    font-weight: 600;
  }

  .card:nth-of-type(3n + 1) > .card-header {
    background: linear-gradient(90deg, #D7EEF1 0%, #F8FCF9 100%);
  }

  .card:nth-of-type(3n + 2) > .card-header {
    background: linear-gradient(90deg, #F9DCC4 0%, #F8FCF9 100%);
  }

  .card:nth-of-type(3n) > .card-header {
    background: linear-gradient(90deg, #E7DDF5 0%, #F8FCF9 100%);
  }

  .btn-primary,
  .btn-success {
    background: #A8DDB8;
    border: var(--thin-border);
    color: #173126;
    font-weight: 650;
  }

  .btn-primary:hover,
  .btn-success:hover {
    background: #95D2A8;
    border-color: #5F6F66;
    color: #12261D;
  }

  .btn-secondary,
  .btn-outline-secondary {
    background: #FFFFFF;
    border: var(--thin-border);
    color: #315342;
  }

  .btn-danger {
    background: #F4A6A6;
    border: var(--thin-border);
    color: #4A1F1F;
  }

  .btn-continue-blue {
    background: #A8DADC;
    border: var(--thin-border);
    color: #1F454B;
    font-weight: 600;
    box-shadow: 0 8px 18px rgba(80, 140, 150, 0.18);
  }

  .btn-continue-blue:hover,
  .btn-continue-blue:focus {
    background: #93CDD3;
    border-color: #5F6F66;
    color: #173A40;
  }

  .form-control,
  .form-select,
  .selectize-input {
    border: var(--thin-border) !important;
    border-radius: 8px;
  }

  /* Botones de seleccion de ficheros y directorios (shinyFiles, fileInput).
     Antes esta regla usaba el selector button[id$='_btn'], que por
     especificidad (0,1,1 frente a 0,1,0) VENCIA a .btn-success y .btn-danger:
     como todos los botones de accion terminan en _btn, el de ejecutar el
     workflow (verde) y el de detener (rojo) se pintaban con este azul de
     selector de ficheros, y se perdia la semantica de color de la accion
     principal y de la destructiva. Ahora se marcan con .btn-picker. */
  .btn-default,
  .btn-file,
  .input-group .btn,
  .btn-picker {
    background: #D7EEF1;
    border: var(--thin-border);
    color: #244A50;
    font-weight: 500;
  }

  .btn-default:hover,
  .btn-file:hover,
  .input-group .btn:hover,
  .btn-picker:hover {
    background: #C6E5EA;
    border-color: #5F6F66;
    color: #1F454B;
  }

  .input-group .form-control {
    background: rgba(255, 255, 255, 0.94);
    color: #20332A;
  }

  .nav-tabs .nav-link.active {
    background: #E7DDF5;
    border-color: #5F6F66 #5F6F66 #E7DDF5;
    color: #20332A;
    font-weight: 600;
  }

  .nav-tabs .nav-link {
    color: #4B3F61;
    font-weight: 500;
  }

  .alert-info {
    background: #D7EEF1;
    border: var(--thin-border);
    color: #244A50;
  }

  .alert-success {
    background: #EAF8EE;
    border: var(--thin-border);
    color: #244B34;
  }

  .alert-warning {
    background: #FFF1C9;
    border: var(--thin-border);
    color: #5C4A16;
  }

  .alert-danger {
    background: #FAD1CF;
    border: var(--thin-border);
    color: #5A2323;
  }

  .metric-card {
    padding: 12px;
    border-radius: 8px;
    min-height: 76px;
    display: flex;
    flex-direction: column;
    justify-content: center;
  }

  .metric-card:nth-child(1) {
    background: linear-gradient(180deg, #FFFFFF 0%, #E7DDF5 100%);
    border: var(--thin-border);
  }

  .metric-card:nth-child(2) {
    background: linear-gradient(180deg, #FFFFFF 0%, #D7EEF1 100%);
    border: var(--thin-border);
  }

  .metric-card:nth-child(3) {
    background: linear-gradient(180deg, #FFFFFF 0%, #F9DCC4 100%);
    border: var(--thin-border);
  }

  .metric-card:nth-child(4) {
    background: linear-gradient(180deg, #FFFFFF 0%, #FADDE1 100%);
    border: var(--thin-border);
  }

  .metric-card-value {
    font-size: 18px;
    font-weight: 650;
    color: #20332A;
  }

  .metric-card-label {
    font-size: 11px;
    color: #60756A;
    margin-top: 6px;
    text-transform: uppercase;
    letter-spacing: .04em;
  }

  table.dataTable thead th {
    background: #E7DDF5;
    color: #20332A;
    font-weight: 600;
    border-bottom: var(--thin-border) !important;
  }

  table.dataTable {
    border: var(--thin-border);
  }

  table.dataTable tbody tr:nth-child(odd) {
    background-color: #F8FCF9 !important;
  }

  table.dataTable tbody tr:nth-child(even) {
    background-color: #F7F1FC !important;
  }

  table.dataTable tbody tr:hover {
    background-color: #FFF1C9 !important;
  }

  pre,
  .shiny-text-output,
  code {
    border-radius: 6px;
    font-family: SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace;
  }

  .card-title-download {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 8px;
    width: 100%;
  }

  .header-download {
    min-width: 30px;
    min-height: 28px;
    padding: 2px 8px;
    line-height: 1;
  }

  .header-download .fa,
  .header-download svg {
    margin: 0;
  }
")
