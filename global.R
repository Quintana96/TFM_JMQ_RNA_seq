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

#' OrgDb instalados en el sistema, para el selector de organismo.
#'
#' La app tenia el OrgDb cableado a org.EcK12.eg.db, de modo que con datos de
#' cualquier otro organismo el enriquecimiento GO no podia funcionar aunque su
#' paquete de anotacion estuviera instalado: el mapeo salia del 0 % y la app
#' avisaba, correctamente, de algo que no tenia forma de arreglar desde la
#' interfaz. Se detectan al arrancar y se ofrecen todos.
available_orgdbs <- function() {
  pkgs <- tryCatch(rownames(utils::installed.packages()), error = function(e) character(0))
  sort(grep("^org[.].*[.]db$", pkgs, value = TRUE))
}
ORGDBS_DISPONIBLES <- available_orgdbs()

#' Nombre legible del organismo de un OrgDb, para no obligar a leer "org.Hs.eg.db".
orgdb_label <- function(pkg) {
  conocidos <- c(org.Hs.eg.db = "Homo sapiens", org.Mm.eg.db = "Mus musculus",
                 org.Rn.eg.db = "Rattus norvegicus", org.Dm.eg.db = "Drosophila melanogaster",
                 org.Ce.eg.db = "Caenorhabditis elegans", org.Sc.sgd.db = "Saccharomyces cerevisiae",
                 org.Dr.eg.db = "Danio rerio", org.At.tair.db = "Arabidopsis thaliana",
                 org.EcK12.eg.db = "Escherichia coli K12",
                 org.EcSakai.eg.db = "Escherichia coli Sakai")
  etiqueta <- unname(conocidos[pkg])
  ifelse(is.na(etiqueta), pkg, paste0(etiqueta, " (", pkg, ")"))
}
# Encogido de log2FC (lfcShrink). apeglm es la opcion recomendada por la
# vinieta de DESeq2; ashr sirve de alternativa. Sin ninguno de los dos se cae a
# type = "normal".
HAS_APEGLM          <- pkg_ok("apeglm")
HAS_ASHR            <- pkg_ok("ashr")
# Fase 2: GSEA, ponderacion de hipotesis y segunda estimacion de pi0. La interfaz
# oculta las opciones cuyo paquete falta, en lugar de ofrecerlas y fallar.
HAS_FGSEA           <- pkg_ok("fgsea")
# Reactome como tercera coleccion de enriquecimiento. Necesita los DOS paquetes:
# ReactomePA corre el test y reactome.db aporta los conjuntos (y es el que
# permite dibujar el running score sin salir a la red). Sin alguno de ellos, la
# coleccion no se ofrece en el selector en lugar de ofrecerse y fallar.
HAS_REACTOMEPA      <- pkg_ok("ReactomePA") && pkg_ok("reactome.db")
HAS_IHW             <- pkg_ok("IHW")
HAS_QVALUE          <- pkg_ok("qvalue")
# Fase 3: variacion no deseada y lectura de la anotacion. rtracklayer es lo que
# permite construir el mapa transcrito-gen que necesita tximport.
HAS_SVA             <- pkg_ok("sva")
HAS_RTRACKLAYER     <- pkg_ok("rtracklayer")
# Fase 4: motores robustos, Swish y calculo de potencia.
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


#' @section Paleta accesible
#' Colores de la codificacion "significativo / no significativo" en volcano y MA.
#'
#' El par anterior (verde #7BBF9A sobre gris #C0C0C0) tenia luminosidades muy
#' parecidas, de modo que con protanopia o deuteranopia —en torno al 8 % de los
#' hombres— los dos grupos resultaban casi indistinguibles. Se usa el par
#' azul/naranja de la paleta de Okabe-Ito, disenada para ser distinguible con
#' las formas habituales de daltonismo y que ademas conserva contraste al
#' imprimir en escala de grises.
DEG_SIG_COLORS <- c("Significativo" = "#0072B2", "No significativo" = "#E69F00")

#' Estados del checklist de configuracion. El color NO viaja solo: cada estado
#' lleva tambien un simbolo y un texto, porque el color por si mismo no es un
#' canal accesible ni lo lee un lector de pantalla.
#'
#' `nivel` es el modificador de la pildora (`.pill-ok`, `.pill-warn`,
#' `.pill-bad`), en lugar de un color literal: asi el checklist, los estados de
#' la portada y los badges de ejecucion comparten exactamente la misma paleta.
CHECKLIST_STATES <- list(
  ok       = list(nivel = "ok",   simbolo = "OK",       etiqueta = "completo"),
  optional = list(nivel = "warn", simbolo = "OPCIONAL", etiqueta = "opcional"),
  missing  = list(nivel = "bad",  simbolo = "FALTA",    etiqueta = "pendiente")
)


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
#'
#' El tema anterior daba a `primary`, `success`, `info`, `warning` y `danger`
#' valores pastel muy claros (#BEE8C8 y similares). Bootstrap deriva de esos
#' colores el fondo de los botones y el de las alertas, de modo que un boton
#' primario y uno de peligro quedaban casi igual de claros y la accion principal
#' no destacaba sobre las secundarias. Aqui los colores de marca son los tonos
#' saturados —que si contrastan como fondo de boton— y Bootstrap genera por su
#' cuenta las versiones claras que usan las alertas, en lugar de fijarlas a mano.
app_theme <- bs_theme(
  version   = 5,
  primary   = "#2E7D5B",
  secondary = "#5C6E64",
  success   = "#2E7D5B",
  info      = "#2A6F87",
  warning   = "#8A6D1C",
  danger    = "#A33A3A",
  base_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Arial", "sans-serif"),
  heading_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Arial", "sans-serif"),
  bg = "#F7FBF8",
  fg = "#1E2F27",
  "border-radius"     = "10px",
  "card-border-color" = "#D7E3DA",
  "card-cap-bg"       = "#F1F8F3",
  "body-font-size"    = "0.95rem"
)

#' @section Hoja de estilos
#'
#' Regla de la que se ha partido: el color solo aparece cuando significa algo.
#'
#' La version anterior pintaba las cabeceras de las tarjetas con
#' `.card:nth-of-type(3n + 1)`, `(3n + 2)` y `(3n)`, es decir, azul, naranja o
#' lavanda segun la POSICION de la tarjeta en el DOM. El color no decia nada del
#' contenido, cambiaba al anadir o quitar una tarjeta, y hacia que dos vistas con
#' las mismas tarjetas en distinto orden se vieran distintas. Lo mismo ocurria
#' con `.metric-card:nth-child(n)`. Ahora la cabecera es una sola y el color de
#' las metricas lo fija un modificador semantico (`is-ok`, `is-warn`, `is-bad`).
app_css <- HTML("
  :root {
    --bg:           #F7FBF8;
    --surface:      #FFFFFF;
    --surface-2:    #F4F9F5;
    --surface-3:    #EDF5EF;
    --ink:          #1E2F27;
    --ink-2:        #3C4F45;
    --muted:        #6B7C73;
    --line:         #D7E3DA;
    --line-strong:  #B6C9BC;
    --accent:       #2E7D5B;
    --accent-soft:  #E7F4EC;
    --info:         #2A6F87;
    --info-soft:    #E4F1F6;
    --warn:         #8A6D1C;
    --warn-soft:    #FBF1D8;
    --bad:          #A33A3A;
    --bad-soft:     #FAE5E4;
    --radius:       10px;
    --radius-sm:    8px;
    --shadow-sm:    0 1px 2px rgba(30, 47, 39, .06);
    --shadow:       0 1px 2px rgba(30, 47, 39, .05), 0 6px 18px rgba(30, 47, 39, .06);
  }

  body,
  .bslib-page-navbar {
    background: var(--bg);
    color: var(--ink);
  }

  /* ── Navbar ──────────────────────────────────────────────────────────────
     Antes era un degradado de tres colores (azul, lavanda, melocoton) que
     competia con el contenido y no se relacionaba con ninguna de las paletas
     de dentro. Ahora es una barra sobria y el unico acento es el subrayado de
     la pestana activa, que es la informacion que hay que ver. */
  .navbar {
    background: var(--surface) !important;
    border-bottom: 1px solid var(--line);
    box-shadow: var(--shadow-sm);
    padding-top: .35rem;
    padding-bottom: .35rem;
  }

  .navbar-brand {
    color: var(--ink) !important;
    font-weight: 700;
    letter-spacing: -.01em;
  }

  .navbar .nav-link {
    color: var(--ink-2) !important;
    font-weight: 500;
    border-radius: var(--radius-sm);
    padding: .35rem .7rem !important;
    border-bottom: 2px solid transparent;
  }

  .navbar .nav-link:hover {
    background: var(--surface-2);
    color: var(--ink) !important;
  }

  .navbar .nav-link.active {
    background: var(--accent-soft);
    color: var(--accent) !important;
    font-weight: 650;
    border-bottom-color: var(--accent);
  }

  /* ── Tarjetas ────────────────────────────────────────────────────────── */
  .card {
    border: 1px solid var(--line);
    border-radius: var(--radius);
    background: var(--surface);
    box-shadow: var(--shadow-sm);
  }

  .card-header {
    background: var(--surface-2);
    border-bottom: 1px solid var(--line);
    color: var(--ink);
    font-weight: 650;
    font-size: .92rem;
    letter-spacing: -.005em;
    padding: .6rem .9rem;
  }

  /* Titulo de seccion dentro de una pestana, para no encadenar cards sueltas
     sin jerarquia. */
  .section-title {
    font-size: 1.02rem;
    font-weight: 700;
    color: var(--ink);
    margin: 1.4rem 0 .1rem;
    letter-spacing: -.01em;
  }

  .section-subtitle {
    font-size: .85rem;
    color: var(--muted);
    margin-bottom: .7rem;
  }

  /* ── Botones ─────────────────────────────────────────────────────────────
     El selector anterior era button[id$='_btn'], que por especificidad ganaba
     a .btn-success y .btn-danger: como todos los botones de accion terminan en
     _btn, el de ejecutar (verde) y el de detener (rojo) se pintaban con el azul
     del selector de ficheros. Los selectores de fichero se marcan con
     .btn-picker y no se toca ningun boton por su id. */
  .btn {
    border-radius: var(--radius-sm);
    font-weight: 550;
  }

  .btn-lg {
    font-size: 1rem;
    padding: .55rem 1.1rem;
  }

  .btn-picker,
  .btn-file,
  .input-group .btn {
    background: var(--surface);
    border: 1px solid var(--line-strong);
    color: var(--info);
    font-weight: 550;
  }

  .btn-picker:hover,
  .btn-file:hover,
  .input-group .btn:hover {
    background: var(--info-soft);
    border-color: var(--info);
    color: var(--info);
  }

  .btn:disabled,
  .btn.disabled {
    opacity: .45;
  }

  /* ── Formularios ─────────────────────────────────────────────────────── */
  .form-control,
  .form-select,
  .selectize-input {
    border: 1px solid var(--line-strong) !important;
    border-radius: var(--radius-sm);
  }

  .form-control:focus,
  .form-select:focus,
  .selectize-input.focus {
    border-color: var(--accent) !important;
    box-shadow: 0 0 0 3px rgba(46, 125, 91, .14) !important;
  }

  .form-label,
  .shiny-input-container > label,
  .control-label {
    font-weight: 600;
    font-size: .86rem;
    color: var(--ink-2);
    margin-bottom: .25rem;
  }

  .form-check-label {
    font-weight: 400;
    font-size: .9rem;
  }

  /* Los textos de ayuda son abundantes en esta app y competian con el propio
     control por ser del mismo tamano. */
  .text-muted,
  small.text-muted {
    color: var(--muted) !important;
  }

  .card .form-group,
  .card .shiny-input-container {
    margin-bottom: .75rem;
  }

  /* ── Sub-pestanas ────────────────────────────────────────────────────── */
  .nav-tabs {
    border-bottom: 1px solid var(--line);
    gap: 2px;
  }

  .nav-tabs .nav-link {
    color: var(--ink-2);
    font-weight: 500;
    border: 1px solid transparent;
    border-radius: var(--radius-sm) var(--radius-sm) 0 0;
  }

  .nav-tabs .nav-link:hover {
    background: var(--surface-2);
    border-color: transparent;
  }

  .nav-tabs .nav-link.active {
    background: var(--surface);
    border-color: var(--line) var(--line) var(--surface);
    color: var(--accent);
    font-weight: 650;
  }

  .nav-pills .nav-link {
    color: var(--ink-2);
    font-weight: 500;
    border-radius: 999px;
    padding: .3rem .8rem;
    font-size: .88rem;
  }

  .nav-pills .nav-link.active {
    background: var(--accent-soft);
    color: var(--accent);
    font-weight: 650;
  }

  /* ── Acordeon ────────────────────────────────────────────────────────── */
  .accordion-button {
    font-weight: 650;
    font-size: .92rem;
    background: var(--surface-2);
    color: var(--ink);
    padding: .6rem .9rem;
  }

  .accordion-button:not(.collapsed) {
    background: var(--accent-soft);
    color: var(--accent);
    box-shadow: none;
  }

  .accordion-button:focus {
    box-shadow: 0 0 0 3px rgba(46, 125, 91, .14);
  }

  /* ── Metricas ────────────────────────────────────────────────────────────
     Una sola forma para todas: el color lo pone el estado (is-ok / is-warn /
     is-bad), no la posicion en la fila. */
  .stat-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 10px;
  }

  .metric-card {
    padding: .7rem .85rem;
    border-radius: var(--radius);
    border: 1px solid var(--line);
    border-left: 3px solid var(--line-strong);
    background: var(--surface);
    box-shadow: var(--shadow-sm);
    min-height: 76px;
    display: flex;
    flex-direction: column;
    justify-content: center;
  }

  .metric-card.is-ok   { border-left-color: var(--accent); }
  .metric-card.is-warn { border-left-color: var(--warn); }
  .metric-card.is-bad  { border-left-color: var(--bad); }
  .metric-card.is-info { border-left-color: var(--info); }

  .metric-card-value {
    font-size: 1.28rem;
    font-weight: 700;
    color: var(--ink);
    line-height: 1.2;
    font-variant-numeric: tabular-nums;
  }

  .metric-card-label {
    font-size: .7rem;
    color: var(--muted);
    margin-top: .35rem;
    text-transform: uppercase;
    letter-spacing: .05em;
    font-weight: 600;
  }

  /* ── Portada ─────────────────────────────────────────────────────────── */
  .home-hero {
    border: 1px solid var(--line);
    border-radius: var(--radius);
    background: linear-gradient(120deg, #EDF7F0 0%, #F4FAF6 55%, #FFFFFF 100%);
    padding: 1.4rem 1.6rem;
    margin-bottom: 1.1rem;
    box-shadow: var(--shadow-sm);
  }

  .home-hero h1 {
    font-size: 1.5rem;
    font-weight: 700;
    margin: 0 0 .3rem;
    letter-spacing: -.02em;
  }

  .home-hero p {
    color: var(--ink-2);
    margin: 0;
    max-width: 68ch;
    font-size: .93rem;
  }

  .home-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(255px, 1fr));
    gap: 14px;
    align-items: stretch;
  }

  /* La tarjeta entera es el area pulsable: obligar a acertar un boton pequeno
     al final de la tarjeta es justo lo que hacia torpe la navegacion. */
  .home-card {
    position: relative;
    display: flex;
    flex-direction: column;
    gap: .5rem;
    padding: 1.05rem 1.1rem 1.1rem;
    border: 1px solid var(--line);
    border-radius: var(--radius);
    background: var(--surface);
    box-shadow: var(--shadow-sm);
    cursor: pointer;
    transition: border-color .12s ease, box-shadow .12s ease, transform .12s ease;
    text-align: left;
    width: 100%;
  }

  .home-card:hover,
  .home-card:focus-visible {
    border-color: var(--accent);
    box-shadow: var(--shadow);
    transform: translateY(-2px);
    outline: none;
  }

  .home-card.is-locked {
    cursor: not-allowed;
    background: var(--surface-2);
    box-shadow: none;
  }

  .home-card.is-locked:hover {
    border-color: var(--line);
    transform: none;
    box-shadow: none;
  }

  .home-card-top {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: .5rem;
  }

  .home-card-step {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 26px;
    height: 26px;
    border-radius: 50%;
    background: var(--accent-soft);
    color: var(--accent);
    font-weight: 700;
    font-size: .8rem;
    flex: 0 0 auto;
  }

  .home-card.is-locked .home-card-step {
    background: var(--surface-3);
    color: var(--muted);
  }

  .home-card-title {
    font-size: 1rem;
    font-weight: 650;
    color: var(--ink);
    margin: 0;
  }

  .home-card-desc {
    font-size: .86rem;
    color: var(--muted);
    margin: 0;
    flex: 1 1 auto;
  }

  .home-card-cta {
    font-size: .84rem;
    font-weight: 650;
    color: var(--accent);
  }

  .home-card.is-locked .home-card-cta {
    color: var(--muted);
  }

  /* ── Pildoras de estado ──────────────────────────────────────────────────
     El estado viaja siempre en texto ademas de en color: un circulo de color no
     lo distingue quien tiene daltonismo ni lo anuncia un lector de pantalla. */
  .pill {
    display: inline-flex;
    align-items: center;
    gap: .3rem;
    font-size: .7rem;
    font-weight: 700;
    letter-spacing: .03em;
    padding: .15rem .5rem;
    border-radius: 999px;
    white-space: nowrap;
    border: 1px solid transparent;
  }

  .pill-ok      { background: var(--accent-soft); color: var(--accent);  border-color: #BFE0CC; }
  .pill-warn    { background: var(--warn-soft);   color: var(--warn);    border-color: #E8D6A4; }
  .pill-bad     { background: var(--bad-soft);    color: var(--bad);     border-color: #EEC3C1; }
  .pill-info    { background: var(--info-soft);   color: var(--info);    border-color: #BFD9E3; }
  .pill-neutral { background: var(--surface-3);   color: var(--muted);   border-color: var(--line); }

  /* ── Checklist ───────────────────────────────────────────────────────── */
  .checklist {
    list-style: none;
    margin: 0;
    padding: 0;
  }

  .checklist li {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: .6rem;
    padding: .4rem 0;
    font-size: .9rem;
    border-bottom: 1px dashed var(--line);
  }

  .checklist li:last-child {
    border-bottom: 0;
  }

  /* ── Tablas ──────────────────────────────────────────────────────────────
     La zebra anterior alternaba mint y lavanda y el hover era amarillo, tres
     tonos compitiendo con el contenido de la celda. */
  table.dataTable thead th {
    background: var(--surface-2);
    color: var(--ink);
    font-weight: 650;
    border-bottom: 1px solid var(--line-strong) !important;
  }

  table.dataTable {
    border: 1px solid var(--line);
    font-size: .88rem;
  }

  table.dataTable tbody tr:nth-child(odd)  { background-color: var(--surface) !important; }
  table.dataTable tbody tr:nth-child(even) { background-color: #FBFDFB !important; }
  table.dataTable tbody tr:hover           { background-color: var(--accent-soft) !important; }

  .dataTables_wrapper .dataTables_filter input {
    border: 1px solid var(--line-strong);
    border-radius: var(--radius-sm);
    padding: .2rem .5rem;
  }

  /* ── Log y codigo ────────────────────────────────────────────────────── */
  pre,
  .shiny-text-output,
  code {
    border-radius: var(--radius-sm);
    font-family: SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace;
  }

  pre {
    font-size: .82rem;
  }

  /* Las rutas absolutas del proyecto son largas y en pantallas estrechas
     desbordaban su columna. */
  code {
    overflow-wrap: anywhere;
    font-size: .82rem;
  }

  .log-pre pre {
    background: #14201A;
    color: #DCEFE3;
    border: 0;
  }

  /* ── Cabecera de tarjeta con controles a la derecha ──────────────────── */
  .card-title-download {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 10px;
    width: 100%;
    flex-wrap: wrap;
  }

  .card-header-tools {
    display: flex;
    gap: 6px;
    align-items: center;
    flex-wrap: wrap;
  }

  /* Los selectores que viven en la cabecera de una tarjeta arrastraban el
     margen inferior del form-group y descuadraban verticalmente el titulo. */
  .card-header-tools .form-group,
  .card-header-tools .shiny-input-container {
    margin-bottom: 0 !important;
  }

  .card-header-tools label {
    display: none;
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

  /* ── Alertas ─────────────────────────────────────────────────────────── */
  .alert {
    border-radius: var(--radius-sm);
    border-width: 1px;
  }

  .alert-secondary {
    background: var(--surface-2);
    border-color: var(--line);
    color: var(--ink-2);
  }

  /* ── Barra lateral de parametros ─────────────────────────────────────── */
  .bslib-sidebar-layout > .sidebar {
    background: var(--surface-2);
  }

  .bslib-sidebar-layout > .sidebar .accordion {
    --bs-accordion-bg: var(--surface);
  }
")
