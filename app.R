# ============================================================
# app.R  —  RNA-seq Workflow Runner  v3.0
# ============================================================
# v2 → v3:
#   [v3-PROC]   processx: ejecucion en tiempo real (no bloqueante)
#   [v3-LAYOUT] Tab 2: progreso arriba + log abajo ancho completo
#   [v3-PROG]   Tiempo transcurrido, ETA, estado por muestra
#   [v3-LOAD]   Carga desde resultados previos (conteos / directorio)
#   [v3-DEG]    Tab 4: DESeq2 — PCA, Volcano, Heatmap, Tabla
#   [v3-FUNC]   Tab 5: ORA y GSEA con clusterProfiler
#   [v3-CONV]   Validacion de nombres de muestras (buenas practicas)
# ============================================================

library(shiny)
library(bslib)
library(shinyjs)
library(DT)

# ── Paquetes opcionales — deteccion en tiempo de carga ───────────────────────
pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)
HAS_PROCESSX   <- pkg_ok("processx")
HAS_DESEQ2     <- pkg_ok("DESeq2")
HAS_TXIMPORT   <- pkg_ok("tximport")
HAS_GGPLOT2    <- pkg_ok("ggplot2")
HAS_PLOTLY     <- pkg_ok("plotly")
HAS_PHEATMAP   <- pkg_ok("pheatmap")
HAS_CLUSTERP   <- pkg_ok("clusterProfiler")
HAS_ENRICHPLOT <- pkg_ok("enrichplot")
HAS_GGREPEL    <- pkg_ok("ggrepel")

if (HAS_GGPLOT2)  suppressPackageStartupMessages(library(ggplot2))
if (HAS_PLOTLY)   suppressPackageStartupMessages(library(plotly))
if (HAS_PHEATMAP) suppressPackageStartupMessages(library(pheatmap))


# ── Operador nulo-coalescente ─────────────────────────────────────────────────
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  HELPERS                                                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── [v2-PRESERVED] Deteccion de muestras FASTQ ───────────────────────────────
detect_samples <- function(dir_path) {
  if (!nzchar(dir_path) || !dir.exists(dir_path)) return(character(0))
  files <- list.files(
    dir_path,
    pattern = "(_1\\.fastq\\.gz$|_R1\\.fastq\\.gz$|_1\\.fastq$|_R1\\.fastq$)",
    full.names = FALSE, ignore.case = TRUE
  )
  gsub("(_1\\.fastq\\.gz|_R1\\.fastq\\.gz|_1\\.fastq|_R1\\.fastq)$",
       "", files, ignore.case = TRUE)
}

missing_r2 <- function(dir_path, samples) {
  if (length(samples) == 0) return(logical(0))
  vapply(samples, function(s) {
    candidates <- file.path(
      dir_path,
      paste0(s, c("_2.fastq.gz", "_R2.fastq.gz", "_2.fastq", "_R2.fastq"))
    )
    !any(file.exists(candidates))
  }, logical(1))
}

fmt_bytes <- function(b) {
  b <- as.numeric(b)
  if (is.na(b))    return("—")
  if (b >= 1e9)    return(sprintf("%.1f GB", b / 1e9))
  if (b >= 1e6)    return(sprintf("%.1f MB", b / 1e6))
  if (b >= 1e3)    return(sprintf("%.1f KB", b / 1e3))
  paste0(b, " B")
}

ts_log <- function(msg) paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

# ── [v3-NEW] Helpers adicionales ─────────────────────────────────────────────

#' Formatea segundos en cadena legible (3s, 2m 15s, 1h 4m)
fmt_elapsed <- function(secs) {
  s <- round(as.numeric(secs))
  if (s < 60)   return(sprintf("%ds", s))
  if (s < 3600) return(sprintf("%dm %ds", s %/% 60, s %% 60))
  sprintf("%dh %dm", s %/% 3600, (s %% 3600) %/% 60)
}

#' Valida que el nombre de muestra no tenga caracteres problematicos
bad_sample_chars <- function(names) {
  grep("[^A-Za-z0-9_.-]", names, value = TRUE)
}

#' Carga la matriz de conteos generada por el workflow
load_counts_from_workflow <- function(output_dir, tool) {
  if (tool == "bowtie2") {
    f <- file.path(output_dir, "04_counts", "count_matrix.tsv")
    if (!file.exists(f)) return(NULL)
    tryCatch(
      read.table(f, header = TRUE, row.names = 1, sep = "\t", comment.char = "#"),
      error = function(e) NULL
    )
  } else if (HAS_TXIMPORT) {
    aln_dir <- file.path(output_dir, "03_alignments", tool)
    if (!dir.exists(aln_dir)) return(NULL)
    sdirs <- list.dirs(aln_dir, recursive = FALSE, full.names = TRUE)
    if (length(sdirs) == 0) return(NULL)
    qfiles <- if (tool == "salmon")
      file.path(sdirs, "quant.sf")
    else
      file.path(sdirs, "abundance.h5")
    if (!all(file.exists(qfiles)))
      qfiles <- file.path(sdirs, "abundance.tsv")
    valid <- file.exists(qfiles)
    if (!any(valid)) return(NULL)
    qfiles <- qfiles[valid]
    names(qfiles) <- basename(sdirs[valid])
    txi <- tryCatch(
      tximport::tximport(qfiles, type = tool, ignoreTxVersion = TRUE),
      error = function(e) NULL
    )
    if (is.null(txi)) return(NULL)
    round(txi$counts)
  } else NULL
}

#' Informe de reproducibilidad
build_repro_text <- function(p, workflow_path) {
  if (length(p) == 0) return("")
  pkgs <- c("shiny", "bslib", "shinyjs", "DT",
            "processx", "DESeq2", "tximport", "clusterProfiler")
  pkg_ver <- sapply(pkgs, function(x)
    tryCatch(as.character(packageVersion(x)), error = function(e) "—"))
  paste(c(
    "=== Informe de reproducibilidad ===",
    paste0("Fecha/hora          : ", format(p$started_at, "%Y-%m-%d %H:%M:%S")),
    paste0("Tipo de analisis    : ", p$analysis_type),
    paste0("Herramienta         : ", p$tool),
    paste0("Num. muestras       : ", p$n_samples),
    paste0("Directorio entrada  : ", p$input_dir),
    paste0("Directorio salida   : ", p$output_dir),
    paste0("Referencia (FASTA)  : ", p$genome_file),
    paste0("Anotacion (GFF/GTF) : ", p$annotation_file),
    paste0("Version R           : ", p$r_version),
    paste0("Paquetes R          : ",
           paste(mapply(paste, pkgs, pkg_ver, sep = " "), collapse = ", ")),
    paste0("workflow.sh         : ", workflow_path)
  ), collapse = "\n")
}

#' Crea un dotplot ggplot2 para resultados ORA/GSEA
make_enrichment_dotplot <- function(result_df, title = "") {
  df <- result_df
  df <- head(df[order(df$p.adjust), ], 20)
  if (nrow(df) == 0) return(NULL)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  df$GR <- vapply(df$GeneRatio, function(x) {
    p <- strsplit(x, "/")[[1]]
    as.numeric(p[1]) / as.numeric(p[2])
  }, numeric(1))
  ggplot(df, aes(x = GR, y = Description,
                 color = p.adjust, size = Count,
                 text = paste0(Description, "\np.adjust: ",
                               formatC(p.adjust, digits = 3, format = "e"),
                               "\nGenes: ", Count))) +
    geom_point() +
    scale_color_gradient(low = "#B2182B", high = "#2166AC",
                         name = "p.adjust",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(name = "N genes") +
    xlab("Gene Ratio") + ylab(NULL) +
    ggtitle(title) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
}


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  UI                                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

ui <- page_navbar(
  title  = tags$span(icon("dna"), " RNA-seq Workflow Runner"),
  id     = "main_nav",
  theme  = bs_theme(bootswatch = "flatly", primary = "#2C7BB6"),
  header = useShinyjs(),

  # ════════════════════════════════════════════════════════════════════════
  # TAB 1  —  CONFIGURACION
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = tagList(icon("sliders"), " 1 \u00b7 Configuracion"),
    value = "tab_config",

    # ── Modo de inicio ────────────────────────────────────────────────────
    card(
      card_header(icon("rocket"), " Modo de inicio"),
      radioButtons(
        "start_mode", label = NULL,
        choices  = c("Ejecutar workflow completo" = "workflow",
                     "Cargar resultados previos"  = "load"),
        selected = "workflow", inline = TRUE
      ),
      # ── Cargar resultados previos ─────────────────────────────────────
      conditionalPanel(
        condition = "input.start_mode === 'load'",
        layout_columns(
          col_widths = c(6, 6),
          div(
            tags$strong("Matriz de conteos (TSV/CSV)"),
            tags$br(),
            tags$small(class="text-muted",
              "Genes como filas, muestras como columnas. Primera columna = ID de gen."),
            fileInput("upload_counts", label = NULL,
                      accept = c(".tsv", ".csv", ".txt")),
            tags$strong("O: directorio de resultados previos"),
            tags$br(),
            textInput("prev_output_dir", label = NULL,
                      placeholder = "/data/results_anteriores")
          ),
          div(
            tags$strong("Metadatos (CSV: sample,condition[,batch])"),
            fileInput("upload_meta", label = NULL, accept = ".csv"),
            tags$small(class = "text-muted",
              icon("circle-info"),
              " Columna 'sample' debe coincidir con nombres de columna de la matriz.")
          )
        ),
        div(
          style = "text-align:right; margin-top:8px;",
          actionButton("btn_load_existing", tagList(icon("upload"), " Cargar y continuar"),
                       class = "btn-success btn-lg")
        )
      )
    ),

    # ── Configuracion del workflow (solo si start_mode == workflow) ────────
    conditionalPanel(
      condition = "input.start_mode === 'workflow'",

      layout_columns(
        col_widths = c(4, 8),

        # Panel izquierdo: tipo de analisis
        card(
          card_header(icon("flask"), " Tipo de analisis"),
          radioButtons(
            "analysis_type", label = NULL,
            choices  = c("Alineamiento"      = "alignment",
                         "Pseudoalineamiento" = "pseudo"),
            selected = "alignment"
          ),
          conditionalPanel(
            "input.analysis_type === 'alignment'",
            tags$p(class = "text-muted small mb-0",
                   icon("circle-info"), " Bowtie2 + featureCounts")
          ),
          conditionalPanel(
            "input.analysis_type === 'pseudo'",
            selectInput("pseudo_tool", "Herramienta",
                        choices  = c("Salmon" = "salmon", "Kallisto" = "kallisto"),
                        selected = "salmon")
          ),
          hr(),
          tags$strong("Tipo de lectura"),
          radioButtons(
            "read_type", label = NULL,
            choices  = c("Paired-end" = "pe", "Single-end" = "se"),
            selected = "pe"
          ),
          conditionalPanel(
            "input.read_type === 'se'",
            div(class = "alert alert-warning py-2 mt-1",
                icon("triangle-exclamation"),
                tags$small(" workflow.sh solo soporta paired-end."))
          )
        ),

        # Panel derecho: rutas
        card(
          card_header(icon("folder-open"), " Rutas y archivos"),
          textInput("input_dir", "Directorio de FASTQs de entrada",
                    placeholder = "/data/reads"),
          uiOutput("genome_label_ui"),
          conditionalPanel(
            "input.analysis_type === 'alignment'",
            textInput("annotation_file",
                      "Archivo de anotacion GFF/GTF (requerido)",
                      placeholder = "/refs/ecoli.gff")
          ),
          conditionalPanel(
            "input.analysis_type === 'pseudo'",
            textInput("annotation_file_pseudo",
                      "Archivo de anotacion GFF/GTF (opcional)",
                      placeholder = "/refs/ecoli.gff"),
            tags$small(class = "text-muted", icon("circle-info"),
                       " Si se omite, se pasa /dev/null al script.")
          ),
          textInput("output_dir", "Directorio de salida",
                    placeholder = "/data/results"),
          hr(),
          tags$small(class = "text-muted", icon("lightbulb"),
                     " Buenas practicas: nombres de muestra sin espacios ni",
                     " caracteres especiales (usar solo letras, numeros, _ - .)")
        )
      ),

      # Deteccion de muestras
      card(
        card_header(icon("microscope"), " Muestras detectadas"),
        uiOutput("sample_preview_ui")
      ),

      # Validacion y navegacion
      card(
        uiOutput("validation_ui"),
        div(
          style = "text-align:right; margin-top:12px;",
          actionButton("btn_to_processing",
                       tagList("Continuar al procesamiento", icon("arrow-right")),
                       class = "btn-primary btn-lg")
        )
      )
    ) # end conditionalPanel workflow
  ), # end tab_config


  # ════════════════════════════════════════════════════════════════════════
  # TAB 2  —  PROCESAMIENTO  (layout rediseñado: progreso arriba, log abajo)
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = tagList(icon("gears"), " 2 \u00b7 Procesamiento"),
    value = "tab_process",
    uiOutput("tab2_content")
  ),


  # ════════════════════════════════════════════════════════════════════════
  # TAB 3  —  RESULTADOS WORKFLOW
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = tagList(icon("folder-open"), " 3 \u00b7 Resultados"),
    value = "tab_results",
    uiOutput("tab3_content")
  ),


  # ════════════════════════════════════════════════════════════════════════
  # TAB 4  —  ANALISIS DE EXPRESION DIFERENCIAL (DEG)
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = tagList(icon("chart-line"), " 4 \u00b7 DEG"),
    value = "tab_deg",
    uiOutput("tab4_gate_ui"),
    conditionalPanel(
      condition = "output.counts_ready_flag === true",
      navset_tab(
        id = "deg_subtabs",
        # ── Configuracion DEG ───────────────────────────────────────────
        nav_panel("Configuracion",
          layout_columns(
            col_widths = c(5, 7),
            card(
              card_header(icon("table"), " Matriz de conteos"),
              uiOutput("counts_summary_ui"),
              hr(),
              tags$strong("Condiciones (separadas por coma)"),
              textInput("condition_labels", label = NULL,
                        placeholder = "control, tratamiento"),
              uiOutput("condition_table_ui"),
              hr(),
              selectInput("deg_ref_cond", "Condicion de referencia (control)",
                          choices = NULL),
              selectInput("deg_treat_cond", "Condicion de tratamiento",
                          choices = NULL)
            ),
            card(
              card_header(icon("sliders"), " Parametros DESeq2"),
              sliderInput("deg_min_counts", "Suma minima de conteos por gen",
                          1, 100, 10, step = 1),
              sliderInput("deg_fdr",   "FDR (alpha)",   0.001, 0.2,  0.05, step = 0.001),
              sliderInput("deg_lfc",   "LFC minimo (|log2FC|)", 0, 3, 1, step = 0.1),
              hr(),
              div(
                style = "text-align:right;",
                actionButton("run_deg", tagList(icon("play"), " Ejecutar DESeq2"),
                             class = "btn-success btn-lg")
              ),
              uiOutput("deg_status_ui")
            )
          )
        ),
        # ── PCA ─────────────────────────────────────────────────────────
        nav_panel("PCA",
          card(
            card_header(icon("circle-dot"), " Analisis de Componentes Principales (PCA)"),
            layout_columns(col_widths = c(3, 9),
              div(
                selectInput("pca_x", "Eje X", choices = c("PC1","PC2","PC3"), selected = "PC1"),
                selectInput("pca_y", "Eje Y", choices = c("PC1","PC2","PC3"), selected = "PC2"),
                checkboxInput("pca_labels", "Mostrar etiquetas", TRUE),
                downloadButton("dl_pca", "Descargar PNG", class = "btn-sm btn-outline-secondary w-100")
              ),
              plotlyOutput("pca_plot", height = "480px")
            )
          )
        ),
        # ── Volcano ──────────────────────────────────────────────────────
        nav_panel("Volcano Plot",
          card(
            card_header(icon("volcano"), " Volcano Plot"),
            layout_columns(col_widths = c(3, 9),
              div(
                sliderInput("vol_fdr",  "FDR cutoff",  0.001, 0.2,  0.05, step = 0.001),
                sliderInput("vol_lfc",  "|LFC| cutoff", 0,    3,    1,    step = 0.1),
                sliderInput("vol_top",  "Top N etiquetados", 0, 30, 10),
                downloadButton("dl_volcano", "Descargar PNG",
                               class = "btn-sm btn-outline-secondary w-100")
              ),
              plotlyOutput("volcano_plot", height = "480px")
            )
          )
        ),
        # ── Heatmap ──────────────────────────────────────────────────────
        nav_panel("Heatmap",
          card(
            card_header(icon("th"), " Heatmap de genes DEG"),
            layout_columns(col_widths = c(3, 9),
              div(
                sliderInput("hm_top_n",  "Top N genes (por padj)",     10, 200, 50),
                sliderInput("hm_fdr",    "FDR maximo de seleccion",  0.001, 0.2, 0.05, step=0.001),
                selectInput("hm_clust",  "Metodo de clustering",
                            choices = c("ward.D2","complete","average","single","none"),
                            selected = "ward.D2"),
                selectInput("hm_palette","Paleta de colores",
                            choices = c("RdBu","RdYlBu","PuOr","Spectral"),
                            selected = "RdBu"),
                selectInput("hm_scale",  "Escala",
                            choices = c("Por gen (z-score)" = "row",
                                        "Sin escala"         = "none"),
                            selected = "row"),
                checkboxInput("hm_rownames", "Mostrar nombres de genes", FALSE),
                downloadButton("dl_heatmap", "Descargar PNG",
                               class = "btn-sm btn-outline-secondary w-100")
              ),
              plotOutput("heatmap_plot",
                         height = "600px", width = "100%")
            )
          )
        ),
        # ── Tabla estadistica ────────────────────────────────────────────
        nav_panel("Tabla estadistica",
          card(
            card_header(icon("table"), " Resultados DESeq2"),
            div(class = "alert alert-info py-2 small",
                icon("circle-info"),
                tags$b(" Columnas: "),
                "log2FoldChange = cambio relativo entre condiciones; ",
                "padj = p-valor ajustado por FDR (Benjamini-Hochberg); ",
                "stat = estadistico de Wald; baseMean = expresion media."),
            layout_columns(col_widths = c(9, 3),
              DTOutput("deg_table"),
              div(class = "d-flex flex-column gap-2",
                  downloadButton("dl_deg_csv", tagList(icon("download"), " CSV completo"),
                                 class = "btn-outline-secondary w-100"),
                  downloadButton("dl_deg_sig", tagList(icon("download"), " Solo significativos"),
                                 class = "btn-outline-secondary w-100"))
            )
          )
        )
      ) # end navset_tab deg
    ) # end conditionalPanel counts_ready
  ), # end tab_deg


  # ════════════════════════════════════════════════════════════════════════
  # TAB 5  —  ANALISIS FUNCIONAL (ORA / GSEA)
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = tagList(icon("magnifying-glass-chart"), " 5 \u00b7 Funcional"),
    value = "tab_func",
    uiOutput("tab5_gate_ui"),
    conditionalPanel(
      condition = "output.deg_done_flag === true",
      # Configuracion comun de organismo
      card(
        card_header(icon("bacterium"), " Organismo y base de datos"),
        layout_columns(col_widths = c(4, 4, 4),
          div(
            tags$strong("Codigo KEGG del organismo"),
            tags$br(),
            tags$small(class="text-muted", "Ej: eco (E.coli K-12), hsa (humano), mmu (raton)"),
            textInput("kegg_org", label = NULL, value = "eco", placeholder = "eco")
          ),
          selectInput("func_db", "Base de datos",
                      choices = c("KEGG" = "kegg", "GO BP" = "gobp",
                                  "GO MF" = "gomf", "GO CC" = "gocc"),
                      selected = "kegg"),
          div(
            tags$strong("Paquete OrgDb (solo GO)"),
            tags$br(),
            tags$small(class="text-muted",
                       "Ej: org.EcK12.eg.db, org.Hs.eg.db, org.Mm.eg.db"),
            textInput("orgdb_pkg", label = NULL, placeholder = "org.EcK12.eg.db"),
            tags$small(class="text-muted",
                       "keyType: "),
            textInput("orgdb_keytype", label = NULL, value = "SYMBOL",
                      placeholder = "SYMBOL")
          )
        )
      ),
      navset_tab(
        id = "func_subtabs",
        # ── ORA ─────────────────────────────────────────────────────────
        nav_panel("ORA",
          card(
            card_header(icon("circle-nodes"), " Analisis de sobrerepresentacion (ORA)"),
            layout_columns(col_widths = c(3, 9),
              div(
                sliderInput("ora_fdr",  "FDR cutoff genes DEG",  0.001, 0.2,  0.05, step=0.001),
                sliderInput("ora_lfc",  "|LFC| cutoff genes DEG", 0,    3,    1,    step=0.1),
                sliderInput("ora_pval", "p.adjust cutoff ORA",    0.001, 0.2,  0.05, step=0.001),
                sliderInput("ora_qval", "q-value cutoff ORA",     0.01,  0.5,  0.2,  step=0.01),
                actionButton("run_ora", tagList(icon("play"), " Ejecutar ORA"),
                             class = "btn-success w-100"),
                br(), br(),
                downloadButton("dl_ora_csv", tagList(icon("download"), " CSV"),
                               class = "btn-sm btn-outline-secondary w-100")
              ),
              div(
                plotlyOutput("ora_dot_plot",  height = "380px"),
                DTOutput("ora_table")
              )
            )
          )
        ),
        # ── GSEA ─────────────────────────────────────────────────────────
        nav_panel("GSEA",
          card(
            card_header(icon("wave-square"), " Analisis de enriquecimiento de conjuntos genicos (GSEA)"),
            layout_columns(col_widths = c(3, 9),
              div(
                selectInput("gsea_rank", "Metrica de ordenacion",
                            choices = c("Estadistico Wald" = "stat",
                                        "log2FoldChange"   = "lfc"),
                            selected = "stat"),
                sliderInput("gsea_pval", "p.adjust cutoff GSEA", 0.001, 0.5, 0.05, step=0.001),
                numericInput("gsea_min", "Tamano minimo gene set", 10,  min=5,  max=100),
                numericInput("gsea_max", "Tamano maximo gene set", 500, min=50, max=2000),
                actionButton("run_gsea", tagList(icon("play"), " Ejecutar GSEA"),
                             class = "btn-success w-100"),
                br(), br(),
                downloadButton("dl_gsea_csv", tagList(icon("download"), " CSV"),
                               class = "btn-sm btn-outline-secondary w-100")
              ),
              div(
                plotOutput("gsea_bar_plot", height = "350px"),
                DTOutput("gsea_table")
              )
            )
          )
        )
      ) # end navset_tab func
    ) # end conditionalPanel deg_done
  ) # end tab_func

) # end page_navbar


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SERVER                                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

server <- function(input, output, session) {

  # ── Estado global ─────────────────────────────────────────────────────────
  workflow_path <- normalizePath("workflow.sh", mustWork = FALSE)

  # Log de ejecucion
  log_text <- reactiveVal(paste0(ts_log("Workflow listo.\n")))

  # Archivos en el directorio de salida
  output_files_rv <- reactiveVal(
    data.frame(Archivo = character(), `Tamano` = character(),
               stringsAsFactors = FALSE, check.names = FALSE)
  )

  # Navegacion entre pestanas
  process_unlocked <- reactiveVal(FALSE)
  analysis_done    <- reactiveVal(FALSE)
  config_snap      <- reactiveVal(list())
  run_params_rv    <- reactiveVal(list())

  # Estado del proceso (real-time con processx)          [v3-PROC]
  proc_rv <- reactiveValues(
    proc       = NULL,
    running    = FALSE,
    start_time = NULL,
    checkpoints = character(0),
    cp_idx     = 0,
    samp_stat  = list(),
    cur_sample = NULL,
    n_total    = 0
  )

  # Datos para DEG y analisis funcional                  [v3-DEG / v3-FUNC]
  data_rv <- reactiveValues(
    count_matrix = NULL,
    counts_ready = FALSE,
    source       = "workflow",
    deg_dds      = NULL,
    deg_res      = NULL,
    deg_done     = FALSE,
    ora_res      = NULL,
    gsea_res     = NULL
  )

  # Flags de estado para conditionalPanel (suspendWhenHidden=FALSE)
  output$counts_ready_flag <- reactive({ data_rv$counts_ready })
  outputOptions(output, "counts_ready_flag", suspendWhenHidden = FALSE)
  output$deg_done_flag <- reactive({ data_rv$deg_done })
  outputOptions(output, "deg_done_flag", suspendWhenHidden = FALSE)

  # ── Timer de polling (processx)                       [v3-PROC] ────────────
  poll_timer <- reactiveTimer(500)

  # ── Observer de polling ────────────────────────────────────────────────────
  observe({
    req(HAS_PROCESSX, proc_rv$running, !is.null(proc_rv$proc))
    poll_timer()

    proc <- proc_rv$proc
    if (!proc$is_alive()) {
      # Leer salida restante
      remaining <- tryCatch(
        c(proc$read_output_lines(), proc$read_error_lines()),
        error = function(e) character(0)
      )
      if (length(remaining) > 0)
        log_text(paste0(log_text(), paste(remaining, collapse = "\n"), "\n"))

      exit_code <- proc$get_exit_status() %||% 0L
      proc_rv$running <- FALSE
      log_text(paste0(log_text(), ts_log(paste0("Codigo de salida: ", exit_code)), "\n"))

      # Marcar muestra actual como terminada
      if (!is.null(proc_rv$cur_sample)) {
        ss <- proc_rv$samp_stat
        ss[[proc_rv$cur_sample]] <- "done"
        proc_rv$samp_stat <- ss
        proc_rv$cur_sample <- NULL
      }
      proc_rv$cp_idx <- length(proc_rv$checkpoints)

      if (exit_code == 0) {
        # Auto-cargar matriz de conteos
        p <- run_params_rv()
        counts <- tryCatch(load_counts_from_workflow(p$output_dir, p$tool), error=function(e) NULL)
        if (!is.null(counts)) {
          data_rv$count_matrix <- counts
          data_rv$counts_ready <- TRUE
          log_text(paste0(log_text(),
            ts_log(sprintf("Matriz cargada: %d genes x %d muestras",
                           nrow(counts), ncol(counts))), "\n"))
        }
        # Actualizar lista de archivos
        files <- list.files(p$output_dir, recursive = TRUE, full.names = FALSE)
        if (length(files) > 0) {
          sz <- file.info(file.path(p$output_dir, files))$size
          output_files_rv(data.frame(
            Archivo = files, `Tamano` = sapply(sz, fmt_bytes),
            stringsAsFactors = FALSE, check.names = FALSE
          ))
        }
        log_text(paste0(log_text(), ts_log("=== Analisis finalizado OK ==="), "\n"))
        analysis_done(TRUE)
        showNotification("Workflow finalizado correctamente.", type = "message")
        updateNavbarPage(session, "main_nav", selected = "tab_results")
      } else {
        log_text(paste0(log_text(), ts_log("=== ERROR en el workflow ==="), "\n"))
        showNotification(
          paste0("Error (codigo ", exit_code, "). Revisa el log de ejecucion."),
          type = "error", duration = 12
        )
      }
    } else {
      # Proceso vivo: leer nueva salida
      poll_res <- proc$poll_io(timeout = 0)
      new_lines <- character(0)
      if (identical(poll_res["output"], "ready"))
        new_lines <- c(new_lines, proc$read_output_lines())
      if (identical(poll_res["error"], "ready"))
        new_lines <- c(new_lines, proc$read_error_lines())

      if (length(new_lines) > 0) {
        log_text(paste0(log_text(), paste(new_lines, collapse = "\n"), "\n"))
        for (line in new_lines) {
          # Avance de checkpoints basado en patrones del workflow.sh
          if (grepl("Building|bowtie2-build|salmon index|kallisto index", line, ignore.case=TRUE))
            proc_rv$cp_idx <- max(proc_rv$cp_idx, 1L)
          if (grepl("fastqc", line, ignore.case=TRUE) && proc_rv$cp_idx < 2)
            proc_rv$cp_idx <- 2L
          if (grepl("^Processing sample:", line, ignore.case=TRUE)) {
            proc_rv$cp_idx <- max(proc_rv$cp_idx, 3L)
            sname <- trimws(sub("^Processing sample:", "", line, ignore.case=TRUE))
            cur   <- proc_rv$cur_sample
            if (!is.null(cur)) {
              ss <- proc_rv$samp_stat; ss[[cur]] <- "done"; proc_rv$samp_stat <- ss
            }
            proc_rv$cur_sample <- sname
            ss <- proc_rv$samp_stat; ss[[sname]] <- "running"; proc_rv$samp_stat <- ss
          }
          if (grepl("featureCounts|Generating count matrix|quant\\.sf", line, ignore.case=TRUE))
            proc_rv$cp_idx <- max(proc_rv$cp_idx, 5L)
          if (grepl("multiqc", line, ignore.case=TRUE))
            proc_rv$cp_idx <- max(proc_rv$cp_idx, 6L)
          if (grepl("Analysis completed", line, ignore.case=TRUE))
            proc_rv$cp_idx <- length(proc_rv$checkpoints)
        }
      }
    }
  })

  # ── UI dinamica: label genoma/transcriptoma ────────────────────────────────
  output$genome_label_ui <- renderUI({
    if (isTRUE(input$analysis_type == "alignment"))
      textInput("genome_file", "Genoma de referencia (FASTA)", placeholder = "/refs/ecoli.fa")
    else
      textInput("genome_file", "Transcriptoma de referencia (FASTA)",
                placeholder = "/refs/ecoli_transcriptome.fa")
  })

  # ── Deteccion de muestras con debounce ────────────────────────────────────
  input_dir_debounced <- debounce(reactive(input$input_dir), 600)
  samples_live <- reactive({
    input_dir_debounced()
    detect_samples(input$input_dir)
  })

  output$sample_preview_ui <- renderUI({
    dir_val <- input$input_dir %||% ""
    if (!nzchar(dir_val))
      return(div(class="text-muted", icon("folder-open"),
                 " Introduce un directorio de entrada."))
    if (!dir.exists(dir_val))
      return(div(class="text-danger", icon("circle-xmark"), " El directorio no existe."))
    samples <- samples_live()
    if (length(samples) == 0)
      return(div(class="alert alert-warning", icon("triangle-exclamation"),
                 " No se encontraron FASTQ (*_1.fastq.gz / *_R1.fastq.gz)."))

    # Buenas practicas: advertir nombres con caracteres problematicos
    bad <- bad_sample_chars(samples)
    r2_miss <- missing_r2(dir_val, samples)

    tagList(
      div(class="alert alert-info py-2",
          icon("circle-check"), sprintf(" %d muestra(s) detectadas.", length(samples))),
      if (length(bad) > 0)
        div(class="alert alert-warning py-2", icon("triangle-exclamation"),
            " Nombres con caracteres especiales (espacios, @, etc.) pueden causar",
            " problemas: ", tags$b(paste(bad, collapse=", "))),
      if (any(r2_miss))
        div(class="alert alert-danger py-2", icon("circle-xmark"),
            " Faltan R2 para: ", tags$b(paste(samples[r2_miss], collapse=", "))),
      tags$table(
        class="table table-sm table-striped",
        tags$thead(tags$tr(tags$th("Muestra"), tags$th("R1"), tags$th("R2"))),
        tags$tbody(lapply(seq_along(samples), function(i) {
          tags$tr(
            tags$td(samples[i]),
            tags$td(tags$span(style="color:green","\u2713")),
            tags$td(if (!r2_miss[[i]]) tags$span(style="color:green","\u2713")
                    else tags$span(style="color:red","\u2717 falta"))
          )
        }))
      )
    )
  })

  # ── Validaciones Tab 1 ────────────────────────────────────────────────────
  val_errors <- reactive({
    errs <- character(0)
    dir_in  <- input$input_dir  %||% ""
    dir_out <- input$output_dir %||% ""
    gf      <- input$genome_file %||% ""
    aln     <- input$analysis_type %||% "alignment"

    if (!nzchar(dir_in))           errs <- c(errs, "Directorio de FASTQs: campo vacio.")
    else if (!dir.exists(dir_in))  errs <- c(errs, "Directorio de FASTQs: no existe.")
    else {
      samps <- detect_samples(dir_in)
      if (length(samps) == 0)
        errs <- c(errs, "No se encontraron archivos FASTQ R1.")
      else {
        bad <- samps[missing_r2(dir_in, samps)]
        if (length(bad)) errs <- c(errs, paste("Faltan R2 para:", paste(bad, collapse=", ")))
      }
    }
    if (!nzchar(gf))              errs <- c(errs, "FASTA de referencia: campo vacio.")
    else if (!file.exists(gf))    errs <- c(errs, "FASTA de referencia: no existe.")
    if (aln == "alignment") {
      af <- input$annotation_file %||% ""
      if (!nzchar(af))            errs <- c(errs, "Anotacion GFF/GTF: campo vacio.")
      else if (!file.exists(af))  errs <- c(errs, "Anotacion GFF/GTF: no existe.")
    }
    if (!nzchar(dir_out))         errs <- c(errs, "Directorio de salida: campo vacio.")
    if (isTRUE(input$read_type == "se"))
      errs <- c(errs, "Single-end no soportado por workflow.sh actual.")
    if (!file.exists(workflow_path))
      errs <- c(errs, paste0("workflow.sh no encontrado en: ", workflow_path))
    errs
  })

  output$validation_ui <- renderUI({
    errs <- val_errors()
    if (length(errs) == 0)
      div(class="alert alert-success mb-0", icon("circle-check"),
          " Todos los campos correctos. Puedes continuar.")
    else
      div(class="alert alert-danger mb-0",
          tags$b(icon("circle-xmark"), " Corrige antes de continuar:"),
          tags$ul(class="mb-0 mt-1", lapply(errs, tags$li)))
  })

  observe({
    if (length(val_errors()) > 0) shinyjs::disable("btn_to_processing")
    else                           shinyjs::enable("btn_to_processing")
  })

  # ── Herramienta y anotacion efectivas ─────────────────────────────────────
  effective_tool <- reactive({
    if (isTRUE(input$analysis_type == "alignment")) "bowtie2"
    else input$pseudo_tool %||% "salmon"
  })

  effective_annotation <- reactive({
    if (isTRUE(input$analysis_type == "alignment")) input$annotation_file %||% ""
    else { ann <- input$annotation_file_pseudo %||% ""; if (nzchar(ann)) ann else "/dev/null" }
  })

  # ── Comando del workflow ──────────────────────────────────────────────────
  workflow_cmd <- reactive({
    req(input$input_dir, input$output_dir)
    gf <- input$genome_file; req(!is.null(gf) && nzchar(gf))
    sprintf(
      "bash %s --INPUT %s --OUTPUT %s --GENOME_FILE %s --ANNOTATION_FILE %s --ALIGNMENT_TYPE %s",
      shQuote(workflow_path),
      shQuote(input$input_dir), shQuote(input$output_dir),
      shQuote(gf), shQuote(effective_annotation()), shQuote(effective_tool())
    )
  })

  # ── Navegar Tab 1 -> Tab 2 ────────────────────────────────────────────────
  observeEvent(input$btn_to_processing, {
    req(length(val_errors()) == 0)
    process_unlocked(TRUE)
    samps <- detect_samples(input$input_dir)
    config_snap(list(
      analysis_type = input$analysis_type, tool = effective_tool(),
      input_dir  = input$input_dir, output_dir = input$output_dir,
      genome_file = input$genome_file %||% "", annotation = effective_annotation(),
      n_samples = length(samps), read_type = "Paired-end"
    ))
    log_text(paste(
      ts_log("=== Configuracion validada ==="),
      ts_log(paste("Analisis:", input$analysis_type, "/", effective_tool())),
      ts_log(paste("Muestras:", length(samps))),
      ts_log(paste("Entrada:", input$input_dir)),
      ts_log(paste("Salida:",  input$output_dir)), "",
      sep = "\n"
    ))
    updateNavbarPage(session, "main_nav", selected = "tab_process")
  })

  observeEvent(input$btn_back, {
    updateNavbarPage(session, "main_nav", selected = "tab_config")
  })

  # ── Carga desde resultados previos                    [v3-LOAD] ────────────
  observeEvent(input$btn_load_existing, {
    counts <- NULL
    # Option A: upload file
    if (!is.null(input$upload_counts)) {
      counts <- tryCatch(
        read.table(input$upload_counts$datapath, header=TRUE, row.names=1,
                   sep="\t", comment.char="#"),
        error = function(e)
          tryCatch(read.csv(input$upload_counts$datapath, row.names=1),
                   error=function(e2) NULL)
      )
    }
    # Option B: directory from previous run
    if (is.null(counts) && nzchar(input$prev_output_dir %||% "")) {
      for (tool in c("bowtie2","salmon","kallisto")) {
        counts <- tryCatch(load_counts_from_workflow(input$prev_output_dir, tool),
                           error=function(e) NULL)
        if (!is.null(counts)) break
      }
    }
    if (is.null(counts)) {
      showNotification("No se pudo cargar la matriz de conteos. Verifica el formato.",
                       type="error"); return()
    }
    counts <- round(counts)
    data_rv$count_matrix <- counts
    data_rv$counts_ready <- TRUE
    data_rv$source       <- "uploaded"

    # Cargar metadata si se proporciono
    if (!is.null(input$upload_meta)) {
      meta <- tryCatch(read.csv(input$upload_meta$datapath, stringsAsFactors=FALSE),
                       error=function(e) NULL)
      if (!is.null(meta)) data_rv$metadata <- meta
    }
    showNotification(
      sprintf("Datos cargados: %d genes x %d muestras. Ve a la pestana DEG.",
              nrow(counts), ncol(counts)), type="message", duration=8
    )
    updateNavbarPage(session, "main_nav", selected="tab_deg")
  })

  # ── Contenido Tab 2 ───────────────────────────────────────────────────────
  output$tab2_content <- renderUI({
    if (!process_unlocked())
      return(div(class="alert alert-info mt-4", icon("arrow-left-long"),
                 " Completa la pestana 1 antes de continuar."))
    cfg <- config_snap()
    if (length(cfg) == 0) return(NULL)

    r1_files <- list.files(cfg$input_dir,
                           pattern="(_1\\.fastq\\.gz$|_R1\\.fastq\\.gz$)",
                           full.names=TRUE, ignore.case=TRUE)
    total_sz <- if (length(r1_files)) fmt_bytes(sum(file.info(r1_files)$size,na.rm=TRUE)*2) else "—"

    tagList(
      # ── Fila superior: resumen + progreso ─────────────────────────────
      layout_columns(
        col_widths = c(4, 8),

        # Izquierda: resumen y controles
        card(
          card_header(icon("list-check"), " Resumen del analisis"),
          tags$dl(class="row small mb-1",
            tags$dt(class="col-6","Tipo:"),
            tags$dd(class="col-6", if(cfg$analysis_type=="alignment")
              "Alineamiento (Bowtie2)" else paste0("Pseudoalineamiento (",cfg$tool,")")),
            tags$dt(class="col-6","Muestras:"), tags$dd(class="col-6", cfg$n_samples),
            tags$dt(class="col-6","Lectura:"),  tags$dd(class="col-6", cfg$read_type),
            tags$dt(class="col-6","Tama\u00f1o est.:"), tags$dd(class="col-6", total_sz),
            tags$dt(class="col-6","Entrada:"),
            tags$dd(class="col-6", tags$code(class="small", cfg$input_dir)),
            tags$dt(class="col-6","Salida:"),
            tags$dd(class="col-6", tags$code(class="small", cfg$output_dir))
          ),
          hr(),
          tags$details(
            tags$summary(tags$small(icon("terminal")," Ver comando")),
            tags$pre(class="small mt-1",
                     style="white-space:pre-wrap;word-break:break-all;background:#f8f9fa;padding:8px;",
                     textOutput("cmd_preview_text", inline=TRUE))
          ),
          hr(),
          div(style="display:flex;gap:8px;align-items:center;",
              actionButton("btn_back", tagList(icon("arrow-left")," Volver"), class="btn-secondary"),
              actionButton("run_btn",  tagList(icon("play")," Ejecutar workflow"), class="btn-success btn-lg")
          )
        ),

        # Derecha: progreso detallado                 [v3-PROG]
        card(
          card_header(icon("gauge"), " Progreso en tiempo real"),

          # Metricas de tiempo y muestras
          div(class="d-flex gap-3 mb-3 p-2 rounded",
              style="background:#f0f4f8; font-size:.9rem;",
              div(icon("clock"), " Transcurrido: ",
                  tags$b(textOutput("elapsed_text", inline=TRUE))),
              div(icon("hourglass-half"), " ETA: ",
                  tags$b(textOutput("eta_text", inline=TRUE))),
              div(icon("vials"), " Muestras: ",
                  tags$b(textOutput("samp_prog_text", inline=TRUE)))
          ),

          # Lista de checkpoints
          uiOutput("checkpoints_ui"),
          hr(),
          # Tabla de estado por muestra
          uiOutput("sample_status_ui")
        )
      ),

      # ── Fila inferior: log ancho completo            [v3-LAYOUT] ───────
      card(
        card_header(
          icon("terminal"), " Log de ejecucion",
          div(class="float-end",
              actionButton("refresh_btn", icon("rotate"),
                           class="btn-sm btn-outline-secondary",
                           title="Actualizar lista de archivos"))
        ),
        verbatimTextOutput("run_log"),
        max_height = "220px",
        style = "overflow-y:auto; width:100%;"
      )
    )
  })

  # ── Outputs de progreso (dependen de poll_timer para actualizarse) ─────────
  output$elapsed_text <- renderText({
    poll_timer()
    if (is.null(proc_rv$start_time)) return("—")
    fmt_elapsed(as.numeric(difftime(Sys.time(), proc_rv$start_time, units="secs")))
  })

  output$eta_text <- renderText({
    poll_timer()
    if (is.null(proc_rv$start_time) || !proc_rv$running) return("—")
    n_done  <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total == 0 || n_done == 0) return("calculando...")
    elapsed <- as.numeric(difftime(Sys.time(), proc_rv$start_time, units="secs"))
    # Estimacion: 15% para el setup, 85% para las muestras
    frac_done <- 0.15 + 0.85 * n_done / n_total
    remaining <- elapsed / frac_done - elapsed
    if (remaining < 0) return("finalizando...")
    paste0("~", fmt_elapsed(remaining))
  })

  output$samp_prog_text <- renderText({
    poll_timer()
    n_done  <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_run   <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "running"))
    n_total <- proc_rv$n_total
    if (n_total == 0) return("—")
    sprintf("%d / %d completadas", n_done, n_total)
  })

  output$checkpoints_ui <- renderUI({
    poll_timer()
    cps    <- proc_rv$checkpoints
    cp_idx <- proc_rv$cp_idx
    if (length(cps) == 0) return(NULL)
    tags$ul(class="list-unstyled mb-0",
      lapply(seq_along(cps), function(i) {
        if (i < cp_idx)
          tags$li(style="color:green;", tags$span("\u2713 "), cps[i])
        else if (i == cp_idx)
          tags$li(style="color:#d97706; font-weight:600;",
                  tags$span("\u27f3 "), cps[i])
        else
          tags$li(style="color:#9ca3af;", tags$span("\u25f7 "), cps[i])
      })
    )
  })

  output$sample_status_ui <- renderUI({
    poll_timer()
    ss <- proc_rv$samp_stat
    if (length(ss) == 0) return(NULL)
    tags$table(class="table table-sm",
      tags$thead(tags$tr(tags$th("Muestra"), tags$th("Estado"))),
      tags$tbody(lapply(names(ss), function(s) {
        st <- ss[[s]]
        icon_el <- switch(st,
          done    = tags$span(style="color:green;",  "\u2713 Completada"),
          running = tags$span(style="color:#d97706;", "\u27f3 Procesando..."),
          tags$span(style="color:#9ca3af;", "\u25f7 Pendiente")
        )
        tags$tr(tags$td(s), tags$td(icon_el))
      }))
    )
  })

  output$cmd_preview_text <- renderText({
    tryCatch(workflow_cmd(),
             error=function(e) "Faltan campos para generar el comando.")
  })

  # ── Ejecutar workflow ──────────────────────────────────────────────────────
  observeEvent(input$run_btn, {
    errs <- val_errors()
    if (length(errs) > 0) {
      showNotification(paste(errs, collapse="\n"), type="error", duration=8); return()
    }
    if (!file.exists(workflow_path)) {
      showNotification("workflow.sh no encontrado.", type="error"); return()
    }
    cmd <- tryCatch(workflow_cmd(), error=function(e) NULL)
    if (is.null(cmd)) {
      showNotification("No se puede construir el comando.", type="error"); return()
    }

    # Snapshot de parametros para Tab 3
    samps <- detect_samples(input$input_dir)
    run_params_rv(list(
      analysis_type = input$analysis_type, tool = effective_tool(),
      input_dir = input$input_dir, output_dir = input$output_dir,
      genome_file = input$genome_file %||% "", annotation_file = effective_annotation(),
      n_samples = length(samps), read_type = "Paired-end",
      started_at = Sys.time(),
      r_version = paste(R.version$major, R.version$minor, sep=".")
    ))

    analysis_done(FALSE)
    data_rv$counts_ready <- FALSE
    shinyjs::disable("run_btn")
    on.exit(shinyjs::enable("run_btn"), add=TRUE)

    # Checkpoints segun tipo de analisis
    cps <- if (input$analysis_type == "alignment")
      c("Construyendo indice Bowtie2",
        "Control de calidad inicial (FastQC)",
        "Alineamiento de muestras (Bowtie2 + fastp)",
        "Procesando BAM (samtools sort + index)",
        "Conteos por gen (featureCounts)",
        "Control de calidad post-trimming (FastQC)",
        "Informe global (MultiQC)")
    else
      c(paste0("Construyendo indice (", effective_tool(), ")"),
        "Control de calidad inicial (FastQC)",
        paste0("Cuantificacion de muestras (", effective_tool(), " + fastp)"),
        "Importando cuantificaciones",
        "Matriz de conteos",
        "Control de calidad post-trimming (FastQC)",
        "Informe global (MultiQC)")

    # Inicializar proc_rv
    proc_rv$checkpoints <- cps
    proc_rv$cp_idx      <- 0L
    proc_rv$n_total     <- length(samps)
    proc_rv$samp_stat   <- setNames(as.list(rep("pending", length(samps))), samps)
    proc_rv$cur_sample  <- NULL
    proc_rv$start_time  <- Sys.time()

    log_text(paste0(log_text(),
      ts_log("=== Iniciando analisis ==="), "\n",
      ts_log(paste0("Lanzando: ", cmd)), "\n"))

    if (HAS_PROCESSX) {
      # ── Ejecucion no bloqueante con processx          [v3-PROC] ─────────
      proc <- tryCatch(
        processx::process$new("bash", c("-lc", cmd), stdout="|", stderr="|",
                              env=c("PATH"=Sys.getenv("PATH"))),
        error=function(e) { showNotification(conditionMessage(e), type="error"); NULL }
      )
      if (!is.null(proc)) {
        proc_rv$proc    <- proc
        proc_rv$running <- TRUE
        shinyjs::enable("run_btn")  # permitir nuevo intento si falla
      }
    } else {
      # ── Fallback bloqueante con system2 ────────────────────────────────
      withProgress(message="Ejecutando analisis RNA-seq...", value=0, {
        for (i in 1:2) {
          incProgress(1/9, detail=cps[i])
          Sys.sleep(0.1)
        }
        incProgress(1/9, detail=cps[3])
        result <- tryCatch(
          system2("bash", c("-lc", cmd), stdout=TRUE, stderr=TRUE),
          error=function(e) structure(conditionMessage(e), status=1L)
        )
        exit_status <- attr(result, "status") %||% 0L
        log_text(paste0(log_text(), paste(result, collapse="\n"), "\n",
                        ts_log(paste0("Codigo de salida: ", exit_status)), "\n"))
        for (i in 4:7) { incProgress(1/9, detail=cps[i]); Sys.sleep(0.1) }

        # Cargar conteos y archivos
        p <- run_params_rv()
        if (exit_status == 0) {
          counts <- tryCatch(load_counts_from_workflow(p$output_dir, p$tool), error=function(e) NULL)
          if (!is.null(counts)) {
            data_rv$count_matrix <- counts
            data_rv$counts_ready <- TRUE
          }
          files <- list.files(p$output_dir, recursive=TRUE, full.names=FALSE)
          if (length(files)) {
            sz <- file.info(file.path(p$output_dir, files))$size
            output_files_rv(data.frame(
              Archivo=files, `Tamano`=sapply(sz,fmt_bytes),
              stringsAsFactors=FALSE, check.names=FALSE
            ))
          }
          log_text(paste0(log_text(), ts_log("=== Analisis finalizado OK ==="), "\n"))
          analysis_done(TRUE)
          showNotification("Workflow finalizado correctamente.", type="message")
          updateNavbarPage(session, "main_nav", selected="tab_results")
        } else {
          showNotification(paste0("Error (codigo ",exit_status,"). Revisa el log."),
                           type="error", duration=12)
        }
      })
    }
  })

  # ── Log y refresh de archivos ─────────────────────────────────────────────
  output$run_log <- renderText({ log_text() })

  observeEvent(input$refresh_btn, {
    out_dir <- input$output_dir %||% ""
    if (!nzchar(out_dir) || !dir.exists(out_dir)) {
      output_files_rv(data.frame(Archivo="Directorio no existe.", `Tamano`="—",
                                 stringsAsFactors=FALSE, check.names=FALSE)); return()
    }
    files <- list.files(out_dir, recursive=TRUE, full.names=FALSE)
    if (!length(files)) {
      output_files_rv(data.frame(Archivo="Sin archivos.", `Tamano`="—",
                                 stringsAsFactors=FALSE, check.names=FALSE))
    } else {
      sz <- file.info(file.path(out_dir, files))$size
      output_files_rv(data.frame(Archivo=files, `Tamano`=sapply(sz,fmt_bytes),
                                 stringsAsFactors=FALSE, check.names=FALSE))
    }
  })

  # ── Contenido Tab 3 ───────────────────────────────────────────────────────
  output$tab3_content <- renderUI({
    if (!analysis_done())
      return(div(class="alert alert-info mt-4", icon("clock"),
                 " Ejecuta el workflow en la pestana 2 para ver los resultados."))
    p <- run_params_rv()
    tagList(
      layout_columns(col_widths=c(3,3,3,3),
        value_box(title="Tipo de analisis",
                  value=if(p$analysis_type=="alignment") "Alineamiento" else "Pseudoalineamiento",
                  showcase=icon("dna"), theme="primary"),
        value_box(title="Muestras", value=p$n_samples, showcase=icon("microscope"), theme="success"),
        value_box(title="Archivos", value=nrow(output_files_rv()), showcase=icon("file-lines"), theme="info"),
        value_box(title="Herramienta", value=toupper(p$tool), showcase=icon("terminal"), theme="secondary")
      ),
      card(card_header(icon("folder")," Archivos generados"),
        layout_columns(col_widths=c(9,3),
          DTOutput("output_files_table"),
          div(class="d-flex flex-column gap-2",
              downloadButton("download_log",  tagList(icon("download")," Log"),
                             class="btn-outline-secondary w-100"),
              downloadButton("download_filelist", tagList(icon("download")," CSV archivos"),
                             class="btn-outline-secondary w-100"),
              downloadButton("download_repro", tagList(icon("recycle")," Reproducibilidad"),
                             class="btn-outline-info w-100"),
              if (data_rv$counts_ready)
                actionButton("go_to_deg", tagList(icon("chart-line")," Ir a DEG"),
                             class="btn-primary w-100")
          )
        )
      ),
      card(card_header(icon("recycle")," Reproducibilidad"), uiOutput("repro_ui"))
    )
  })

  output$output_files_table <- renderDT({
    datatable(output_files_rv(), options=list(pageLength=20, scrollX=TRUE, dom="ftip"),
              rownames=FALSE)
  })

  repro_text_rv <- reactive({ build_repro_text(run_params_rv(), workflow_path) })
  output$repro_ui <- renderUI({
    txt <- repro_text_rv(); if (!nzchar(txt)) return(NULL)
    tags$pre(class="small", style="background:#f8f9fa;padding:14px;border-radius:4px;", txt)
  })

  observeEvent(input$go_to_deg, {
    updateNavbarPage(session, "main_nav", selected="tab_deg")
  })

  output$download_log <- downloadHandler(
    filename=function() paste0("rnaseq_log_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt"),
    content=function(f) writeLines(log_text(), f))
  output$download_filelist <- downloadHandler(
    filename=function() paste0("output_files_",format(Sys.time(),"%Y%m%d_%H%M%S"),".csv"),
    content=function(f) write.csv(output_files_rv(), f, row.names=FALSE))
  output$download_repro <- downloadHandler(
    filename=function() paste0("rnaseq_repro_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt"),
    content=function(f) writeLines(repro_text_rv(), f))

  # ╔══════════════════════════════════════════════════════════════════════════
  # ║  TAB 4 — DEG                                             [v3-DEG]       ║
  # ╚══════════════════════════════════════════════════════════════════════════

  output$tab4_gate_ui <- renderUI({
    if (!data_rv$counts_ready)
      div(class="alert alert-info mt-3", icon("clock"),
          " Necesitas datos de conteos. Ejecuta el workflow (Tab 2) o carga",
          " resultados previos (Tab 1 \u2192 Cargar resultados previos).")
  })

  # Resumen de la matriz de conteos
  output$counts_summary_ui <- renderUI({
    req(data_rv$counts_ready)
    mat <- data_rv$count_matrix
    tags$dl(class="row small mb-0",
      tags$dt(class="col-6","Genes:"),    tags$dd(class="col-6", nrow(mat)),
      tags$dt(class="col-6","Muestras:"), tags$dd(class="col-6", ncol(mat)),
      tags$dt(class="col-6","Fuente:"),   tags$dd(class="col-6", data_rv$source)
    )
  })

  # Tabla de asignacion condicion por muestra
  output$condition_table_ui <- renderUI({
    req(data_rv$counts_ready)
    mat   <- data_rv$count_matrix
    conds <- trimws(strsplit(input$condition_labels %||% "", ",")[[1]])
    conds <- conds[nzchar(conds)]
    if (length(conds) < 2)
      return(div(class="alert alert-warning py-2",
                 "Define al menos 2 condiciones (separadas por coma)."))
    # Intentar prerellenar desde metadata si existe
    pre <- if (!is.null(data_rv$metadata)) {
      meta <- data_rv$metadata
      setNames(meta$condition[match(colnames(mat), meta$sample)], colnames(mat))
    } else setNames(rep(conds[1], ncol(mat)), colnames(mat))

    tags$table(class="table table-sm",
      tags$thead(tags$tr(tags$th("Muestra"), tags$th("Condicion"))),
      tags$tbody(lapply(colnames(mat), function(s) {
        sel <- pre[[s]] %||% conds[1]
        sel <- if (sel %in% conds) sel else conds[1]
        tags$tr(tags$td(s),
                tags$td(selectInput(paste0("cond_", s), label=NULL,
                                    choices=conds, selected=sel, width="100%")))
      }))
    )
  })

  # Actualizar selects de referencia y tratamiento segun condiciones
  observe({
    conds <- trimws(strsplit(input$condition_labels %||% "", ",")[[1]])
    conds <- conds[nzchar(conds)]
    if (length(conds) >= 2) {
      updateSelectInput(session, "deg_ref_cond",   choices=conds, selected=conds[1])
      updateSelectInput(session, "deg_treat_cond", choices=conds, selected=conds[2])
    }
  })

  # ── Ejecutar DESeq2 ────────────────────────────────────────────────────────
  observeEvent(input$run_deg, {
    if (!HAS_DESEQ2) {
      showNotification("Instala DESeq2: BiocManager::install('DESeq2')",
                       type="error", duration=10); return()
    }
    req(data_rv$counts_ready)
    counts <- data_rv$count_matrix
    samples <- colnames(counts)
    conds <- vapply(samples, function(s) input[[paste0("cond_", s)]] %||% "", character(1))
    if (any(!nzchar(conds))) {
      showNotification("Asigna condicion a todas las muestras.", type="error"); return()
    }
    ref_cond   <- input$deg_ref_cond
    treat_cond <- input$deg_treat_cond
    if (ref_cond == treat_cond) {
      showNotification("Referencia y tratamiento deben ser diferentes.", type="error"); return()
    }
    meta <- data.frame(sample=samples, condition=conds, stringsAsFactors=FALSE)
    meta$condition <- factor(meta$condition, levels=c(ref_cond, treat_cond))

    data_rv$deg_done <- FALSE

    withProgress(message="Ejecutando DESeq2...", value=0, {
      incProgress(0.1, detail="Filtrando genes de baja expresion...")
      cf <- counts[rowSums(counts >= 1) >= 2, , drop=FALSE]
      cf <- cf[rowSums(cf) >= input$deg_min_counts, , drop=FALSE]
      if (nrow(cf) == 0) {
        showNotification("Sin genes tras filtrado. Reduce el umbral.", type="error"); return()
      }

      incProgress(0.2, detail="Construyendo DESeqDataSet...")
      meta2 <- meta[match(colnames(cf), meta$sample), ]
      dds <- tryCatch(
        DESeq2::DESeqDataSetFromMatrix(countData=cf, colData=meta2, design= ~condition),
        error=function(e) { showNotification(conditionMessage(e), type="error"); NULL }
      )
      if (is.null(dds)) return()

      incProgress(0.5, detail="Ajustando modelo (puede tardar)...")
      dds <- tryCatch(
        DESeq2::DESeq(dds),
        error=function(e) { showNotification(conditionMessage(e), type="error"); NULL }
      )
      if (is.null(dds)) return()

      incProgress(0.85, detail="Extrayendo resultados...")
      res <- DESeq2::results(dds, contrast=c("condition", treat_cond, ref_cond),
                             alpha=input$deg_fdr)
      df_res            <- as.data.frame(res)
      df_res$gene       <- rownames(df_res)
      data_rv$deg_dds   <- dds
      data_rv$deg_res   <- df_res
      data_rv$deg_done  <- TRUE

      n_sig <- sum(!is.na(df_res$padj) & df_res$padj < input$deg_fdr &
                     abs(df_res$log2FoldChange) > input$deg_lfc, na.rm=TRUE)
      incProgress(1.0, detail="Completado.")
      showNotification(sprintf("DESeq2 OK: %d genes DEG (FDR<%.3f, |LFC|>%.1f)",
                               n_sig, input$deg_fdr, input$deg_lfc), type="message")
    })
  })

  output$deg_status_ui <- renderUI({
    if (!data_rv$deg_done) return(NULL)
    df  <- data_rv$deg_res
    n_up   <- sum(!is.na(df$padj) & df$padj<input$deg_fdr & df$log2FoldChange>input$deg_lfc, na.rm=TRUE)
    n_down <- sum(!is.na(df$padj) & df$padj<input$deg_fdr & df$log2FoldChange < -input$deg_lfc, na.rm=TRUE)
    div(class="alert alert-success mt-2 py-2",
        icon("circle-check"),
        sprintf(" DESeq2 completado: %d UP, %d DOWN, %d ns",
                n_up, n_down, nrow(df) - n_up - n_down))
  })

  # ── PCA ────────────────────────────────────────────────────────────────────
  output$pca_plot <- renderPlotly({
    req(data_rv$deg_done, HAS_PLOTLY)
    dds <- data_rv$deg_dds; req(!is.null(dds))
    vst_d <- tryCatch(DESeq2::vst(dds, blind=TRUE), error=function(e) NULL)
    if (is.null(vst_d)) return(plotly_empty())

    pca_data <- DESeq2::plotPCA(vst_d, intgroup="condition", returnData=TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

    x_col <- input$pca_x %||% "PC1"; y_col <- input$pca_y %||% "PC2"
    pc_cols <- grep("^PC", names(pca_data), value=TRUE)
    if (!x_col %in% pc_cols) x_col <- pc_cols[1]
    if (!y_col %in% pc_cols) y_col <- pc_cols[min(2, length(pc_cols))]
    x_v <- pct_var[which(pc_cols == x_col) %||% 1] %||% 0
    y_v <- pct_var[which(pc_cols == y_col) %||% 2] %||% 0

    p <- ggplot(pca_data, aes(x=.data[[x_col]], y=.data[[y_col]],
                               color=condition, label=name,
                               text=paste0("Muestra: ", name, "\nCondicion: ", condition))) +
      geom_point(size=5, alpha=0.85) +
      xlab(sprintf("%s  (%.1f%% varianza)", x_col, x_v)) +
      ylab(sprintf("%s  (%.1f%% varianza)", y_col, y_v)) +
      theme_bw(base_size=13) +
      scale_color_brewer(palette="Set1")

    if (isTRUE(input$pca_labels) && HAS_GGREPEL)
      p <- p + ggrepel::geom_text_repel(show.legend=FALSE, size=3.5)

    ggplotly(p, tooltip="text")
  })

  output$dl_pca <- downloadHandler(
    filename=function() paste0("PCA_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"),
    content=function(f) {
      req(data_rv$deg_done)
      vst_d <- DESeq2::vst(data_rv$deg_dds, blind=TRUE)
      pca_data <- DESeq2::plotPCA(vst_d, intgroup="condition", returnData=TRUE)
      pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
      p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=condition, label=name)) +
        geom_point(size=5, alpha=0.85) +
        xlab(paste0("PC1 (",pct_var[1],"% varianza)")) +
        ylab(paste0("PC2 (",pct_var[2],"% varianza)")) +
        theme_bw(base_size=13) + scale_color_brewer(palette="Set1")
      if (HAS_GGREPEL) p <- p + ggrepel::geom_text_repel(show.legend=FALSE, size=3)
      ggsave(f, p, width=8, height=6, dpi=150)
    })

  # ── Volcano Plot ───────────────────────────────────────────────────────────
  output$volcano_plot <- renderPlotly({
    req(data_rv$deg_done, HAS_PLOTLY)
    df <- data_rv$deg_res; req(!is.null(df))
    df$nlp <- -log10(pmax(df$padj, 1e-300))
    df$cat  <- "ns"
    df$cat[!is.na(df$padj) & df$padj < input$vol_fdr & df$log2FoldChange >  input$vol_lfc] <- "up"
    df$cat[!is.na(df$padj) & df$padj < input$vol_fdr & df$log2FoldChange < -input$vol_lfc] <- "down"
    df$tooltip <- paste0("Gen: ", df$gene,
                         "\nLFC: ", round(df$log2FoldChange, 3),
                         "\npadj: ", formatC(df$padj, format="e", digits=2),
                         "\nCategoria: ", df$cat)
    colors <- c(up="#D73027", down="#4575B4", ns="#AAAAAA")

    p <- ggplot(df, aes(x=log2FoldChange, y=nlp, color=cat, text=tooltip)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(values=colors, name="") +
      geom_vline(xintercept=c(-input$vol_lfc, input$vol_lfc),
                 linetype="dashed", color="gray50", alpha=0.7) +
      geom_hline(yintercept=-log10(input$vol_fdr),
                 linetype="dashed", color="gray50", alpha=0.7) +
      xlab("log2 Fold Change") + ylab("-log10(p.adj)") +
      theme_bw(base_size=13)

    # Etiquetas top N genes
    n_top <- input$vol_top %||% 0
    if (n_top > 0 && HAS_GGREPEL) {
      sig_df <- df[df$cat != "ns" & !is.na(df$cat), ]
      if (nrow(sig_df) > 0) {
        top_df <- head(sig_df[order(sig_df$padj), ], n_top)
        p <- p + ggrepel::geom_text_repel(
          data=top_df, aes(label=gene), show.legend=FALSE, size=3, max.overlaps=20)
      }
    }
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.1))
  })

  output$dl_volcano <- downloadHandler(
    filename=function() paste0("Volcano_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"),
    content=function(f) {
      req(data_rv$deg_done)
      df <- data_rv$deg_res
      df$nlp <- -log10(pmax(df$padj, 1e-300))
      df$cat <- "ns"
      df$cat[!is.na(df$padj) & df$padj<input$vol_fdr & df$log2FoldChange>input$vol_lfc]  <- "up"
      df$cat[!is.na(df$padj) & df$padj<input$vol_fdr & df$log2FoldChange< -input$vol_lfc] <- "down"
      p <- ggplot(df, aes(x=log2FoldChange, y=nlp, color=cat)) +
        geom_point(alpha=0.6, size=1.5) +
        scale_color_manual(values=c(up="#D73027",down="#4575B4",ns="#AAAAAA"),name="") +
        geom_vline(xintercept=c(-input$vol_lfc,input$vol_lfc), linetype="dashed", color="gray50") +
        geom_hline(yintercept=-log10(input$vol_fdr), linetype="dashed", color="gray50") +
        xlab("log2 Fold Change") + ylab("-log10(p.adj)") + theme_bw(base_size=13)
      ggsave(f, p, width=8, height=6, dpi=150)
    })

  # ── Heatmap ────────────────────────────────────────────────────────────────
  output$heatmap_plot <- renderPlot({
    req(data_rv$deg_done, HAS_PHEATMAP)
    res <- data_rv$deg_res; dds <- data_rv$deg_dds
    req(!is.null(res), !is.null(dds))

    res_sig <- res[!is.na(res$padj) & res$padj < input$hm_fdr, , drop=FALSE]
    if (nrow(res_sig) == 0) {
      plot.new(); text(0.5, 0.5, "Sin genes significativos con el FDR actual.", cex=1.2)
      return()
    }
    top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), input$hm_top_n)

    vst_d <- tryCatch(DESeq2::vst(dds, blind=FALSE), error=function(e) NULL)
    if (is.null(vst_d)) return()
    mat <- DESeq2::assay(vst_d)[top_genes, , drop=FALSE]
    if (isTRUE(input$hm_scale == "row")) mat <- t(scale(t(mat)))

    col_ann <- data.frame(Condicion=as.character(dds$condition), row.names=colnames(mat))
    row_ann <- data.frame(
      DEG = ifelse(res[top_genes,"log2FoldChange"] > 0, "Up", "Down"),
      row.names = top_genes
    )
    ann_colors <- list(DEG=c(Up="#D73027", Down="#4575B4"))

    pal_name <- input$hm_palette %||% "RdBu"
    pal_raw  <- tryCatch(RColorBrewer::brewer.pal(11, pal_name), error=function(e) NULL)
    pal <- if (!is.null(pal_raw)) colorRampPalette(rev(pal_raw))(100) else
      colorRampPalette(c("#2166AC","white","#B2182B"))(100)

    clust_method <- input$hm_clust %||% "ward.D2"
    do_clust <- clust_method != "none"
    if (!do_clust) clust_method <- "complete"

    pheatmap(mat, color=pal,
             annotation_col=col_ann, annotation_row=row_ann, annotation_colors=ann_colors,
             clustering_method=clust_method, cluster_rows=do_clust, cluster_cols=do_clust,
             show_rownames = isTRUE(input$hm_rownames) || nrow(mat) <= 30,
             show_colnames = TRUE,
             fontsize_row  = max(5, min(10, 300/nrow(mat))),
             main = sprintf("Top %d genes DEG  (FDR < %.3f)", nrow(mat), input$hm_fdr),
             silent = TRUE)
  }, height=function() max(400, min(1400, (input$hm_top_n %||% 50) * 14 + 120)))

  output$dl_heatmap <- downloadHandler(
    filename=function() paste0("Heatmap_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"),
    content=function(f) {
      req(data_rv$deg_done, HAS_PHEATMAP)
      res <- data_rv$deg_res; dds <- data_rv$deg_dds
      res_sig   <- res[!is.na(res$padj) & res$padj < input$hm_fdr, , drop=FALSE]
      top_genes <- head(rownames(res_sig[order(res_sig$padj),]), input$hm_top_n)
      vst_d     <- DESeq2::vst(dds, blind=FALSE)
      mat       <- DESeq2::assay(vst_d)[top_genes, , drop=FALSE]
      if (isTRUE(input$hm_scale=="row")) mat <- t(scale(t(mat)))
      col_ann   <- data.frame(Condicion=as.character(dds$condition), row.names=colnames(mat))
      pal_raw   <- tryCatch(RColorBrewer::brewer.pal(11,input$hm_palette%||%"RdBu"),error=function(e)NULL)
      pal <- if(!is.null(pal_raw)) colorRampPalette(rev(pal_raw))(100) else
        colorRampPalette(c("#2166AC","white","#B2182B"))(100)
      png(f, width=1800, height=max(800, nrow(mat)*22+200), res=150)
      pheatmap(mat, color=pal, annotation_col=col_ann,
               clustering_method=input$hm_clust%||%"ward.D2",
               show_rownames=TRUE, silent=TRUE)
      dev.off()
    })

  # ── Tabla estadistica DEG ─────────────────────────────────────────────────
  output$deg_table <- renderDT({
    req(data_rv$deg_done)
    df <- data_rv$deg_res[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    df <- df[order(df$padj %||% Inf, na.last=TRUE), ]
    for (col in c("baseMean","log2FoldChange","lfcSE","stat"))
      df[[col]] <- round(df[[col]], 3)
    for (col in c("pvalue","padj"))
      df[[col]] <- formatC(df[[col]], format="e", digits=3)
    datatable(df, rownames=FALSE,
              options=list(pageLength=20, scrollX=TRUE, dom="ftip"),
              filter="top") %>%
      formatStyle("log2FoldChange",
                  color = styleInterval(c(-input$deg_lfc, input$deg_lfc),
                                        c("#4575B4","#888","#D73027")))
  })

  output$dl_deg_csv <- downloadHandler(
    filename=function() paste0("DEG_full_",format(Sys.time(),"%Y%m%d_%H%M%S"),".csv"),
    content=function(f) { req(data_rv$deg_done); write.csv(data_rv$deg_res, f, row.names=FALSE) })
  output$dl_deg_sig <- downloadHandler(
    filename=function() paste0("DEG_sig_",format(Sys.time(),"%Y%m%d_%H%M%S"),".csv"),
    content=function(f) {
      req(data_rv$deg_done)
      df <- data_rv$deg_res
      sig <- df[!is.na(df$padj) & df$padj<input$deg_fdr & abs(df$log2FoldChange)>input$deg_lfc, ]
      write.csv(sig, f, row.names=FALSE)
    })

  # ╔══════════════════════════════════════════════════════════════════════════
  # ║  TAB 5 — ANALISIS FUNCIONAL                          [v3-FUNC]          ║
  # ╚══════════════════════════════════════════════════════════════════════════

  output$tab5_gate_ui <- renderUI({
    if (!data_rv$deg_done)
      div(class="alert alert-info mt-3", icon("clock"),
          " Ejecuta el analisis DEG (Tab 4) antes de continuar.")
  })

  # ── ORA ───────────────────────────────────────────────────────────────────
  observeEvent(input$run_ora, {
    if (!HAS_CLUSTERP) {
      showNotification("Instala: BiocManager::install('clusterProfiler')",
                       type="error", duration=10); return()
    }
    req(data_rv$deg_done)
    res <- data_rv$deg_res
    sig <- rownames(res[!is.na(res$padj) & res$padj < input$ora_fdr &
                          abs(res$log2FoldChange) > input$ora_lfc, ])
    if (!length(sig)) {
      showNotification("Sin genes significativos con los umbrales actuales.", type="warning"); return()
    }

    withProgress(message="Ejecutando ORA...", value=0, {
      db     <- input$func_db %||% "kegg"
      result <- NULL

      if (db == "kegg") {
        incProgress(0.3, detail="Consultando KEGG (requiere internet)...")
        result <- tryCatch(
          clusterProfiler::enrichKEGG(
            gene=sig, organism=input$kegg_org %||% "eco",
            pvalueCutoff=input$ora_pval, qvalueCutoff=input$ora_qval),
          error=function(e){ showNotification(conditionMessage(e), type="error", duration=10); NULL })
      } else {
        ont <- switch(db, gobp="BP", gomf="MF", gocc="CC", "BP")
        org_pkg <- input$orgdb_pkg %||% ""
        if (!nzchar(org_pkg) || !pkg_ok(org_pkg)) {
          showNotification(
            paste0("Instala el paquete OrgDb: BiocManager::install('", org_pkg, "')"),
            type="error", duration=10); return()
        }
        org_db <- tryCatch(getExportedValue(org_pkg, org_pkg), error=function(e) NULL)
        if (is.null(org_db)) {
          suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
          org_db <- get(org_pkg)
        }
        incProgress(0.3, detail=paste("enrichGO", ont, "..."))
        result <- tryCatch(
          clusterProfiler::enrichGO(
            gene=sig, OrgDb=org_db, keyType=input$orgdb_keytype %||% "SYMBOL",
            ont=ont, pAdjustMethod="BH",
            pvalueCutoff=input$ora_pval, qvalueCutoff=input$ora_qval, readable=FALSE),
          error=function(e){ showNotification(conditionMessage(e), type="error", duration=10); NULL })
      }

      incProgress(1.0)
      data_rv$ora_res <- result
      if (!is.null(result))
        showNotification(sprintf("ORA: %d terminos enriquecidos.", nrow(as.data.frame(result))),
                         type="message")
    })
  })

  output$ora_dot_plot <- renderPlotly({
    req(!is.null(data_rv$ora_res), HAS_PLOTLY)
    df <- as.data.frame(data_rv$ora_res)
    if (nrow(df) == 0) return(plotly_empty(type="scatter", mode="markers") %>%
                                 layout(title="Sin terminos enriquecidos."))
    p <- make_enrichment_dotplot(df, title="ORA — Top 20")
    ggplotly(p, tooltip="text")
  })

  output$ora_table <- renderDT({
    req(!is.null(data_rv$ora_res))
    df <- as.data.frame(data_rv$ora_res)
    if (nrow(df) == 0) return(datatable(data.frame(Mensaje="Sin terminos enriquecidos.")))
    df <- df[, intersect(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count"),
                         names(df)), drop=FALSE]
    df$pvalue   <- formatC(df$pvalue,   format="e", digits=3)
    df$p.adjust <- formatC(df$p.adjust, format="e", digits=3)
    datatable(df, rownames=FALSE, options=list(pageLength=10, scrollX=TRUE, dom="ftip"))
  })

  output$dl_ora_csv <- downloadHandler(
    filename=function() paste0("ORA_",format(Sys.time(),"%Y%m%d_%H%M%S"),".csv"),
    content=function(f) {
      req(!is.null(data_rv$ora_res))
      write.csv(as.data.frame(data_rv$ora_res), f, row.names=FALSE)
    })

  # ── GSEA ──────────────────────────────────────────────────────────────────
  observeEvent(input$run_gsea, {
    if (!HAS_CLUSTERP) {
      showNotification("Instala: BiocManager::install('clusterProfiler')",
                       type="error", duration=10); return()
    }
    req(data_rv$deg_done)
    res <- data_rv$deg_res
    rank_col  <- if (isTRUE(input$gsea_rank == "stat")) "stat" else "log2FoldChange"
    gene_list <- res[[rank_col]]
    names(gene_list) <- res$gene %||% rownames(res)
    gene_list <- sort(gene_list[!is.na(gene_list)], decreasing=TRUE)

    withProgress(message="Ejecutando GSEA...", value=0, {
      db     <- input$func_db %||% "kegg"
      result <- NULL

      if (db == "kegg") {
        incProgress(0.3, detail="GSEA con KEGG (requiere internet)...")
        result <- tryCatch(
          clusterProfiler::gseKEGG(
            geneList=gene_list, organism=input$kegg_org %||% "eco",
            minGSSize=input$gsea_min %||% 10, maxGSSize=input$gsea_max %||% 500,
            pvalueCutoff=input$gsea_pval, verbose=FALSE),
          error=function(e){ showNotification(conditionMessage(e), type="error", duration=10); NULL })
      } else {
        ont    <- switch(db, gobp="BP", gomf="MF", gocc="CC", "BP")
        org_pkg <- input$orgdb_pkg %||% ""
        if (!nzchar(org_pkg) || !pkg_ok(org_pkg)) {
          showNotification(paste0("Instala: BiocManager::install('", org_pkg, "')"),
                           type="error", duration=10); return()
        }
        org_db <- tryCatch(getExportedValue(org_pkg, org_pkg), error=function(e){
          suppressPackageStartupMessages(library(org_pkg, character.only=TRUE)); get(org_pkg)
        })
        incProgress(0.3, detail=paste("GSEA GO", ont, "..."))
        result <- tryCatch(
          clusterProfiler::gseGO(
            geneList=gene_list, OrgDb=org_db,
            ont=ont, keyType=input$orgdb_keytype %||% "SYMBOL",
            minGSSize=input$gsea_min %||% 10, maxGSSize=input$gsea_max %||% 500,
            pvalueCutoff=input$gsea_pval, verbose=FALSE),
          error=function(e){ showNotification(conditionMessage(e), type="error", duration=10); NULL })
      }

      incProgress(1.0)
      data_rv$gsea_res <- result
      if (!is.null(result))
        showNotification(sprintf("GSEA: %d conjuntos enriquecidos.", nrow(as.data.frame(result))),
                         type="message")
    })
  })

  output$gsea_bar_plot <- renderPlot({
    req(!is.null(data_rv$gsea_res))
    df <- as.data.frame(data_rv$gsea_res)
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "Sin conjuntos genicos significativos.", cex=1.2)
      return()
    }
    df <- head(df[order(df$NES, decreasing=TRUE), ], 20)
    df$Description <- factor(df$Description, levels=df$Description[order(df$NES)])
    df$Tipo <- ifelse(df$NES > 0, "Activado", "Suprimido")
    ggplot(df, aes(x=NES, y=Description, fill=Tipo)) +
      geom_col(alpha=0.85) +
      scale_fill_manual(values=c(Activado="#D73027", Suprimido="#4575B4")) +
      geom_vline(xintercept=0, linewidth=0.4, color="gray30") +
      xlab("Normalized Enrichment Score (NES)") + ylab(NULL) +
      theme_bw(base_size=11) +
      theme(legend.position="bottom",
            axis.text.y=element_text(size=9))
  })

  output$gsea_table <- renderDT({
    req(!is.null(data_rv$gsea_res))
    df <- as.data.frame(data_rv$gsea_res)
    if (nrow(df) == 0) return(datatable(data.frame(Mensaje="Sin conjuntos significativos.")))
    cols <- intersect(c("ID","Description","setSize","NES","pvalue","p.adjust"),names(df))
    df   <- df[, cols, drop=FALSE]
    for (col in intersect(c("pvalue","p.adjust"), names(df)))
      df[[col]] <- formatC(df[[col]], format="e", digits=3)
    for (col in intersect("NES", names(df)))
      df[[col]] <- round(df[[col]], 3)
    datatable(df, rownames=FALSE, options=list(pageLength=10, scrollX=TRUE, dom="ftip"))
  })

  output$dl_gsea_csv <- downloadHandler(
    filename=function() paste0("GSEA_",format(Sys.time(),"%Y%m%d_%H%M%S"),".csv"),
    content=function(f) {
      req(!is.null(data_rv$gsea_res))
      write.csv(as.data.frame(data_rv$gsea_res), f, row.names=FALSE)
    })

}

shinyApp(ui, server)
