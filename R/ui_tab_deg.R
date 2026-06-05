#' ui_tab_deg.R
#' Tab 4: Expresion diferencial (DEG).
#' Tres motores (DESeq2 / edgeR / limma-voom), metadatos editables y filtros
#' interactivos. Se renderiza via uiOutput("tab_deg_content") desde server_tab_deg.

#' Aviso cuando faltan paquetes Bioconductor necesarios.
ui_tab_deg_missing_pkgs <- function() {
  pkgs <- c()
  if (!isTRUE(HAS_DESEQ2))  pkgs <- c(pkgs, "DESeq2")
  if (!isTRUE(HAS_EDGER))   pkgs <- c(pkgs, "edgeR")
  if (!isTRUE(HAS_LIMMA))   pkgs <- c(pkgs, "limma")
  if (!length(pkgs)) return(NULL)
  div(class = "alert alert-warning mt-3",
      icon("triangle-exclamation"),
      tags$b(" Paquetes ausentes: "),
      paste(pkgs, collapse = ", "),
      tags$p(class = "mb-0 small",
             "Instala los paquetes indicados para activar el motor correspondiente."))
}

#' UI del Tab 4 (cards de configuracion + zona de resultados).
ui_tab_deg <- function() {
  tagList(
    ui_tab_deg_missing_pkgs(),
    tags$div(
      style = paste(
        "display:grid;",
        "grid-template-columns:repeat(3, minmax(280px, 1fr));",
        "grid-auto-rows:minmax(220px, auto);",
        "gap:12px; align-items:stretch; margin-top:6px;"
      ),

      # Card 1: Datos
      card(
        card_header("1. Datos"),
        radioButtons("deg_source", label = NULL,
                     choices = c("Ejecucion actual" = "current",
                                 "Ejecucion guardada" = "saved",
                                 "Matriz subida" = "upload"),
                     selected = "current"),
        conditionalPanel(
          "input.deg_source === 'saved'",
          uiOutput("deg_saved_run_ui")
        ),
        conditionalPanel(
          "input.deg_source === 'upload'",
          fileInput("deg_counts_upload", "Matriz de conteos (TSV/CSV)",
                    accept = c(".tsv", ".csv", ".txt"))
        ),
        tags$small(class = "text-muted",
                   "Genes por filas, muestras por columnas. Se autodetectan los IDs de muestra.")
      ),

      # Card 2: Metadatos
      card(
        card_header("2. Metadatos"),
        fileInput("deg_meta_upload", "Samplesheet (CSV/TSV)",
                  accept = c(".csv", ".tsv", ".txt")),
        div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-bottom:6px;",
            actionButton("deg_meta_sync_btn",
                         tagList(icon("rotate"), " Sincronizar con muestras detectadas"),
                         class = "btn-sm"),
            actionButton("deg_meta_add_row_btn",
                         tagList(icon("plus"), " Anadir fila"),
                         class = "btn-sm")
        ),
        tags$small(class = "text-muted",
                   "Columnas requeridas: sample_id, condition. Opcional: batch."),
        DTOutput("deg_meta_table")
      ),

      # Card 3: Diseno
      card(
        card_header("3. Diseno experimental"),
        selectInput("deg_condition_col", "Columna de condicion",
                    choices = c("condition"), selected = "condition"),
        selectInput("deg_ref_level", "Nivel de referencia",
                    choices = NULL),
        checkboxInput("deg_use_batch", "Incluir variable batch", FALSE),
        conditionalPanel(
          "input.deg_use_batch === true",
          selectInput("deg_batch_col", "Columna de batch",
                      choices = c("batch"), selected = "batch")
        )
      ),

      # Card 4: Analisis
      card(
        card_header("4. Analisis DEG"),
        selectInput("deg_method", "Motor",
                    choices = c("DESeq2", "edgeR", "limma-voom"),
                    selected = "DESeq2"),
        layout_columns(
          col_widths = c(6, 6),
          numericInput("deg_min_count", "Min. cuenta por fila", value = 10, min = 0, step = 1),
          numericInput("deg_min_samples", "Min. muestras", value = 2, min = 1, step = 1)
        ),
        actionButton("run_deg_btn",
                     tagList(icon("flask"), " Lanzar DEG"),
                     class = "btn-success btn-lg"),
        tags$div(style = "margin-top:8px;",
                 verbatimTextOutput("deg_status_text", placeholder = TRUE))
      ),

      # Card 5: Filtros interactivos
      card(
        card_header("5. Filtros (recalculan al vuelo)"),
        sliderInput("deg_fdr_cutoff", "FDR maximo",
                    min = 0, max = 1, value = 0.05, step = 0.01),
        sliderInput("deg_log2fc_cutoff", "|log2FC| minimo",
                    min = 0, max = 5, value = 1, step = 0.1),
        sliderInput("deg_basemean_cutoff", "baseMean minimo",
                    min = 0, max = 1000, value = 0, step = 5)
      )
    ),

    # Zona de resultados (cargada al lanzar run_deg)
    tags$div(class = "mt-3", uiOutput("deg_results_ui"))
  )
}

#' tagList con el navset de resultados (tabla, plots, enriquecimiento).
ui_tab_deg_results <- function() {
  tagList(
    navset_tab(
      id = "deg_result_tabs",
      nav_panel(
        "Tabla",
        card(
          download_header("Tabla DEG (filtrada)", "download_deg_table"),
          DTOutput("deg_table")
        )
      ),
      nav_panel(
        "Volcano",
        card(
          download_header("Volcano plot", "download_deg_volcano_plot"),
          plotly::plotlyOutput("deg_volcano_plot", height = "480px")
        )
      ),
      nav_panel(
        "MA",
        card(
          download_header("MA plot", "download_deg_ma_plot"),
          plotly::plotlyOutput("deg_ma_plot", height = "480px")
        )
      ),
      nav_panel(
        "PCA",
        card(
          download_header("PCA (vst/rlog)", "download_deg_pca_plot"),
          plotly::plotlyOutput("deg_pca_plot", height = "480px")
        )
      ),
      nav_panel(
        "Heatmap top-N",
        card(
          card_header(tags$div(
            class = "card-title-download",
            tags$span("Heatmap top-N genes"),
            div(style = "display:flex;gap:6px;align-items:center;",
              numericInput("deg_heatmap_topn", NULL, value = 30, min = 5, max = 200, step = 5, width = "100px"),
              downloadButton("download_deg_heatmap", label = NULL, icon = icon("download"),
                             class = "btn-sm btn-outline-secondary header-download",
                             title = "Descargar heatmap como PNG")
            )
          )),
          plotOutput("deg_heatmap", height = "560px")
        )
      ),
      nav_panel(
        "Distancia entre muestras",
        card(
          download_header("Matriz de distancias entre muestras", "download_deg_dist_heatmap"),
          plotOutput("deg_dist_heatmap", height = "520px")
        )
      ),
      nav_panel(
        "Enriquecimiento",
        card(
          card_header(tags$div(
            class = "card-title-download",
            tags$span("Enriquecimiento funcional"),
            div(style = "display:flex;gap:6px;align-items:flex-end;flex-wrap:wrap;",
                selectInput("deg_ontology", NULL,
                            choices = c("GO: Procesos biologicos" = "BP",
                                        "GO: Funcion molecular" = "MF",
                                        "GO: Componente celular" = "CC",
                                        "KEGG" = "KEGG"),
                            selected = "BP", width = "200px"),
                conditionalPanel(
                  "input.deg_ontology === 'KEGG'",
                  textInput("deg_kegg_organism", NULL, value = "eco",
                            placeholder = "eco, hsa, mmu...", width = "110px")
                ),
                actionButton("deg_run_enrich_btn",
                             tagList(icon("play"), " Calcular"),
                             class = "btn-sm"))
          )),
          plotly::plotlyOutput("deg_enrich_dotplot", height = "440px"),
          tags$hr(),
          tags$div(
            class = "card-title-download",
            style = "margin-bottom:6px;",
            tags$span("Tabla de terminos enriquecidos"),
            downloadButton("download_enrich_table", label = NULL, icon = icon("download"),
                           class = "btn-sm btn-outline-secondary header-download",
                           title = "Descargar tabla de enriquecimiento")
          ),
          DTOutput("deg_enrich_table")
        )
      )
    )
  )
}
