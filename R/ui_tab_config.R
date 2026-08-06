#' ui_tab_config.R
#' Contenido del nav_panel "1 . Configuracion" (6 cards en grid 3x2).

#' Devuelve el contenido (tagList) del Tab 1
ui_tab_config <- function() {
  layout_columns(
    col_widths = c(12),
    div(style = "display:flex; justify-content:flex-end; margin-bottom:12px;",
        actionButton("btn_to_processing",
                     tagList("Continuar al procesamiento", icon("arrow-right")),
                     class = "btn-continue-blue btn-lg")
    ),
    tags$div(
      style = paste(
        "display:grid;",
        "grid-template-columns:repeat(3, minmax(240px, 1fr));",
        "grid-auto-rows:minmax(240px, auto);",
        "gap:12px; align-items:stretch;"
      ),

      # Card 1: Modo de inicio
      card(
        card_header("Modo de inicio"),
        style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between; overflow:auto;",
        radioButtons(
          "start_mode", label = NULL,
          choices  = c("Ejecutar workflow completo" = "workflow",
                       "Análisis a partir de matriz de conteos" = "load"),
          selected = "workflow", inline = TRUE
        ),
        conditionalPanel(
          condition = "input.start_mode === 'load'",
          layout_columns(
            col_widths = c(6, 6),
            div(
              tags$strong("Matriz de conteos (TSV/CSV)"),
              tags$br(),
              tags$small(class = "text-muted",
                "Genes como filas, muestras como columnas. Primera columna = ID de gen."),
              fileInput("upload_counts", label = NULL,
                        accept = c(".tsv", ".csv", ".txt"))
            ),
            div(
              tags$strong("Uso"),
              tags$p(class = "text-muted small",
                     "Carga una matriz externa o usa la pestaña Resultados para revisar ejecuciones previas guardadas automaticamente en outputs/.")
            )
          ),
          div(
            style = "text-align:right; margin-top:8px;",
            actionButton("btn_load_existing",
                         tagList(icon("upload"), " Análisis a partir de matriz de conteos"),
                         class = "btn-success btn-sm")
          )
        )
      ),

      # Card 2: Tipo de analisis
      card(
        card_header("Tipo de analisis"),
        style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between;",
        radioButtons(
          "analysis_type", label = NULL,
          choices  = c("Alineamiento"      = "alignment",
                       "Pseudoalineamiento" = "pseudo"),
          selected = "alignment"
        ),
        radioButtons(
          "read_type", "Tipo de lectura",
          choices = c("Paired-end" = "pe", "Single-end" = "se"),
          selected = "pe",
          inline = TRUE
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
                      selected = "salmon"),
          conditionalPanel(
            "input.pseudo_tool === 'kallisto' && input.read_type === 'se'",
            layout_columns(
              col_widths = c(6, 6),
              numericInput("fragment_length", "Longitud media fragmento", value = 200, min = 1, step = 1),
              numericInput("fragment_sd", "SD fragmento", value = 20, min = 1, step = 1)
            )
          )
        )
      ),

      # Card 3: Rutas y archivos
      card(
        card_header("Rutas y archivos"),
        style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between; overflow:auto;",
        tags$div(style = "display:flex;gap:12px;align-items:flex-start;",
          tags$div(style = "flex:1;",
            shinyDirButton("input_dir_btn", "Seleccionar directorio de FASTQs", "Seleccionar..."),
            tags$div(style = "margin-top:6px;", textOutput("input_dir_path")),
            uiOutput("genome_label_ui"),
            conditionalPanel(
              "input.analysis_type === 'alignment'",
              fileInput("annotation_file_upload", "Archivo de anotacion GFF/GTF (requerido)",
                        accept = c(".gff", ".gtf", ".gff3", ".gz"), multiple = FALSE)
            ),
            conditionalPanel(
              "input.analysis_type === 'pseudo'",
              tagList(
                fileInput("annotation_file_pseudo_upload", "Archivo de anotacion GFF/GTF (opcional)",
                          accept = c(".gff", ".gtf", ".gff3", ".gz"), multiple = FALSE),
                tags$small(class = "text-muted", icon("circle-info"),
                           " Si se omite, se pasa /dev/null al script.")
              )
            ),
            tags$strong("Directorio de salida automatico"),
            tags$div(style = "margin-top:6px;", textOutput("output_base_path"))
          )
        )
      ),

      # Card 4: Resumen
      card(
        card_header("Resumen"),
        style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between;",
        tags$dl(class = "row small mb-0",
          tags$dt(class = "col-6", "Entrada:"),
          tags$dd(class = "col-6", tags$code(class = "small", textOutput("input_dir_path", inline = TRUE))),
          tags$dt(class = "col-6", "Salida (base):"),
          tags$dd(class = "col-6", tags$code(textOutput("output_base_path", inline = TRUE))),
          tags$dt(class = "col-6", "Workflow:"),
          tags$dd(class = "col-6", tags$code(textOutput("workflow_path_text", inline = TRUE)))
        )
      ),

      # Card 5: Checklist
      card(
        card_header("Checklist"),
        style = "min-height:240px; display:flex; flex-direction:column; justify-content:center; overflow:auto;",
        tags$div(style = "flex:1; display:flex; align-items:center; justify-content:center;",
          uiOutput("checklist_ui")
        )
      ),

      # Card 6: Muestras detectadas
      card(
        card_header("Muestras detectadas"),
        style = "min-height:240px; overflow:auto; max-height:420px;",
        uiOutput("sample_preview_ui")
      ),

      # Card 7: Potencia a priori. Va en configuracion porque la decision que
      # informa (cuantas replicas) se toma ANTES de secuenciar, no despues.
      if (isTRUE(HAS_RNASEQPOWER)) card(
        card_header("Potencia y tamano muestral"),
        style = "grid-column: 1 / -1;",
        tags$p(class = "small text-muted mb-2",
               paste("El tamano muestral es el determinante mas fuerte de la",
                     "calidad del resultado. Esta calculadora orienta, no",
                     "garantiza: ninguna herramienta es fiable cuando se exigen",
                     "efectos pequenos y potencias altas, porque los parametros",
                     "no se pueden fijar bien desde datos piloto limitados.")),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          numericInput("pw_n", "Replicas por grupo", value = 3, min = 2, max = 100, step = 1),
          numericInput("pw_effect", "Fold-change a detectar", value = 2, min = 1.1,
                       max = 10, step = 0.1),
          numericInput("pw_cv", "CV biologico", value = 0.4, min = 0.05, max = 1.5,
                       step = 0.05),
          numericInput("pw_depth", "Profundidad media por gen", value = 20, min = 1,
                       max = 1000, step = 5)
        ),
        tags$small(class = "text-muted d-block mb-2",
                   paste("CV tipico: 0,1 en lineas celulares, 0,4 en muestras",
                         "humanas. La profundidad es el numero medio de conteos por",
                         "gen, no el total de lecturas.")),
        uiOutput("pw_verdict"),
        plotly::plotlyOutput("pw_curve", height = "300px")
      ) else NULL
    )
  )
}
