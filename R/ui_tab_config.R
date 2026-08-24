#' ui_tab_config.R
#' Contenido del nav_panel "1 . Configuración".
#'
#' Reorganizacion respecto de la versión anterior (un grid de seis a ocho
#' tarjetas de igual peso):
#'
#'   - Se separa lo que el usuario RELLENA (columna izquierda) de lo que la
#'     aplicación le RESPONDE (columna derecha: checklist, errores, muestras).
#'     Antes ambas cosas se alternaban en el mismo grid, de modo que el
#'     checklist podia quedar a dos tarjetas de distancia del campo que lo
#'     dejaba en rojo.
#'   - La tarjeta "Resumen" se ha eliminado: repetia `textOutput("input_dir_path")`
#'     y `textOutput("output_base_path")`, que ya estaban en "Rutas y archivos".
#'     Shiny vincula cada valor de salida a UN solo elemento del DOM, así que la
#'     segunda copia se quedaba permanentemente en blanco; en la pantalla se veia
#'     un "Directorio de salida automático" vacio con la ruta correcta a dos
#'     tarjetas de distancia. Cada salida aparece ahora una única vez.
#'   - `uiOutput("validation_ui")` estaba definido en el server pero no existía
#'     en ninguna parte de la interfaz: la lista de errores concretos no se
#'     mostraba nunca y el boton de continuar se quedaba gris sin explicar por
#'     que. Ahora acompaña al checklist.
#'   - Las opciones avanzadas y la calculadora de potencia bajan a un acordeon
#'     plegado. Son útiles, pero no forman parte de configurar una ejecución, y
#'     ocupaban dos terceras partes del alto de la pestana.

#' Bloque de rutas efectivas (cada salida aparece una sola vez en toda la app).
ui_config_paths <- function() {
  tags$dl(class = "row small mb-0 mt-2",
    tags$dt(class = "col-4 text-muted", "Entrada"),
    tags$dd(class = "col-8 mb-1", tags$code(textOutput("input_dir_path", inline = TRUE))),
    tags$dt(class = "col-4 text-muted", "Salida"),
    tags$dd(class = "col-8 mb-1", tags$code(textOutput("output_base_path", inline = TRUE))),
    tags$dt(class = "col-4 text-muted", "Workflow"),
    tags$dd(class = "col-8 mb-0", tags$code(textOutput("workflow_path_text", inline = TRUE)))
  )
}

#' Acordeon con lo que no forma parte de configurar la ejecución en curso.
ui_config_extras <- function() {
  paneles <- list(
    accordion_panel(
      "Opciones avanzadas del pipeline",
      icon = icon("sliders"),
      value = "adv",
      tags$p(class = "small text-muted",
             "Orientación de la libreria, hilos y replicas inferenciales."),
      layout_columns(
        col_widths = c(4, 4, 4),
        selectInput("adv_strandedness", "Orientación de la libreria",
                    choices = c("Inferir de los datos" = "auto",
                                "Sin orientar (-s 0)" = "0",
                                "Directa (-s 1)" = "1",
                                "Inversa, protocolos dUTP (-s 2)" = "2"),
                    selected = "auto"),
        numericInput("adv_threads", "Hilos",
                     value = max(1, parallel::detectCores() - 1),
                     min = 1, max = 64, step = 1),
        numericInput("adv_inferential_reps", "Replicas inferenciales",
                     value = 20, min = 0, max = 100, step = 5)
      ),
      tags$small(class = "text-muted d-block",
                 paste("Contar como no orientada una libreria stranded suma las",
                       "lecturas antisentido y degrada la especificidad, sobre",
                       "todo en genomas densos como los procariotas. Las replicas",
                       "inferenciales (salmon/kallisto) cuantifican la incertidumbre",
                       "de asignación de lecturas; 0 las desactiva."))
    )
  )

  # La calculadora de potencia va en configuración porque la decisión que
  # informa (cuantas replicas) se toma ANTES de secuenciar, no después.
  if (isTRUE(HAS_RNASEQPOWER)) {
    paneles <- c(paneles, list(accordion_panel(
      "Potencia y tamaño muestral",
      icon = icon("chart-line"),
      value = "power",
      tags$p(class = "small text-muted mb-2",
             paste("El tamaño muestral es el determinante más fuerte de la",
                   "calidad del resultado. Esta calculadora orienta, no",
                   "garantiza: ninguna herramienta es fiable cuando se exigen",
                   "efectos pequeños y potencias altas, porque los parámetros",
                   "no se pueden fijar bien desde datos piloto limitados.")),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        numericInput("pw_n", "Replicas por grupo", value = 3, min = 2, max = 100, step = 1),
        numericInput("pw_effect", "Fold-change a detectar", value = 2, min = 1.1,
                     max = 10, step = 0.1),
        numericInput("pw_cv", "CV biológico", value = 0.4, min = 0.05, max = 1.5,
                     step = 0.05),
        numericInput("pw_depth", "Profundidad media por gen", value = 20, min = 1,
                     max = 1000, step = 5)
      ),
      # El CV y la profundidad son justo los parámetros que quien usa la app no
      # conoce, y de los que depende todo el resultado. Si hay una matriz
      # cargada se pueden MEDIR en lugar de adivinarse.
      tags$div(class = "d-flex align-items-center gap-2 mb-2 flex-wrap",
        actionButton("pw_estimate_btn",
                     tagList(icon("wand-magic-sparkles"), " Estimar de mis datos"),
                     class = "btn-outline-secondary btn-sm"),
        tags$small(class = "text-muted", uiOutput("pw_estimate_note", inline = TRUE))
      ),
      tags$small(class = "text-muted d-block mb-2",
                 paste("CV típico: 0,1 en líneas celulares, 0,4 en muestras",
                       "humanas. La profundidad es el número medio de conteos por",
                       "gen, no el total de lecturas.")),
      uiOutput("pw_verdict"),
      plotly::plotlyOutput("pw_curve", height = "300px")
    )))
  }

  tagList(
    section_title("Herramientas de apoyo",
                  "No hacen falta para lanzar el análisis; se despliegan cuando se necesitan."),
    do.call(accordion, c(paneles, list(open = FALSE, multiple = TRUE)))
  )
}

#' Devuelve el contenido (tagList) del Tab 1
ui_tab_config <- function() {
  tagList(

    # ── Modo de inicio: gobierna todo lo que se ve debajo ────────────────────
    card(
      card_header("Modo de inicio"),
      radioButtons(
        "start_mode", label = NULL,
        choices  = c("Ejecutar workflow completo" = "workflow",
                     "Análisis a partir de matriz de conteos" = "load"),
        selected = "workflow", inline = TRUE
      ),
      # En modo "matriz" no hay pipeline que configurar. Antes se seguian
      # mostrando las tarjetas de FASTQ, genoma y anotación junto a un checklist
      # en rojo que exigia rellenarlas: requisitos que en ese modo no aplican.
      conditionalPanel(
        condition = "input.start_mode === 'load'",
        layout_columns(
          col_widths = c(7, 5),
          div(
            tags$strong("Matriz de conteos (TSV/CSV)"),
            tags$br(),
            tags$small(class = "text-muted",
              "Genes como filas, muestras como columnas. Primera columna = ID de gen."),
            fileInput("upload_counts", label = NULL,
                      accept = c(".tsv", ".csv", ".txt"), width = "100%")
          ),
          div(
            tags$strong("Uso"),
            tags$p(class = "text-muted small",
                   paste("Carga una matriz externa o usa la pestana Resultados",
                         "para revisar ejecuciones previas guardadas",
                         "automaticamente en outputs/.")),
            actionButton("btn_load_existing",
                         tagList(icon("upload"), " Cargar y continuar"),
                         class = "btn-primary")
          )
        )
      )
    ),

    # ── Configuración del pipeline (solo en modo workflow) ──────────────────
    conditionalPanel(
      condition = "input.start_mode === 'workflow'",

      layout_columns(
        col_widths = c(7, 5),

        # Columna izquierda: lo que el usuario rellena
        tags$div(
          card(
            card_header("Tipo de análisis"),
            layout_columns(
              col_widths = c(6, 6),
              radioButtons(
                "analysis_type", "Estrategia",
                choices  = c("Alineamiento"       = "alignment",
                             "Pseudoalineamiento" = "pseudo"),
                selected = "alignment"
              ),
              radioButtons(
                "read_type", "Tipo de lectura",
                choices = c("Paired-end" = "pe", "Single-end" = "se"),
                selected = "pe"
              )
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
                  numericInput("fragment_length", "Longitud media fragmento",
                               value = 200, min = 1, step = 1),
                  numericInput("fragment_sd", "SD fragmento", value = 20, min = 1, step = 1)
                )
              )
            )
          ),

          card(
            card_header("Rutas y archivos"),
            # Dos caminos hacía la misma ruta. El boton abre el dialogo del
            # sistema —Finder, Explorador o el del escritorio— que es donde la
            # gente sabe moverse; el campo de texto permite pegarla directamente,
            # que es más rápido cuando ya se conoce. El campo muestra además la
            # selección actual, así que no hace falta una etiqueta aparte.
            tags$label(class = "control-label", "Directorio de FASTQs"),
            div(class = "d-flex gap-2 align-items-start mb-1",
                uiOutput("input_dir_btn_ui", inline = TRUE),
                div(class = "flex-grow-1",
                    textInput("input_dir_texto", label = NULL, value = "",
                              placeholder = "/ruta/a/los/fastq  (o pulsa Examinar)",
                              width = "100%"))),
            uiOutput("input_dir_aviso"),
            uiOutput("genome_label_ui"),
            conditionalPanel(
              "input.analysis_type === 'alignment'",
              fileInput("annotation_file_upload",
                        "Archivo de anotación GFF/GTF (requerido)",
                        accept = c(".gff", ".gtf", ".gff3", ".gz"), multiple = FALSE,
                        width = "100%"),
              # Aviso si la anotación tiene genes con varios exones: bowtie2 no es
              # splice-aware y perdería las uniones exon-exon.
              uiOutput("splice_warning_ui")
            ),
            conditionalPanel(
              "input.analysis_type === 'pseudo'",
              fileInput("annotation_file_pseudo_upload",
                        "Archivo de anotación GFF/GTF (opcional)",
                        accept = c(".gff", ".gtf", ".gff3", ".gz"), multiple = FALSE,
                        width = "100%"),
              tags$small(class = "text-muted", icon("circle-info"),
                         " Si se omite, se pasa /dev/null al script.")
            ),
            tags$hr(class = "my-2"),
            ui_config_paths()
          )
        ),

        # Columna derecha: lo que la aplicación responde
        tags$div(
          card(
            card_header("Estado de la configuración"),
            uiOutput("checklist_ui"),
            tags$hr(class = "my-2"),
            uiOutput("validation_ui"),
            tags$div(class = "mt-3 d-grid",
              actionButton("btn_to_processing",
                           tagList("Continuar al procesamiento", icon("arrow-right")),
                           class = "btn-primary btn-lg")
            )
          ),
          card(
            card_header("Muestras detectadas"),
            max_height = "460px",
            uiOutput("sample_preview_ui")
          )
        )
      )
    ),

    # Fuera del conditionalPanel: la calculadora de potencia sirve igual en modo
    # matriz, donde además es cuando puede estimar el CV y la profundidad de los
    # datos reales en lugar de pedirlos a ojo.
    ui_config_extras()
  )
}
