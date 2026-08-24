#' ui_tab_results.R
#' Helper que devuelve el tagList para el renderUI("tab3_content") (Tab 3).
#'
#' Reorganizacion respecto de la versión anterior:
#'
#'   - Había DOS navset_tab apilados en la misma pantalla. El primero
#'     (Resumen / Calidad / Conteos / Informes / Log) y, debajo, una sección
#'     titulada "Resultados adicionales de alineamiento y pseudoalineamiento"
#'     con su propia barra de pestanas SIEMPRE visible, con independencia de la
#'     pestana elegida arriba. Dos barras de pestanas en una pagina, una de las
#'     cuales no responde a la otra, es la principal razón de que esta vista
#'     resultara confusa. Ahora hay un único navset y el control de calidad del
#'     método es una pestana más.
#'   - Ese bloque mostraba además las metricas de alineamiento Y las de
#'     pseudoalineamiento a la vez. Una ejecución solo puede ser de un tipo, así
#'     que la mitad de las tarjetas decía siempre "No se encontraron metricas
#'     para este análisis": diez gráficos vacios como estado normal. Ahora se
#'     muestra la pestana que corresponde a la herramienta de la ejecución.
#'   - `download_log` y `download_qc_alerts` estaban definidos en el server pero
#'     no tenían ningun boton en la interfaz: eran descargas inalcanzables. Se
#'     colocan en la pestana Log y en el resumen, respectivamente.

#' Metricas de cabecera. El color del filo lo fija el estado de la metrica, no
#' su posición en la fila.
ui_results_tiles <- function(s) {
  mapped <- if (!is.null(s)) s$mean_mapped else NA_real_
  nivel_mapeo <- if (is.null(s) || is.na(mapped)) "neutral"
                 else if (mapped < 50) "bad"
                 else if (mapped < 70) "warn"
                 else "ok"
  estado <- if (!is.null(s)) s$status else NULL
  nivel_estado <- switch(estado %||% "",
                         completado = "ok", error = "bad",
                         incompleto = "warn", "neutral")

  tags$div(class = "stat-grid",
    stat_tile(if (!is.null(s)) status_badge(s$status) else "—", "Estado", nivel_estado),
    stat_tile(if (!is.null(s)) paste0((s$n_samples %||% "—"), " / ", (s$n_features %||% "—")) else "—",
              "Muestras / genes", if (is.null(s)) "neutral" else "info"),
    stat_tile(if (!is.null(s)) pct_label(mapped) else "—", "Mapeo medio", nivel_mapeo),
    # Antes esta tarjeta repetia n_features, que es el número de FILAS de la
    # matriz y no los genes detectados. Mostrar dos veces el mismo número con
    # etiquetas distintas induce a error: ahora son los genes con al menos una
    # lectura asignada.
    stat_tile(if (!is.null(s)) (s$n_detected %||% "—") else "—",
              "Genes con conteo > 0", if (is.null(s)) "neutral" else "info"),
    # Coste de la ejecución. Es lo que permite responder a "esto cabe en mi
    # portatil o necesito un servidor", que es una pregunta que se hace antes
    # de lanzar un conjunto de datos grande, no después.
    stat_tile(if (!is.null(s)) fmt_duracion(s$duration_seconds) else "—",
              "Duración", if (is.null(s)) "neutral" else "info"),
    stat_tile(if (!is.null(s)) fmt_memoria(s$peak_rss_mb) else "—",
              "RAM máxima", if (is.null(s)) "neutral" else "info")
  )
}

#' Pestana de control de calidad del alineamiento clasico.
ui_results_align_qc <- function() {
  tagList(
    layout_columns(
      col_widths = c(6, 6),
      card(download_header("Resumen de alineamiento", "download_align_qc_summary"),
           DTOutput("align_qc_summary_table")),
      card(download_header("Alertas interpretativas", "download_align_qc_alerts"),
           DTOutput("align_qc_alerts_table"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(download_header("Lecturas únicas, multimapeadas y no alineadas", "download_align_qc_mapping_plot"),
           plotly::plotlyOutput("align_qc_mapping_plot", height = "320px")),
      card(download_header("Lecturas asignadas y no asignadas", "download_align_qc_assignment_plot"),
           plotly::plotlyOutput("align_qc_assignment_plot", height = "320px"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(download_header("Distribución exonica, intronica e intergenica", "download_align_qc_region_plot"),
           plotly::plotlyOutput("align_qc_region_plot", height = "280px")),
      card(download_header("Cobertura 5'-3' del cuerpo genico", "download_align_qc_gene_body_plot"),
           plotly::plotlyOutput("align_qc_gene_body_plot", height = "280px"))
    )
  )
}

#' Pestana de control de calidad del pseudoalineamiento.
ui_results_pseudo_qc <- function() {
  tagList(
    layout_columns(
      col_widths = c(6, 6),
      card(download_header("Resumen de pseudoalineamiento", "download_pseudo_qc_summary"),
           DTOutput("pseudo_qc_summary_table")),
      card(download_header("Alertas interpretativas", "download_pseudo_qc_alerts"),
           DTOutput("pseudo_qc_alerts_table"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(download_header("Pseudoalignment rate por muestra", "download_pseudo_qc_rate_plot"),
           plotly::plotlyOutput("pseudo_qc_rate_plot", height = "300px")),
      card(download_header("Distribución de TPM por muestra", "download_pseudo_qc_tpm_plot"),
           plotly::plotlyOutput("pseudo_qc_tpm_plot", height = "300px"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(download_header("Transcritos detectados por muestra", "download_pseudo_qc_detected_plot"),
           plotly::plotlyOutput("pseudo_qc_detected_plot", height = "300px")),
      card(download_header("TPM vs NumReads", "download_pseudo_qc_scatter_plot"),
           plotly::plotlyOutput("pseudo_qc_scatter_plot", height = "300px"))
    ),
    card(download_header("Cuantificación", "download_pseudo_quant"),
         DTOutput("pseudo_qc_quant_table"))
  )
}

#' Contenido principal del Tab 3 cuando hay datos o ejecuciones disponibles.
#'
#' @param s Resumen de la ejecución seleccionada (`summarise_result`), o NULL.
#' @param p Parámetros inferidos de la ejecución. Solo se usa `p$tool`, para
#'   decidir que pestana de control de calidad del método tiene sentido mostrar.
ui_tab_results_content <- function(s, p = list()) {
  tool      <- tolower(p$tool %||% "")
  es_pseudo <- tool %in% c("salmon", "kallisto")
  es_align  <- tool %in% c("bowtie2", "alignment")
  # Herramienta desconocida (por ejemplo, una matriz subida): se ofrecen las dos
  # y que decida quien mira, en lugar de esconder la que si tiene datos.
  if (!es_pseudo && !es_align) { es_pseudo <- TRUE; es_align <- TRUE }

  paneles_qc <- c(
    if (es_align)  list(nav_panel("QC alineamiento", ui_results_align_qc())),
    if (es_pseudo) list(nav_panel("QC pseudoalineamiento", ui_results_pseudo_qc()))
  )

  tagList(
    # El selector va primero y a ancho completo: es el control que gobierna todo
    # lo que hay debajo, y estaba metido en una tarjeta de media anchura a la
    # derecha de las metricas que el mismo determina.
    card(
      card_header("Ejecución analizada"),
      layout_columns(
        col_widths = c(8, 4),
        uiOutput("result_selector_ui"),
        uiOutput("multiqc_open_ui")
      )
    ),

    ui_results_tiles(s),

    tags$div(class = "mt-3",
      do.call(navset_tab, c(
        list(id = "results_tabs"),
        list(nav_panel(
          "Resumen",
          layout_columns(
            col_widths = c(5, 7),
            card(
              card_header("Interpretación rápida"),
              uiOutput("result_interpretation_ui"),
              tags$hr(class = "my-2"),
              # Reune las alertas de los dos bloques de QC en un único CSV. El
              # handler existía desde el principio; lo que faltaba era el boton.
              downloadButton("download_qc_alerts", "Alertas de calidad (CSV)",
                             icon = icon("download"),
                             class = "btn-sm btn-outline-secondary")
            ),
            card(
              download_header("Estadísticas principales MultiQC", "download_run_stats"),
              DTOutput("run_stats_table")
            )
          )
        )),
        list(nav_panel(
          "Calidad",
          layout_columns(
            col_widths = c(5, 7),
            card(
              download_header("Estado FastQC", "download_fastqc"),
              DTOutput("fastqc_table")
            ),
            card(
              download_header("Trimming y cuantificación/alineamiento", "download_alignment"),
              DTOutput("alignment_table")
            )
          ),
          # Semaforo por muestra con los modulos de FastQC, que aplican los cortes
          # del manual (aviso si el cuartil inferior baja de 10 o la mediana de 25;
          # fallo si bajan de 5 y 20). Se indica QUE modulo falla, porque un fallo
          # de calidad por base y uno de adaptadores no se arreglan igual.
          card(
            card_header("Calidad por muestra (FastQC)"),
            uiOutput("fastqc_light_note"),
            DTOutput("fastqc_light_table")
          ),
          layout_columns(
            col_widths = c(6, 6),
            # rRNA como metrica de primer nivel: en procariotas supera el 85 % del
            # RNA celular y la depleccion es imperfecta, así que un arrastre
            # desigual entre muestras sesga los factores de tamaño de TODA la
            # normalización.
            card(
              download_header("Lecturas asignadas a rRNA por muestra", "download_rrna_plot"),
              uiOutput("rrna_note"),
              plotly::plotlyOutput("rrna_plot", height = "320px")
            ),
            # Saturación: si al 50 % de la profundidad ya se detectan casi los
            # mismos genes, secuenciar más no aporta y el límite son las replicas.
            card(
              download_header("Saturación de la libreria", "download_saturation_plot"),
              uiOutput("saturation_note"),
              plotly::plotlyOutput("saturation_plot", height = "320px")
            )
          )
        )),
        list(nav_panel(
          "Conteos",
          layout_columns(
            col_widths = c(5, 7),
            card(
              download_header("Librerias por muestra", "download_count_lib"),
              DTOutput("count_lib_table")
            ),
            card(
              download_header("Genes/transcritos más abundantes", "download_count_top"),
              DTOutput("count_top_table")
            )
          )
        )),
        paneles_qc,
        list(nav_panel(
          "Informes y archivos",
          card(
            download_header("Informes principales", "download_artifacts"),
            DTOutput("artifact_table")
          ),
          card(
            download_header("Todos los archivos generados", "download_filelist"),
            DTOutput("output_files_table")
          )
        )),
        list(nav_panel(
          "Log",
          card(
            download_header("Últimas líneas del workflow", "download_log"),
            class = "log-pre",
            verbatimTextOutput("selected_log_tail")
          )
        ))
      ))
    )
  )
}

#' Aviso mostrado cuando no hay resultados (Tab 3 sin run y sin outputs guardados).
#' Incluye boton para refrescar la lista de outputs/.
ui_tab_results_empty <- function() {
  div(class = "alert alert-info mt-4", icon("clock"),
      " Ejecuta el workflow en la pestana 2 para ver los resultados o usa una ejecución previa guardada en outputs/.",
      div(style = "margin-top:10px;",
          actionButton("refresh_results_btn",
                       tagList(icon("rotate"), " Refrescar lista de outputs/"),
                       class = "btn-outline-secondary btn-sm")))
}
