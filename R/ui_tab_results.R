#' ui_tab_results.R
#' Helper que devuelve el tagList para el renderUI("tab3_content") (Tab 3).

#' Contenido principal del Tab 3 cuando hay datos o ejecuciones disponibles
ui_tab_results_content <- function(s) {
  tagList(
    layout_columns(
      col_widths = c(6, 6),
      # Izquierda: cajas 2x2 con metricas resumen
      div(style = "display:grid; grid-template-columns:repeat(2,1fr); gap:10px;",
        div(class = "metric-card",
          div(class = "metric-card-value",
            if (!is.null(s)) status_badge(s$status) else "—"),
          div(class = "metric-card-label", "Estado")
        ),
        div(class = "metric-card",
          div(class = "metric-card-value",
            if (!is.null(s)) paste0((s$n_samples %||% "—"), " / ", (s$n_features %||% "—")) else "—"),
          div(class = "metric-card-label", "Muestras / genes")
        ),
        div(class = "metric-card",
          div(class = "metric-card-value",
            if (!is.null(s)) pct_label(s$mean_mapped) else "—"),
          div(class = "metric-card-label", "Mapeo medio")
        ),
        div(class = "metric-card",
          div(class = "metric-card-value",
            if (!is.null(s)) (s$n_features %||% "—") else "—"),
          div(class = "metric-card-label", "Genes detectados")
        )
      ),
      # Derecha: selector vertical de resultados
      card(
        card_header("Resultados"),
        div(style = "min-height:160px; display:flex; flex-direction:column; justify-content:center; gap:10px;",
            uiOutput("result_selector_ui"),
            uiOutput("multiqc_open_ui")
        )
      )
    ),
    navset_tab(
      id = "results_tabs",
      nav_panel(
        "Resumen",
        layout_columns(
          col_widths = c(5, 7),
          card(
            card_header("Interpretacion rapida"),
            uiOutput("result_interpretation_ui")
          ),
          card(
            download_header("Estadisticas principales MultiQC", "download_run_stats"),
            DTOutput("run_stats_table")
          )
        )
      ),
      nav_panel(
        "Calidad",
        layout_columns(
          col_widths = c(5, 7),
          card(
            download_header("Estado FastQC", "download_fastqc"),
            DTOutput("fastqc_table")
          ),
          card(
            download_header("Trimming y cuantificacion/alineamiento", "download_alignment"),
            DTOutput("alignment_table")
          )
        ),
        # rRNA como metrica de primer nivel: en procariotas supera el 85 % del RNA
        # celular y la depleccion es imperfecta, asi que un arrastre desigual
        # entre muestras sesga los factores de tamano de TODA la normalizacion.
        card(
          download_header("Lecturas asignadas a rRNA por muestra", "download_rrna_plot"),
          uiOutput("rrna_note"),
          plotly::plotlyOutput("rrna_plot", height = "320px")
        )
      ),
      nav_panel(
        "Conteos",
        layout_columns(
          col_widths = c(5, 7),
          card(
            download_header("Librerias por muestra", "download_count_lib"),
            DTOutput("count_lib_table")
          ),
          card(
            download_header("Genes/transcritos mas abundantes", "download_count_top"),
            DTOutput("count_top_table")
          )
        )
      ),
      nav_panel(
        "Informes y archivos",
        card(
          download_header("Informes principales", "download_artifacts"),
          DTOutput("artifact_table")
        ),
        card(
          download_header("Todos los archivos generados", "download_filelist"),
          DTOutput("output_files_table")
        )
      ),
      nav_panel(
        "Log",
        card(
          card_header("Ultimas lineas del workflow"),
          verbatimTextOutput("selected_log_tail")
        )
      )
    ),

    # Resultados adicionales de alineamiento y pseudoalineamiento
    tags$div(
      class = "mt-4",
      tags$h3("Resultados adicionales de alineamiento y pseudoalineamiento"),
      navset_tab(
        id = "method_qc_tabs",
        nav_panel(
          "Alineamiento",
          tags$h4("Control de calidad del alineamiento"),
          layout_columns(
            col_widths = c(6, 6),
            card(download_header("Resumen de alineamiento", "download_align_qc_summary"),
                 DTOutput("align_qc_summary_table")),
            card(download_header("Alertas interpretativas", "download_align_qc_alerts"),
                 DTOutput("align_qc_alerts_table"))
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(download_header("Lecturas unicas, multimapeadas y no alineadas", "download_align_qc_mapping_plot"),
                 plotly::plotlyOutput("align_qc_mapping_plot", height = "320px")),
            card(download_header("Lecturas asignadas y no asignadas", "download_align_qc_assignment_plot"),
                 plotly::plotlyOutput("align_qc_assignment_plot", height = "320px"))
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(download_header("Distribucion exonica, intronica e intergenica", "download_align_qc_region_plot"),
                 plotly::plotlyOutput("align_qc_region_plot", height = "280px")),
            card(download_header("Cobertura 5'-3' del cuerpo genico", "download_align_qc_gene_body_plot"),
                 plotly::plotlyOutput("align_qc_gene_body_plot", height = "280px"))
          )
        ),
        nav_panel(
          "Pseudoalineamiento",
          tags$h4("Control de calidad del pseudoalineamiento"),
          layout_columns(
            col_widths = c(12),
            card(download_header("Resumen de pseudoalineamiento", "download_pseudo_qc_summary"),
                 DTOutput("pseudo_qc_summary_table"))
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(download_header("Pseudoalignment rate por muestra", "download_pseudo_qc_rate_plot"),
                 plotly::plotlyOutput("pseudo_qc_rate_plot", height = "300px")),
            card(download_header("Distribucion de TPM por muestra", "download_pseudo_qc_tpm_plot"),
                 plotly::plotlyOutput("pseudo_qc_tpm_plot", height = "300px"))
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(download_header("Transcritos detectados por muestra", "download_pseudo_qc_detected_plot"),
                 plotly::plotlyOutput("pseudo_qc_detected_plot", height = "300px")),
            card(download_header("TPM vs NumReads", "download_pseudo_qc_scatter_plot"),
                 plotly::plotlyOutput("pseudo_qc_scatter_plot", height = "300px"))
          ),
          layout_columns(
            col_widths = c(12),
            card(download_header("Cuantificacion", "download_pseudo_quant"),
                 DTOutput("pseudo_qc_quant_table"))
          ),
          card(download_header("Alertas interpretativas", "download_pseudo_qc_alerts"),
               DTOutput("pseudo_qc_alerts_table"))
        )
      )
    )
  )
}

#' Aviso mostrado cuando no hay resultados (Tab 3 sin run y sin outputs guardados).
#' Incluye boton para refrescar la lista de outputs/.
ui_tab_results_empty <- function() {
  div(class = "alert alert-info mt-4", icon("clock"),
      " Ejecuta el workflow en la pestaña 2 para ver los resultados o usa una ejecucion previa guardada en outputs/.",
      div(style = "margin-top:10px;",
          actionButton("refresh_results_btn",
                       tagList(icon("rotate"), " Refrescar lista de outputs/"),
                       class = "btn-secondary btn-sm")))
}
