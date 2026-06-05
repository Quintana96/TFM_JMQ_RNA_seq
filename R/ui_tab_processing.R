#' ui_tab_processing.R
#' Helper que devuelve el tagList para el renderUI("tab2_content") (Tab 2).

#' Construye el contenido del Tab 2 a partir del snapshot de configuracion `cfg`.
#' Si process_unlocked es FALSE, muestra un aviso con un boton de retorno a Tab 1.
ui_tab_processing_content <- function(cfg, total_sz) {
  tagList(
    # Fila superior: resumen + progreso
    layout_columns(
      col_widths = c(4, 8),

      # Izquierda: resumen y controles
      card(
        card_header("Resumen del analisis"),
        tags$dl(class = "row small mb-1",
          tags$dt(class = "col-6", "Tipo:"),
          tags$dd(class = "col-6", if (cfg$analysis_type == "alignment")
            "Alineamiento (Bowtie2)" else paste0("Pseudoalineamiento (", cfg$tool, ")")),
          tags$dt(class = "col-6", "Muestras:"), tags$dd(class = "col-6", cfg$n_samples),
          tags$dt(class = "col-6", "Lectura:"),  tags$dd(class = "col-6", cfg$read_type),
          if (identical(cfg$tool, "kallisto") && identical(cfg$read_type, "Single-end"))
            tagList(
              tags$dt(class = "col-6", "Fragmento:"),
              tags$dd(class = "col-6", paste0(cfg$fragment_length, " ± ", cfg$fragment_sd))
            ),
          tags$dt(class = "col-6", "Tamaño est.:"), tags$dd(class = "col-6", total_sz),
          tags$dt(class = "col-6", "Entrada:"),
          tags$dd(class = "col-6", tags$code(class = "small", cfg$input_dir)),
          tags$dt(class = "col-6", "Salida:"),
          tags$dd(class = "col-6", tags$code(class = "small", cfg$output_dir))
        ),
        hr(),
        tags$details(
          tags$summary(tags$small(icon("terminal"), " Ver comando")),
          tags$pre(class = "small mt-1",
                   style = "white-space:pre-wrap;word-break:break-all;background:#f8f9fa;padding:8px;",
                   textOutput("cmd_preview_text", inline = TRUE))
        ),
        hr(),
        div(style = "display:flex;gap:8px;align-items:center;",
            actionButton("btn_back", tagList(icon("arrow-left"), " Volver"), class = "btn-secondary"),
            actionButton("stop_btn", tagList(icon("stop"), " Detener"), class = "btn-danger", disabled = "disabled"),
            actionButton("run_btn",  tagList(icon("play"), " Ejecutar workflow"), class = "btn-success btn-lg")
        )
      ),

      # Derecha: progreso detallado
      card(
        card_header("Progreso en tiempo real"),

        # Metricas de tiempo y muestras
        div(class = "d-flex gap-3 mb-3 p-2 rounded",
            style = "background:#F7F1FC; border:1px solid #E3D6F2; font-size:.9rem;",
            div(icon("clock"), " Transcurrido: ",
                tags$b(textOutput("elapsed_text", inline = TRUE))),
            div(icon("hourglass-half"), " Tiempo restante: ",
                tags$b(textOutput("eta_text", inline = TRUE))),
            div(icon("percent"), " % completado: ",
                tags$b(textOutput("remaining_pct", inline = TRUE))),
            div(icon("vials"), " Muestras: ",
                tags$b(textOutput("samp_prog_text", inline = TRUE)))
        ),

        uiOutput("checkpoints_ui"),
        hr(),
        uiOutput("sample_status_ui")
      )
    ),

    # Fila inferior: log ancho completo
    card(
      card_header(
        "Log de terminal",
        div(class = "float-end",
            actionButton("refresh_btn", icon("rotate"),
                         class = "btn-sm btn-outline-secondary",
                         title = "Actualizar lista de archivos"))
      ),
      tags$small(class = "text-muted",
                 "Se muestra stdout/stderr en vivo. Copia completa: workflow_live.log en el directorio de salida."),
      verbatimTextOutput("run_log"),
      max_height = "420px",
      style = "overflow-y:auto; width:100%; background:#ffffff; color:#111827;"
    )
  )
}

#' Aviso mostrado cuando process_unlocked = FALSE (con boton para volver a Tab 1).
ui_tab_processing_locked <- function() {
  div(class = "alert alert-info mt-4",
      icon("arrow-left-long"),
      " Completa la pestaña 1 antes de continuar.",
      div(style = "margin-top:10px;",
          actionButton("btn_goto_config",
                       tagList(icon("gear"), " Ir a configuracion"),
                       class = "btn-secondary btn-sm")))
}
