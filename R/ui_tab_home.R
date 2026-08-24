#' ui_tab_home.R
#' Tab 0: portada / panel de la aplicación.
#'
#' La aplicación no tenía punto de entrada: se abria directamente sobre la
#' pestana de configuración, que además mezclaba la configuración propiamente
#' dicha con el estado de la sesión y con una calculadora de potencia. Quien la
#' abria por primera vez no tenía forma de saber en cuantos pasos consistia el
#' análisis ni por donde iba.
#'
#' Esta pestana solo hace dos cosas: enumerar los cuatro pasos con su estado
#' actual y resumir la sesión. No contiene ningun control del análisis, a
#' propósito: es un índice, no un formulario.

#' Tarjeta-paso de la portada.
#'
#' Es un `<button class="action-button">`, no un div con un boton dentro: Shiny
#' vincula cualquier elemento con la clase `action-button` como si fuera un
#' actionButton, de modo que toda la tarjeta es el área pulsable. Obligar a
#' acertar un boton pequeño al pie de la tarjeta era parte de lo que hacía
#' torpe la navegación.
home_step_card <- function(id, step, title, desc, pill, cta = "Abrir") {
  tags$button(
    id = id,
    class = "action-button home-card",
    type = "button",
    tags$div(class = "home-card-top",
      tags$div(style = "display:flex;align-items:center;gap:.55rem;",
        tags$span(class = "home-card-step", step),
        tags$span(class = "home-card-title", title)
      ),
      pill
    ),
    tags$p(class = "home-card-desc", desc),
    tags$span(class = "home-card-cta", cta, HTML("&nbsp;→"))
  )
}

#' Contenido completo de la portada.
#'
#' @param steps Lista de cuatro listas con `pill` (tag) y `cta` (texto).
#' @param tiles Lista de `stat_tile()` ya construidos.
#' @param hint  Aviso contextual opcional (tag) bajo las tarjetas.
ui_tab_home_content <- function(steps, tiles, hint = NULL) {
  tagList(
    # El título no repite el nombre de la aplicación, que ya está en la barra de
    # navegación a cuarenta pixeles de aquí: dice lo que la aplicación hace.
    tags$div(class = "home-hero",
      tags$h1("Del FASTQ al gen diferencial"),
      tags$p(paste(
        "Pipeline completo de RNA-seq: control de calidad, cuantificación,",
        "expresión diferencial con varios motores y enriquecimiento funcional.",
        "Elige un paso para empezar; cada tarjeta indica en que estado está."
      ))
    ),

    tags$div(class = "home-grid",
      home_step_card("home_go_config", "1", "Configuración",
                     paste("Modo de inicio, tipo de análisis, FASTQ, genoma y",
                           "anotación. Válida los parámetros antes de ejecutar."),
                     steps$config$pill, steps$config$cta),
      home_step_card("home_go_process", "2", "Procesamiento",
                     paste("Ejecuta workflow.sh y sigue el progreso en vivo:",
                           "checkpoints, muestras y log de terminal."),
                     steps$process$pill, steps$process$cta),
      home_step_card("home_go_results", "3", "Resultados",
                     paste("Calidad, alineamiento y conteos de cualquier",
                           "ejecución guardada en outputs/."),
                     steps$results$pill, steps$results$cta),
      home_step_card("home_go_deg", "4", "Expresión diferencial",
                     paste("DESeq2, edgeR y limma-voom, con",
                           "diagnósticos, enriquecimiento e informe."),
                     steps$deg$pill, steps$deg$cta)
    ),

    if (!is.null(hint)) tags$div(class = "mt-3", hint) else NULL,

    section_title("Estado de la sesión",
                  "Lo que la aplicación tiene cargado ahora mismo."),
    tags$div(class = "stat-grid", tiles)
  )
}
