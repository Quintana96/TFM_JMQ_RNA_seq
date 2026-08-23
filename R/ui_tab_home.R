#' ui_tab_home.R
#' Tab 0: portada / panel de la aplicacion.
#'
#' La aplicacion no tenia punto de entrada: se abria directamente sobre la
#' pestana de configuracion, que ademas mezclaba la configuracion propiamente
#' dicha con el estado de la sesion y con una calculadora de potencia. Quien la
#' abria por primera vez no tenia forma de saber en cuantos pasos consistia el
#' analisis ni por donde iba.
#'
#' Esta pestana solo hace dos cosas: enumerar los cuatro pasos con su estado
#' actual y resumir la sesion. No contiene ningun control del analisis, a
#' proposito: es un indice, no un formulario.

#' Tarjeta-paso de la portada.
#'
#' Es un `<button class="action-button">`, no un div con un boton dentro: Shiny
#' vincula cualquier elemento con la clase `action-button` como si fuera un
#' actionButton, de modo que toda la tarjeta es el area pulsable. Obligar a
#' acertar un boton pequeno al pie de la tarjeta era parte de lo que hacia
#' torpe la navegacion.
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
    # El titulo no repite el nombre de la aplicacion, que ya esta en la barra de
    # navegacion a cuarenta pixeles de aqui: dice lo que la aplicacion hace.
    tags$div(class = "home-hero",
      tags$h1("Del FASTQ al gen diferencial"),
      tags$p(paste(
        "Pipeline completo de RNA-seq: control de calidad, cuantificacion,",
        "expresion diferencial con varios motores y enriquecimiento funcional.",
        "Elige un paso para empezar; cada tarjeta indica en que estado esta."
      ))
    ),

    tags$div(class = "home-grid",
      home_step_card("home_go_config", "1", "Configuracion",
                     paste("Modo de inicio, tipo de analisis, FASTQ, genoma y",
                           "anotacion. Valida los parametros antes de ejecutar."),
                     steps$config$pill, steps$config$cta),
      home_step_card("home_go_process", "2", "Procesamiento",
                     paste("Ejecuta workflow.sh y sigue el progreso en vivo:",
                           "checkpoints, muestras y log de terminal."),
                     steps$process$pill, steps$process$cta),
      home_step_card("home_go_results", "3", "Resultados",
                     paste("Calidad, alineamiento y conteos de cualquier",
                           "ejecucion guardada en outputs/."),
                     steps$results$pill, steps$results$cta),
      home_step_card("home_go_deg", "4", "Expresion diferencial",
                     paste("DESeq2, edgeR y limma-voom, con",
                           "diagnosticos, enriquecimiento e informe."),
                     steps$deg$pill, steps$deg$cta)
    ),

    if (!is.null(hint)) tags$div(class = "mt-3", hint) else NULL,

    section_title("Estado de la sesion",
                  "Lo que la aplicacion tiene cargado ahora mismo."),
    tags$div(class = "stat-grid", tiles)
  )
}
