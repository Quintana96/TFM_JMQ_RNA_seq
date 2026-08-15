#' server.R
#' Server principal: instancia el estado compartido y delega en las
#' funciones server_tab_*. No usa moduleServer() / NS() para preservar
#' los IDs Shiny originales (input$run_btn, output$tab2_content, etc.).

server <- function(input, output, session) {
  state <- create_app_state(session)

  server_tab_config(input, output, session, state)
  server_tab_processing(input, output, session, state)
  server_tab_results(input, output, session, state)
  server_tab_deg(input, output, session, state)
  # La portada va la ultima: lee `state$shared`, que server_tab_config rellena
  # al final de su ejecucion.
  server_tab_home(input, output, session, state)
}
