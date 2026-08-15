#' ui.R
#' UI principal: navbar con la portada mas las cuatro pestanas del analisis.
#' El contenido detallado de cada pestana vive en los helpers ui_tab_*. Todas
#' salvo la de configuracion se renderizan dinamicamente (uiOutput) porque
#' dependen del estado del server.

ui <- page_navbar(
  title  = tags$span("RNA-seq Workflow Runner"),
  id     = "main_nav",
  theme  = app_theme,
  header = tagList(useShinyjs(), tags$style(app_css)),
  fillable = FALSE,

  # Tab 0 — Portada. Es un indice de los cuatro pasos con su estado, no un
  # formulario: la aplicacion no tenia punto de entrada y se abria directamente
  # sobre un tablero de seis tarjetas de configuracion.
  nav_panel(
    title = tagList(icon("house"), " Inicio"),
    value = "tab_home",
    uiOutput("home_content")
  ),

  # Tab 1 — Configuracion (estatica)
  nav_panel(
    title = "1 · Configuracion",
    value = "tab_config",
    ui_tab_config()
  ),

  # Tab 2 — Procesamiento (dinamica via renderUI)
  nav_panel(
    title = "2 · Procesamiento",
    value = "tab_process",
    uiOutput("tab2_content")
  ),

  # Tab 3 — Resultados (dinamica via renderUI)
  nav_panel(
    title = "3 · Resultados",
    value = "tab_results",
    uiOutput("tab3_content")
  ),

  # Tab 4 — Expresion diferencial (dinamica via renderUI)
  nav_panel(
    title = "4 · Expresion diferencial",
    value = "tab_deg",
    uiOutput("tab_deg_content")
  )
)
