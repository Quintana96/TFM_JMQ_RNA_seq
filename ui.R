#' ui.R
#' UI principal: navbar de tres pestanas. El contenido detallado de cada
#' pestana vive en los helpers ui_tab_*. Tab 2 y Tab 3 se renderizan
#' dinamicamente (uiOutput) porque dependen del estado del server.

ui <- page_navbar(
  title  = tags$span("RNA-seq Workflow Runner"),
  id     = "main_nav",
  theme  = app_theme,
  header = tagList(useShinyjs(), tags$style(app_css)),

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
