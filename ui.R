#' ui.R
#' UI principal: navbar con la portada mas las cuatro pestanas del analisis.
#' El contenido detallado de cada pestana vive en los helpers ui_tab_*. Todas
#' salvo la de configuracion se renderizan dinamicamente (uiOutput) porque
#' dependen del estado del server.

ui <- page_navbar(
  # La marca en la barra: icono, nombre y expansion del acronimo. El icono va
  # como <img> y no en linea porque asi el navegador lo cachea una vez y no se
  # reenvia el SVG entero en cada render de la UI.
  title  = tags$span(
    class = "app-brand",
    tags$img(src = "sara_icono.svg", alt = "", height = "26", class = "app-logo"),
    tags$span(class = "app-name", "SARA"),
    tags$span(class = "app-subtitle", "Shiny App for RNA-seq Analysis")
  ),
  id     = "main_nav",
  theme  = app_theme,
  header = tagList(
    useShinyjs(),
    tags$head(
      tags$link(rel = "icon", type = "image/png", sizes = "32x32",
                href = "sara_icono_32.png"),
      tags$link(rel = "icon", type = "image/png", sizes = "64x64",
                href = "sara_icono_64.png"),
      tags$link(rel = "apple-touch-icon", href = "sara_icono_180.png")
    ),
    tags$style(app_css)
  ),
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
