#' app.R
#' Entrypoint de la aplicacion. Cuando se ejecuta con `Rscript app.R`,
#' cargamos explicitamente los archivos que Shiny normalmente espera en
#' una app multi-archivo.

#' Localiza la raiz de la aplicacion a partir del working directory y del
#' `--file=` de la linea de comandos.
#'
#' Ninguna de las dos pistas vale por si sola:
#'   - Shiny fija el wd al directorio de la app antes de sourcear app.R, asi que
#'     cuando ahi esta global.R esa es la respuesta correcta AUNQUE `--file`
#'     apunte a otra cosa (p. ej. `Rscript lanzador.R` que llama a runApp()).
#'   - Con `Rscript app.R` desde otro directorio pasa lo contrario: el wd es el
#'     del usuario y la pista buena es `--file`.
#' Se usa global.R como fichero ancla, igual que hace el helper de los tests.
resolve_app_dir <- function(file_arg, wd = getwd()) {
  tiene_app <- function(d) !is.null(d) && nzchar(d) && file.exists(file.path(d, "global.R"))
  wd <- normalizePath(wd, mustWork = FALSE)
  desde_file <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else NULL
  if (tiene_app(wd)) return(wd)
  if (tiene_app(desde_file)) return(desde_file)
  wd
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
app_dir <- resolve_app_dir(file_arg)

# `runApp()` solo fija el working directory cuando recibe una RUTA; aqui recibe
# un objeto app, asi que con `Rscript app.R` desde otro directorio el wd seguia
# siendo el del usuario. Como `outputs_base_dir()` es `file.path(getwd(),
# "outputs")` y la ruta de workflow.sh se resuelve igual, la app creaba outputs/
# en el sitio equivocado y no encontraba el workflow. No se restaura al salir a
# proposito: el wd tiene que seguir siendo este mientras la app viva.
setwd(app_dir)

source(file.path(app_dir, "global.R"), local = TRUE)

for (r_file in sort(list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(r_file, local = TRUE)
}

source(file.path(app_dir, "ui.R"), local = TRUE)
source(file.path(app_dir, "server.R"), local = TRUE)

app <- shiny::shinyApp(ui, server)

is_direct_rscript <- length(file_arg) > 0L

if (is_direct_rscript) {
  # Host y puerto configurables por entorno. Fuera de un contenedor no cambia
  # nada (se mantienen los valores por defecto de Shiny); dentro son
  # imprescindibles: con 127.0.0.1 el servidor arranca y queda inalcanzable
  # desde fuera del contenedor, y con un puerto aleatorio no hay nada que
  # publicar. Un puerto mal escrito se ignora en lugar de tumbar el arranque.
  # SARA_PORT y no SHINY_PORT: Shiny reserva esa variable para deducir que se
  # ejecuta bajo Shiny Server, y con ella definida avisa en cada arranque de que
  # "Shiny Server v0.3.4 or later is required", que no viene a cuento y asusta.
  # Se siguen leyendo las antiguas como respaldo para no romper despliegues ya
  # hechos con ellas.
  primera_no_vacia <- function(...) {
    for (v in c(...)) if (nzchar(v)) return(v)
    ""
  }
  host <- primera_no_vacia(Sys.getenv("SARA_HOST", ""), Sys.getenv("SHINY_HOST", ""))
  port <- suppressWarnings(as.integer(
    primera_no_vacia(Sys.getenv("SARA_PORT", ""), Sys.getenv("SHINY_PORT", ""))))
  if (nzchar(host)) options(shiny.host = host)
  if (!is.na(port) && port > 0) options(shiny.port = port)
  # Con `Rscript app.R` la sesion no es interactiva, de modo que hasta ahora no
  # se abria nada y habia que ir al navegador a mano. `lanzador_de_interfaz()`
  # decide que hacer segun la maquina: una ventana de aplicacion sin adornos de
  # navegador en un escritorio, nada dentro de un contenedor. Se puede forzar
  # con SARA_UI=app|navegador|ninguno. Ver R/utils_launch.R.
  shiny::runApp(app, launch.browser = lanzador_de_interfaz())
} else {
  app
}
