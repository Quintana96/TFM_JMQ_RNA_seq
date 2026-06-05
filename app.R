#' app.R
#' Entrypoint de la aplicacion. Cuando se ejecuta con `Rscript app.R`,
#' cargamos explicitamente los archivos que Shiny normalmente espera en
#' una app multi-archivo.

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
app_file <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else NULL
app_dir <- if (!is.null(app_file)) {
  dirname(normalizePath(app_file, mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(app_dir, "global.R"), local = TRUE)

for (r_file in sort(list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(r_file, local = TRUE)
}

source(file.path(app_dir, "ui.R"), local = TRUE)
source(file.path(app_dir, "server.R"), local = TRUE)

app <- shiny::shinyApp(ui, server)

is_direct_rscript <- length(file_arg) > 0L

if (is_direct_rscript) {
  shiny::runApp(app, launch.browser = interactive())
} else {
  app
}
