#' test-arranque-cableado.R
#' Dos guardas sobre la estructura de la app, no sobre su estadística:
#'
#'   1. La raíz de la aplicación. `runApp()` solo fija el working directory
#'      cuando recibe una RUTA; app.R le pasa un objeto app, así que con
#'      `Rscript app.R` desde otro directorio el wd seguía siendo el del usuario.
#'      Como `outputs_base_dir()` es `file.path(getwd(), "outputs")` y la ruta de
#'      workflow.sh se resuelve igual, la app creaba outputs/ en el sitio
#'      equivocado y no encontraba el workflow.
#'   2. El cableado interfaz-servidor. Un `observeEvent(input$x, ...)` cuyo `x`
#'      no lo declara ningun control es código muerto silencioso: no falla, no
#'      avisa, y la funcionalidad que anuncia su comentario simplemente no
#'      existe. Era el caso de `input$btn_home`.

# ── Raíz de la aplicación ───────────────────────────────────────────────────

#' Extrae `resolve_app_dir()` de app.R sin ejecutar el resto del fichero (que
#' arrancaria el servidor). Se evalua solo la expresión que la define.
cargar_resolve_app_dir <- function() {
  exprs <- parse(file.path(app_root, "app.R"))
  env <- new.env(parent = globalenv())
  encontrada <- FALSE
  for (e in exprs) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1]])[1] %in% c("<-", "=") &&
        identical(as.character(e[[2]])[1], "resolve_app_dir")) {
      eval(e, envir = env)
      encontrada <- TRUE
    }
  }
  if (!encontrada) return(NULL)
  get("resolve_app_dir", envir = env)
}

test_that("app.R expone resolve_app_dir() y fija el wd con su resultado", {
  resolver <- cargar_resolve_app_dir()
  expect_false(is.null(resolver))

  # El resultado tiene que USARSE: sin el setwd, resolver bien la ruta no
  # arregla nada, porque las rutas de salida se resuelven contra getwd().
  fuente <- paste(readLines(file.path(app_root, "app.R"), warn = FALSE),
                  collapse = "\n")
  codigo <- grep("^\\s*#", strsplit(fuente, "\n", fixed = TRUE)[[1]],
                 value = TRUE, invert = TRUE)
  expect_true(any(grepl("setwd(app_dir)", codigo, fixed = TRUE)))
})

test_that("con `Rscript app.R` desde otro directorio se resuelve la raíz real", {
  resolver <- cargar_resolve_app_dir()
  skip_if(is.null(resolver), "resolve_app_dir() no está en app.R")

  ajeno <- withr::local_tempdir()   # un wd sin global.R: el caso que fallaba
  obtenido <- resolver(paste0("--file=", file.path(app_root, "app.R")), wd = ajeno)

  expect_equal(normalizePath(obtenido, mustWork = FALSE),
               normalizePath(app_root, mustWork = FALSE))
  # Y con esa raíz, las rutas derivadas caen dentro del proyecto y no del wd.
  expect_true(file.exists(file.path(obtenido, "global.R")))
  expect_false(startsWith(normalizePath(obtenido, mustWork = FALSE),
                          normalizePath(ajeno, mustWork = FALSE)))
})

test_that("cuando Shiny ya ha fijado el wd, manda el wd y no --file", {
  resolver <- cargar_resolve_app_dir()
  skip_if(is.null(resolver), "resolve_app_dir() no está en app.R")

  # `Rscript lanzador.R` que llama a runApp("<app>"): Shiny fija el wd a la app,
  # pero --file apunta al lanzador. Fiarse de --file sourcearia global.R del
  # sitio equivocado.
  lanzador <- withr::local_tempdir()
  obtenido <- resolver(paste0("--file=", file.path(lanzador, "lanzador.R")),
                       wd = app_root)

  expect_equal(normalizePath(obtenido, mustWork = FALSE),
               normalizePath(app_root, mustWork = FALSE))
})

test_that("sin ninguna pista utilizable se cae al working directory", {
  resolver <- cargar_resolve_app_dir()
  skip_if(is.null(resolver), "resolve_app_dir() no está en app.R")

  ajeno <- withr::local_tempdir()
  expect_equal(normalizePath(resolver(character(0), wd = ajeno), mustWork = FALSE),
               normalizePath(ajeno, mustWork = FALSE))
})

# ── Cableado interfaz-servidor ──────────────────────────────────────────────

#' Ids de input que la aplicación lee pero que no crea ningun control de la
#' interfaz. La lista blanca recoge los que Shiny o DT sintetizan solos.
inputs_huerfanos <- function() {
  ficheros <- c(list.files(file.path(app_root, "R"), pattern = "[.]R$",
                           full.names = TRUE),
                file.path(app_root, c("ui.R", "server.R", "global.R")))
  lineas <- unlist(lapply(ficheros[file.exists(ficheros)], readLines, warn = FALSE))
  lineas <- lineas[!grepl("^\\s*#", lineas)]
  todo <- paste(lineas, collapse = "\n")

  leidos <- unique(gsub(".*input\\$([A-Za-z0-9_.]+).*", "\\1",
                        regmatches(todo, gregexpr("input\\$[A-Za-z0-9_.]+", todo))[[1]]))

  controles <- c("selectInput", "selectizeInput", "checkboxInput", "numericInput",
                 "sliderInput", "actionButton", "actionLink", "textInput",
                 "textAreaInput", "radioButtons", "fileInput", "dateInput",
                 "checkboxGroupInput", "shinyDirButton", "shinyFilesButton")
  patron <- paste0("(", paste(controles, collapse = "|"), ")\\(\\s*\"[A-Za-z0-9_.]+\"")
  declarados <- unique(gsub(".*\"([A-Za-z0-9_.]+)\"", "\\1",
                            regmatches(todo, gregexpr(patron, todo))[[1]]))

  # Los genera el cliente, no un control declarado en la interfaz:
  #   deg_meta_table_cell_edit  lo emite DT al editar una celda (id de la tabla
  #                             + "_cell_edit").
  sinteticos <- c("deg_meta_table_cell_edit")

  setdiff(leidos, c(declarados, sinteticos))
}

test_that("no hay observers escuchando inputs que ningun control declara", {
  # `input$btn_home` era exactamente esto: un observer con su comentario
  # ("enlace de vuelta a la portada desde cualquier pestana") para un boton que
  # no existía en ninguna parte de la interfaz.
  expect_equal(inputs_huerfanos(), character(0))
})
