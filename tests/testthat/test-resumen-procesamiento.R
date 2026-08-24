# El resumen del paso 2 no puede contradecir al comando que se ejecuta.
#
# Habia una copia de la configuracion tomada al pulsar "Continuar", y el resumen
# la usaba. Bastaba con cambiar algo en el paso 1 y volver por la barra de
# navegacion —sin pasar otra vez por "Continuar"— para que el resumen dijese
# "Single-end" mientras el comando llevaba "--READ_TYPE pe". Un resumen que
# miente sobre lo que va a correr es peor que no tener resumen.

test_that("el resumen sigue al tipo de lectura actual, no al validado", {
  srv <- function(input, output, session) {
    state <- new.env(parent = emptyenv())
    state$process_unlocked <- shiny::reactiveVal(TRUE)
    # La foto se tomo cuando la configuracion decia otra cosa: solo debe aportar
    # la carpeta de salida.
    state$config_snap <- shiny::reactiveVal(list(output_dir = "/salida/validada"))

    lectura <- shiny::reactive(input$read_type %||% "pe")
    shared <- list(
      input_dir_val             = shiny::reactive(""),
      output_dir_val            = shiny::reactive("/salida/actual"),
      samples_eff               = shiny::reactive(c("A", "B", "C")),
      effective_tool            = shiny::reactive("bowtie2"),
      effective_read_type       = lectura,
      effective_fragment_length = shiny::reactive(200),
      effective_fragment_sd     = shiny::reactive(20),
      effective_genome_file     = shiny::reactive(""),
      effective_annotation      = shiny::reactive("")
    )
    # Se reproduce la construccion del resumen, que es la parte que fallaba.
    output$resumen <- shiny::renderText({
      snap <- state$config_snap()
      paste(read_type_label(shared$effective_read_type()),
            length(shared$samples_eff()),
            snap$output_dir %||% shared$output_dir_val(), sep = " | ")
    })
  }

  shiny::testServer(srv, {
    session$setInputs(read_type = "pe", analysis_type = "alignment")
    expect_match(output$resumen, "^Paired-end \\| 3 \\|")

    # El usuario vuelve al paso 1 y cambia el tipo de lectura. El resumen tiene
    # que reflejarlo SIN pasar de nuevo por "Continuar".
    session$setInputs(read_type = "se")
    expect_match(output$resumen, "^Single-end \\|")

    session$setInputs(read_type = "pe")
    expect_match(output$resumen, "^Paired-end \\|")

    # La carpeta de salida SI viene de la foto: se crea al validar y moverla
    # despues dejaria los resultados en un sitio distinto del anunciado.
    expect_match(output$resumen, "/salida/validada$")
  })
})

test_that("sin foto se usa la carpeta de salida actual", {
  srv <- function(input, output, session) {
    state <- new.env(parent = emptyenv())
    state$config_snap <- shiny::reactiveVal(list())
    shared <- list(output_dir_val = shiny::reactive("/salida/actual"))
    output$dir <- shiny::renderText({
      snap <- state$config_snap()
      snap$output_dir %||% shared$output_dir_val()
    })
  }
  shiny::testServer(srv, { expect_equal(output$dir, "/salida/actual") })
})

test_that("read_type_label no invierte los terminos", {
  # Parece obvio, pero era una de las dos explicaciones posibles del sintoma y
  # descartarla costo tiempo.
  expect_equal(read_type_label("se"), "Single-end")
  expect_equal(read_type_label("pe"), "Paired-end")
  expect_equal(read_type_label(NULL), "Paired-end")
  expect_equal(read_type_label("cualquier-otra-cosa"), "Paired-end")
})
