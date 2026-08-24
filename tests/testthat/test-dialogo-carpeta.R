# Dialogo nativo para elegir carpeta (R/utils_dialogos.R).
#
# No se abre ninguna ventana: se comprueba que la ruta devuelta se interpreta
# bien y, sobre todo, que CANCELAR no se confunde con un fallo. Es el caso más
# frecuente después de elegir, y tratarlo como error llenaria la interfaz de
# avisos por pulsar el boton de cerrar.

test_that("una ruta elegida se devuelve limpia y normalizada", {
  d <- withr::local_tempdir()
  # macOS devuelve las carpetas con barra final y un salto de línea.
  falso <- function(...) structure(paste0(d, "/"), status = 0L)
  expect_equal(elegir_directorio_nativo(ejecutar = falso),
               normalizePath(d, mustWork = FALSE))
})

test_that("cancelar devuelve NULL y no un error", {
  # osascript sale con estado 1 y el texto "User canceled." cuando se cierra.
  cancelado <- function(...) structure("User canceled. (-128)", status = 1L)
  expect_null(elegir_directorio_nativo(ejecutar = cancelado))

  # zenity y kdialog simplemente no escriben nada.
  vacio <- function(...) structure(character(0), status = 0L)
  expect_null(elegir_directorio_nativo(ejecutar = vacio))
  expect_null(elegir_directorio_nativo(ejecutar = function(...) structure("", status = 0L)))
})

test_that("una ruta que ya no existe no se da por buena", {
  # Entre que se elige y se lee puede haberse borrado, y devolverla haría que el
  # pipeline fallase mucho más adelante con un error mucho menos claro.
  fantasma <- function(...) structure("/no/existe/esta/carpeta", status = 0L)
  expect_null(elegir_directorio_nativo(ejecutar = fantasma))
})

test_that("un fallo del propio dialogo se traga sin romper la aplicación", {
  expect_null(elegir_directorio_nativo(ejecutar = function(...) stop("sin escritorio")))
  expect_null(elegir_directorio_nativo(ejecutar = function(...) warning("nada")))
})

test_that("el punto de partida cae al home cuando no sirve", {
  d <- withr::local_tempdir()
  recibidos <- NULL
  espia <- function(cmd, args, ...) { recibidos <<- args; structure(paste0(d, "/"), status = 0L) }

  elegir_directorio_nativo(inicio = "/no/existe", ejecutar = espia)
  expect_true(any(grepl(path.expand("~"), recibidos, fixed = TRUE)))

  elegir_directorio_nativo(inicio = "", ejecutar = espia)
  expect_true(any(grepl(path.expand("~"), recibidos, fixed = TRUE)))

  # Y cuando si sirve, se respeta.
  elegir_directorio_nativo(inicio = d, ejecutar = espia)
  expect_true(any(grepl(d, recibidos, fixed = TRUE)))
})

# ── Ruta escrita o pegada a mano ────────────────────────────────────────────

test_that("una carpeta válida se acepta y se normaliza", {
  d <- withr::local_tempdir()
  v <- validar_directorio(d)
  expect_null(v$error)
  expect_equal(v$ruta, normalizePath(d, mustWork = FALSE))

  # `~` es lo que la gente teclea.
  v2 <- validar_directorio("~")
  expect_null(v2$error)
  expect_true(nzchar(v2$ruta))
})

test_that("se quitan las comillas que deja arrastrar una carpeta al terminal", {
  d <- withr::local_tempdir()
  for (t in c(sprintf('"%s"', d), sprintf("'%s'", d), sprintf("  %s  ", d))) {
    v <- validar_directorio(t)
    expect_null(v$error, info = t)
    expect_equal(v$ruta, normalizePath(d, mustWork = FALSE), info = t)
  }
})

test_that("el campo vacio no es un error, pero una ruta mala si lo explica", {
  # Vacio es el estado inicial: marcarlo en rojo antes de que el usuario escriba
  # nada sería un reproche por no haber empezado.
  v <- validar_directorio("")
  expect_null(v$error); expect_equal(v$ruta, "")
  expect_null(validar_directorio(NULL)$error)

  expect_match(validar_directorio("/no/existe/nada")$error, "no existe")

  # Un fichero no es una carpeta, y decirlo así ahorra el "pero si existe".
  f <- withr::local_tempfile(); file.create(f)
  v2 <- validar_directorio(f)
  expect_equal(v2$ruta, "")
  expect_match(v2$error, "fichero, no una carpeta")
})

test_that("dentro de un contenedor no se ofrece el dialogo", {
  # No hay escritorio: ofrecer el boton solo llevaría a una aplicación
  # congelada esperando una ventana que nunca va a aparecer. Se cae entonces al
  # selector de shinyFiles, que si funciona con el servidor en otra maquina.
  expect_false(dialogo_nativo_disponible(contenedor = TRUE, grafico = TRUE))
})

test_that("en un escritorio de verdad si se ofrece", {
  # En macOS y Windows siempre hay forma; en Linux depende de zenity o kdialog,
  # así que ahí se comprueba solo que la respuesta es coherente consigo misma.
  disponible <- dialogo_nativo_disponible(contenedor = FALSE, grafico = TRUE)
  expect_type(disponible, "logical")
  if (Sys.info()[["sysname"]] == "Darwin") {
    expect_true(disponible)
  }
})

test_that("sin entorno gráfico en Linux no se ofrece", {
  skip_if(Sys.info()[["sysname"]] == "Darwin" || .Platform$OS.type == "windows",
          "macOS y Windows siempre tienen escritorio cuando hay sesión")
  expect_false(dialogo_nativo_disponible(contenedor = FALSE, grafico = FALSE))
})
