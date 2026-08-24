# Apertura de la interfaz en una ventana sin adornos de navegador
# (R/utils_launch.R).
#
# Nada de esto lanza un navegador de verdad: se comprueba la DECISION —que modo
# se elige y con que argumentos— porque es donde un cambio silencioso rompe el
# arranque en un sitio (el contenedor) sin que se note en el otro.

test_that("el modo se lee del entorno y tolera un valor sin sentido", {
  withr::local_envvar(c(SARA_UI = "app"));      expect_equal(modo_de_apertura(), "app")
  withr::local_envvar(c(SARA_UI = "navegador"));expect_equal(modo_de_apertura(), "navegador")
  withr::local_envvar(c(SARA_UI = "ninguno"));  expect_equal(modo_de_apertura(), "ninguno")
  # Sinonimos en ingles, porque es lo que se teclea sin pensar.
  withr::local_envvar(c(SARA_UI = "browser"));  expect_equal(modo_de_apertura(), "navegador")
  withr::local_envvar(c(SARA_UI = "none"));     expect_equal(modo_de_apertura(), "ninguno")
  # Mayusculas y espacios sobrantes no deberian cambiar nada.
  withr::local_envvar(c(SARA_UI = "  APP  "));  expect_equal(modo_de_apertura(), "app")
  # Un valor invalido cae al comportamiento por defecto en lugar de fallar: el
  # arranque de la aplicacion no puede depender de una errata en una variable.
  withr::local_envvar(c(SARA_UI = "azul"))
  expect_true(modo_de_apertura() %in% c("app", "ninguno"))
})

test_that("dentro de un contenedor no se intenta abrir nada", {
  # Es la unica decision con consecuencias reales: en la imagen de contenedor no
  # hay escritorio, y lanzar un navegador solo ensucia el log de arranque con un
  # error que nadie puede resolver.
  withr::local_envvar(c(SARA_UI = ""))
  expect_equal(modo_de_apertura(contenedor = TRUE, grafico = TRUE), "ninguno")
  expect_identical(lanzador_de_interfaz(modo = "ninguno"), FALSE)
})

test_that("sin entorno grafico tampoco", {
  withr::local_envvar(c(SARA_UI = ""))
  expect_equal(modo_de_apertura(contenedor = FALSE, grafico = FALSE), "ninguno")
})

test_that("una peticion explicita gana sobre la deteccion automatica", {
  # Quien pone SARA_UI=app sabe lo que quiere, aunque estemos en un contenedor:
  # puede haber reenvio de X11 o un navegador dentro de la imagen.
  withr::local_envvar(c(SARA_UI = "app"))
  expect_equal(modo_de_apertura(contenedor = TRUE, grafico = FALSE), "app")
})

test_that("se elige el primer navegador que exista, en orden de preferencia", {
  d <- withr::local_tempdir()
  segundo <- file.path(d, "segundo"); file.create(segundo)
  tercero <- file.path(d, "tercero"); file.create(tercero)
  rutas <- c(file.path(d, "no_existe"), segundo, tercero)
  expect_equal(buscar_navegador_de_aplicacion(rutas), segundo)
  expect_null(buscar_navegador_de_aplicacion(file.path(d, c("nada", "tampoco"))))
  expect_null(buscar_navegador_de_aplicacion(character(0)))
})

test_that("sin navegador compatible se avisa y se cae al de siempre", {
  expect_false(abrir_ventana_de_aplicacion("http://127.0.0.1:1", navegador = NULL))

  # El lanzador no puede quedarse callado: sin ventana y sin mensaje, el usuario
  # se queda mirando una terminal sin saber que la aplicacion ya esta servida.
  llamada <- NULL
  local_mocked_bindings(browseURL = function(url, ...) { llamada <<- url; invisible() },
                        .package = "utils")
  lb <- lanzador_de_interfaz(modo = "app", abrir = function(...) FALSE)
  expect_true(is.function(lb))
  expect_message(lb("http://127.0.0.1:3838"), "navegador por defecto")
  expect_equal(llamada, "http://127.0.0.1:3838")
})

test_that("cuando la ventana se abre no se recurre al navegador por defecto", {
  llamada <- NULL
  local_mocked_bindings(browseURL = function(url, ...) { llamada <<- url; invisible() },
                        .package = "utils")
  lb <- lanzador_de_interfaz(modo = "app", abrir = function(...) TRUE)
  expect_silent(lb("http://127.0.0.1:3838"))
  expect_null(llamada)
})

test_that("se pasan los argumentos que quitan los adornos del navegador", {
  d <- withr::local_tempdir()
  falso <- file.path(d, "navegador"); file.create(falso)
  recibidos <- NULL
  expect_true(abrir_ventana_de_aplicacion(
    "http://127.0.0.1:3838", navegador = falso, tamano = "1280,800",
    ejecutar = function(command, args, ...) { recibidos <<- args; 0L }))
  # --app es lo que convierte la pestana en ventana: sin el, sale un navegador
  # normal con barra de direcciones y el efecto se pierde entero.
  expect_true(any(grepl("^--app=http://127\\.0\\.0\\.1:3838$", recibidos)))
  expect_true(any(grepl("^--window-size=1280,800$", recibidos)))
  # Sin estos dos, un perfil recien creado abre pestanas de bienvenida encima.
  expect_true("--no-first-run" %in% recibidos)
  expect_true("--no-default-browser-check" %in% recibidos)
})

test_that("el navegador se lanza sin esperarlo", {
  # Si se esperase a que el navegador termine, runApp() no llegaria a servir
  # nada y la ventana se abriria contra un servidor que aun no existe.
  d <- withr::local_tempdir()
  falso <- file.path(d, "navegador"); file.create(falso)
  espera <- NULL
  abrir_ventana_de_aplicacion(
    "http://127.0.0.1:3838", navegador = falso,
    ejecutar = function(command, args, wait = TRUE, ...) { espera <<- wait; 0L })
  expect_false(espera)
})

test_that("la deteccion de contenedor no emite avisos fuera de Linux", {
  # readLines() sobre /proc/1/cgroup avisa cuando el fichero no existe, y ese
  # aviso acababa en la consola en cada arranque en macOS.
  expect_silent(en_contenedor())
})
