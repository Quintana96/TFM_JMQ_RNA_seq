#' testthat.R
#' Lanzador de la batería de tests.
#'
#' La app no es un paquete de R, así que no se usa `test_check()`. Desde la raíz
#' del proyecto:
#'
#'   Rscript tests/testthat.R
#'
#' Los ficheros `helper-*.R` de tests/testthat cargan `global.R` y `R/*.R` antes
#' de ejecutar los tests, replicando el orden en que Shiny los sourcea.

library(testthat)

test_dir_path <- if (dir.exists("tests/testthat")) "tests/testthat" else "testthat"
if (!dir.exists(test_dir_path)) {
  stop("Ejecuta este script desde la raíz del proyecto: Rscript tests/testthat.R")
}

testthat::test_dir(test_dir_path, reporter = "summary", stop_on_failure = TRUE)
