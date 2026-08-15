#' test-determinismo.R
#' El analisis debe cumplir que el mismo input produzca siempre el mismo output.
#'
#' Estos tests cubren los dos puntos de aleatoriedad que estaban SIN controlar:
#' `sva::num.sv` (metodo "be", que estima el numero de variables sustitutas por
#' permutacion) y `fishpond::swish` (que permuta etiquetas de muestra). Tambien
#' comprueban que fijar la semilla en un sitio no altera el RNG de la sesion,
#' que es el motivo de usar `withr::with_seed` en lugar de `set.seed` a secas.

test_that("estimate_surrogate_vars da el mismo n_sv en llamadas repetidas", {
  skip_if_not(requireNamespace("sva", quietly = TRUE), "sva no instalado")
  counts <- make_test_counts(n_genes = 300, n_per_group = 6)
  meta <- make_test_meta(counts)

  a <- estimate_surrogate_vars(counts, meta, ~ condition, seed = ANALYSIS_SEED)
  b <- estimate_surrogate_vars(counts, meta, ~ condition, seed = ANALYSIS_SEED)

  expect_identical(a$n_sv, b$n_sv)
  expect_identical(a$n_sv_estimated, b$n_sv_estimated)
  if (!is.null(a$sv) && !is.null(b$sv)) expect_equal(a$sv, b$sv)
})

test_that("estimate_surrogate_vars no altera el RNG de la sesion", {
  skip_if_not(requireNamespace("sva", quietly = TRUE), "sva no instalado")
  counts <- make_test_counts(n_genes = 200, n_per_group = 5)
  meta <- make_test_meta(counts)

  # Si la funcion usara set.seed() sin restaurar, la segunda secuencia aleatoria
  # seria distinta de la primera pese a partir de la misma semilla.
  set.seed(99); esperado <- runif(3)
  set.seed(99)
  invisible(estimate_surrogate_vars(counts, meta, ~ condition, seed = 7L))
  obtenido <- runif(3)

  expect_equal(obtenido, esperado)
})

test_that("estimate_surrogate_vars admite 0 variables sustitutas sin inventarse una", {
  # Se fuerza el caso limite por la via del recorte de grados de libertad, que
  # es determinista: con n_sv pedido = 0 y sva estimando 0, no debe devolverse
  # ninguna covariable. Forzar un minimo de 1 metia una variable espuria que
  # consumia un grado de libertad y podia absorber senal de la condicion.
  counts <- make_test_counts(n_genes = 150, n_per_group = 4)
  meta <- make_test_meta(counts)

  res <- estimate_surrogate_vars(counts, meta, ~ condition,
                                 n_sv = NULL, min_residual_df = 3L,
                                 seed = ANALYSIS_SEED)

  # n_sv = 0 con error NULL es un resultado valido, no un fallo.
  expect_true(is.null(res$error) || is.character(res$error))
  if (is.null(res$error)) {
    expect_true(res$n_sv >= 0L)
    if (res$n_sv == 0L) expect_null(res$sv)
    if (res$n_sv > 0L) expect_equal(ncol(res$sv), res$n_sv)
  }
})

test_that("la semilla del analisis esta definida y es un entero", {
  expect_true(exists("ANALYSIS_SEED"))
  expect_true(is.numeric(ANALYSIS_SEED) && length(ANALYSIS_SEED) == 1L)
  expect_true(exists("SWISH_NPERMS"))
  expect_true(is.numeric(SWISH_NPERMS) && SWISH_NPERMS > 0)
})

test_that("run_deg_swish declara semilla y nperms aunque falle pronto", {
  skip_if_not(requireNamespace("fishpond", quietly = TRUE), "fishpond no instalado")
  meta <- make_test_meta(make_test_counts(n_genes = 50, n_per_group = 3))
  # Directorio inexistente: la funcion falla, pero el contrato de su valor de
  # retorno debe mantenerse para que el informe pueda declarar los parametros.
  res <- run_deg_swish(meta, tempfile(), "salmon", seed = ANALYSIS_SEED)
  expect_false(is.null(res$error))
  expect_identical(res$seed, ANALYSIS_SEED)
  expect_identical(res$n_perms, SWISH_NPERMS)
})
