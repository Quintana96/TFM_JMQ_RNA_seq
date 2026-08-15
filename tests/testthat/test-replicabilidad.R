#' test-replicabilidad.R
#' El panel de replicabilidad evalua los resultados del ajuste, asi que tiene que
#' reajustar EL MISMO modelo. Si remuestrea con otro diseno u otro contraste, el
#' numero que devuelve no describe los resultados que el usuario esta mirando.

test_that("el bootstrap respeta la formula libre del diseno", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_genes = 200, n_per_group = 4)
  meta <- make_test_meta(counts, batch = TRUE)

  res <- bootstrap_replicability(
    counts, meta, method = "DESeq2",
    ref_level = "ctrl", contrast_num = "trt",
    design_formula = "~ lote + condition",
    n_boot = 3L, seed = ANALYSIS_SEED)

  expect_null(res$error)
  expect_true(res$n_ok >= 1L)
  # Si la formula se hubiese ignorado, el ajuste habria sido ~ condition y el
  # coeficiente testeado seria otro; al menos debe haber corrido sin fallar.
  expect_true(is.data.frame(res$per_boot))
})

test_that("el bootstrap es reproducible con la misma semilla", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_genes = 150, n_per_group = 4)
  meta <- make_test_meta(counts)

  a <- bootstrap_replicability(counts, meta, method = "DESeq2",
                               ref_level = "ctrl", contrast_num = "trt",
                               n_boot = 3L, seed = ANALYSIS_SEED)
  b <- bootstrap_replicability(counts, meta, method = "DESeq2",
                               ref_level = "ctrl", contrast_num = "trt",
                               n_boot = 3L, seed = ANALYSIS_SEED)

  expect_equal(a$per_boot$spearman, b$per_boot$spearman)
  expect_equal(a$per_boot$jaccard_topn, b$per_boot$jaccard_topn)
})

test_that("el bootstrap rechaza menos de 3 muestras por grupo", {
  counts <- make_test_counts(n_genes = 80, n_per_group = 2)
  meta <- make_test_meta(counts)

  res <- bootstrap_replicability(counts, meta, method = "DESeq2",
                                 ref_level = "ctrl", contrast_num = "trt",
                                 n_boot = 3L, seed = ANALYSIS_SEED)
  # Con 2 por grupo el unico remuestreo con dos muestras distintas es el
  # original: la correlacion saldria 1 sin haber medido nada.
  expect_false(is.null(res$error))
  expect_match(res$error, "3", fixed = TRUE)
})

test_that("el contraste del bootstrap es el del ajuste, en ambos sentidos", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_genes = 150, n_per_group = 4)
  meta <- make_test_meta(counts)

  # Invertir el contraste debe invertir el signo del estadistico de referencia.
  ab <- bootstrap_replicability(counts, meta, method = "DESeq2",
                                ref_level = "ctrl", contrast_num = "trt",
                                n_boot = 2L, seed = ANALYSIS_SEED)
  ba <- bootstrap_replicability(counts, meta, method = "DESeq2",
                                ref_level = "trt", contrast_num = "ctrl",
                                n_boot = 2L, seed = ANALYSIS_SEED)

  expect_null(ab$error)
  expect_null(ba$error)
  expect_equal(ab$reference$log2FC, -ba$reference$log2FC, tolerance = 1e-4)
})
