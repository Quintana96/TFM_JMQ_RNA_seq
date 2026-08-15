#' test-potencia-outliers.R
#' Dos decisiones que la app dejaba fuera del alcance del usuario: que hacer con
#' los outliers de Cook, y de donde salen los parametros del calculo de potencia.

test_that("estimate_power_params mide el CV y la profundidad de los datos", {
  skip_if_not(requireNamespace("edgeR", quietly = TRUE), "edgeR no instalado")
  counts <- make_test_counts(n_genes = 400, n_per_group = 4)
  meta <- make_test_meta(counts)

  p <- estimate_power_params(counts, meta)
  expect_false(is.null(p))
  expect_true(is.finite(p$cv))
  expect_gt(p$cv, 0)
  expect_true(is.finite(p$depth))
  expect_gt(p$depth, 0)
  expect_equal(p$n_por_grupo, 4)
})

test_that("los parametros medidos cambian la potencia frente a los supuestos", {
  skip_if_not(requireNamespace("edgeR", quietly = TRUE), "edgeR no instalado")
  skip_if_not(requireNamespace("RNASeqPower", quietly = TRUE), "RNASeqPower no instalado")
  counts <- make_test_counts(n_genes = 400, n_per_group = 4)
  meta <- make_test_meta(counts)
  p <- estimate_power_params(counts, meta)

  medida <- power_for_n(4, cv = p$cv, effect = 2, depth = p$depth)$power
  supuesta <- power_for_n(4, cv = 0.4, effect = 2, depth = 20)$power
  # No se exige una direccion concreta, solo que la diferencia sea material:
  # adivinar estos dos valores cambia la conclusion del calculo.
  expect_true(is.finite(medida) && is.finite(supuesta))
  expect_gt(abs(medida - supuesta), 0.05)
})

test_that("estimate_power_params tolera entradas degeneradas", {
  expect_null(estimate_power_params(NULL))
  expect_null(estimate_power_params(matrix(1:4, nrow = 4, ncol = 1)))
})

test_that("el tratamiento de outliers de Cook cambia el resultado", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_genes = 300, n_per_group = 4)
  meta <- make_test_meta(counts)
  # Un valor extremo en una sola muestra: DESeq2 lo marca y deja el gen sin padj.
  counts["gene0001", 5] <- 500000L

  na   <- suppressMessages(run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                                   contrast_num = "trt", fdr = 0.05, shrink = FALSE,
                                   outliers = "na"))
  keep <- suppressMessages(run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                                   contrast_num = "trt", fdr = 0.05, shrink = FALSE,
                                   outliers = "keep"))

  expect_false(is.null(na$table))
  expect_false(is.null(keep$table))
  # Ignorar el filtro de Cook devuelve al test genes que quedaban fuera.
  expect_lte(sum(is.na(keep$table$padj)), sum(is.na(na$table$padj)))
})

test_that("el modo de outliers se valida", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_genes = 100, n_per_group = 3)
  meta <- make_test_meta(counts)
  expect_error(run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                       contrast_num = "trt", outliers = "inventado"))
})
