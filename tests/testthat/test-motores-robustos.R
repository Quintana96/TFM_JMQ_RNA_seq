#' test-motores-robustos.R
#' Los motores robustos (Wilcoxon, dearseq) no ajustan el mismo modelo que los
#' parametricos, y eso tiene que quedar declarado en lugar de que el banner y el
#' informe describan un ajuste que no ocurrio.

test_that("normalized_cpm corrige por composicion y no solo por libreria", {
  counts <- make_test_counts(n_genes = 300, n_per_group = 4)
  # Se induce un sesgo de composicion: unos pocos genes enormes en un grupo se
  # llevan las lecturas y hacen parecer reprimido al resto.
  counts[1:5, 5:8] <- counts[1:5, 5:8] * 200L

  tmm <- normalized_cpm(counts)
  libs <- colSums(counts); crudo <- t(t(counts) / libs) * 1e6

  expect_false(isTRUE(all.equal(tmm, crudo)))
  # El CPM crudo desplaza la mediana del grupo afectado; TMM lo corrige, asi que
  # la diferencia de medianas entre grupos debe ser MENOR tras TMM.
  med_dif <- function(m) abs(stats::median(m[-(1:5), 1:4]) - stats::median(m[-(1:5), 5:8]))
  expect_lt(med_dif(tmm), med_dif(crudo))
})

test_that("Wilcoxon declara que no ajusta un modelo", {
  counts <- make_test_counts(n_genes = 120, n_per_group = 4)
  meta <- make_test_meta(counts)
  res <- run_deg(counts, meta, method = "Wilcoxon", ref_level = "ctrl",
                 contrast_num = "trt", fdr = 0.05)
  expect_false(is.null(res$table))
  expect_match(res$design, "sin modelo", fixed = TRUE)
})

test_that("una formula libre con motor robusto avisa en lugar de fingir", {
  counts <- make_test_counts(n_genes = 120, n_per_group = 4)
  meta <- make_test_meta(counts, batch = TRUE)
  res <- run_deg(counts, meta, method = "Wilcoxon", ref_level = "ctrl",
                 contrast_num = "trt", fdr = 0.05,
                 design_formula = "~ lote + condition")
  # El diseno reportado NO puede ser la formula libre, porque no se ajusto.
  expect_false(identical(res$design, "~ lote + condition"))
  expect_false(is.null(res$design_warning))
  expect_match(res$design_warning, "no admite formulas", fixed = TRUE)
})

test_that("los motores parametricos siguen reportando su formula", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_genes = 150, n_per_group = 4)
  meta <- make_test_meta(counts, batch = TRUE)
  res <- suppressMessages(run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                                  contrast_num = "trt", fdr = 0.05, shrink = FALSE,
                                  design_formula = "~ lote + condition"))
  # deparse1 normaliza el espacio tras la virgulilla, asi que se compara el
  # contenido y no la cadena exacta.
  expect_match(res$design, "lote \\+ condition")
  expect_null(res$design_warning)
})

test_that("el bootstrap estratifica por condicion x batch cuando hay batch", {
  # Diseno exigente: una sola muestra por celda condicion x lote. Sin
  # estratificar por lote, un remuestreo puede perder un lote entero y ese
  # ajuste falla; como los fallos solo se cuentan, la estimacion queda apoyada
  # en los remuestreos faciles y sale optimista.
  meta <- data.frame(sample_id = paste0("s", 1:6),
                     condition = rep(c("ctrl", "trt"), each = 3),
                     lote = rep(c("A", "B", "C"), times = 2),
                     stringsAsFactors = FALSE)
  degenerado <- function(i) {
    is.null(i) || length(unique(meta$lote[i])) < 3 || any(table(meta$condition[i]) < 2)
  }
  withr::with_seed(1, {
    sin <- mean(replicate(200, degenerado(bootstrap_sample_indices(meta))))
    con <- mean(replicate(200, degenerado(bootstrap_sample_indices(meta, batch = "lote"))))
  })
  expect_gt(sin, 0.05)
  expect_equal(con, 0)
})

test_that("bootstrap_sample_indices no confunde un indice con un rango", {
  # Trampa clasica de R: sample(pos, 1) con pos de longitud 1 interpreta pos
  # como 1:pos. Con estratos de tamano 1 (frecuentes al cruzar condicion x
  # batch) eso devolveria indices inventados.
  meta <- data.frame(sample_id = paste0("s", 1:4),
                     condition = rep(c("ctrl", "trt"), each = 2),
                     lote = c("A", "B", "A", "B"), stringsAsFactors = FALSE)
  withr::with_seed(3, {
    for (k in 1:50) {
      i <- bootstrap_sample_indices(meta, batch = "lote")
      if (!is.null(i)) expect_true(all(i %in% seq_len(nrow(meta))))
    }
  })
})

test_that("el RLE se calcula sobre conteos normalizados por composicion", {
  counts <- make_test_counts(n_genes = 200, n_per_group = 4)
  counts[1:5, 5:8] <- counts[1:5, 5:8] * 200L

  norm <- rle_summary(counts, normalized = TRUE)
  crudo <- rle_summary(counts, normalized = FALSE)
  expect_false(isTRUE(all.equal(norm$med, crudo$med)))
  # El RLE diagnostica la normalizacion APLICADA: sobre conteos normalizados,
  # las medianas deben quedar mas cerca de cero que sobre CPM crudo.
  expect_lt(stats::median(abs(norm$med)), stats::median(abs(crudo$med)))
})
