#' test-auditoria.R
#' Cada análisis diferencial debe dejar rastro en disco. Antes solo lo dejaba si
#' el usuario pulsaba descargar, de modo que una ejecución del pipeline tenía su
#' log y sus parámetros pero los análisis hechos sobre ella no dejaban ninguno.

rv_para_guardar <- function() {
  counts <- make_test_counts(n_genes = 150, n_per_group = 4)
  meta <- make_test_meta(counts)
  res <- run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                 contrast_num = "trt", fdr = 0.05, shrink = FALSE)
  list(results = res$table, meta = meta, method = "DESeq2", design = res$design,
       coef = res$coef, contrast = res$contrast, fdr = 0.05, lfc_threshold = 0,
       shrink = "ninguno", padj_method = "BH", prefilter = NULL,
       run_at = Sys.time(), seeds = list(sva = 1L),
       ref_level = "ctrl", contrast_num = "trt",
       counts_origin = list(tipo = "Matriz subida", ruta = "x.tsv", md5 = "abc"),
       counts_source = list(method = "tximport", ok = TRUE))
}

test_that("un análisis se persiste completo en 05_deg/<timestamp>", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  base <- withr::local_tempdir()
  rv <- suppressMessages(rv_para_guardar())

  d <- suppressMessages(persist_deg_analysis(rv, base_dir = base, outputs_dir = base))
  expect_false(is.null(d))
  expect_true(dir.exists(d))
  expect_true(grepl("05_deg", d, fixed = TRUE))

  for (f in c("deg_params.tsv", "resultados.tsv", "samplesheet.tsv",
              "informe.html", "analisis.R")) {
    expect_true(file.exists(file.path(d, f)), info = f)
  }

  # La tabla guardada es la COMPLETA, no solo los significativos: recortarla
  # impediria recalcular otros umbrales o rehacer el enriquecimiento después.
  tab <- utils::read.delim(file.path(d, "resultados.tsv"))
  expect_equal(nrow(tab), nrow(rv$results))
})

test_that("deg_params.tsv es legible con el mismo lector que run_params.tsv", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  base <- withr::local_tempdir()
  rv <- suppressMessages(rv_para_guardar())
  d <- suppressMessages(persist_deg_analysis(rv, base_dir = base, outputs_dir = base))

  # Mismo formato clave-valor que run_params.tsv, mismo lector.
  p <- read_run_params_file(d, "deg_params.tsv")
  expect_true(length(p) > 5)
  expect_identical(p[["Motor"]], "DESeq2")
  expect_identical(p[["numerador"]], "trt")
  expect_identical(p[["denominador"]], "ctrl")
})

test_that("el registro de auditoria acumula eventos sin perder los previos", {
  base <- withr::local_tempdir()
  expect_true(append_audit_log("deg_run", list(motor = "DESeq2", fdr = 0.05),
                               outputs_dir = base))
  expect_true(append_audit_log("deg_export", list(formato = "html"),
                               outputs_dir = base))

  al <- read_audit_log(base)
  expect_equal(nrow(al), 2)
  expect_identical(al$accion, c("deg_run", "deg_export"))
  expect_true(all(c("timestamp", "usuario", "accion", "detalles") %in% names(al)))
  expect_match(al$detalles[1], "motor=DESeq2")
})

test_that("los detalles no rompen el formato TSV", {
  base <- withr::local_tempdir()
  # Un valor con tabulador y salto de línea debe quedar en UNA línea.
  append_audit_log("prueba", list(nota = "linea1\nlinea2\tcon tab"), outputs_dir = base)
  al <- read_audit_log(base)
  expect_equal(nrow(al), 1)
  expect_false(grepl("\t", al$detalles[1], fixed = TRUE))
})

test_that("un fallo al escribir el registro no tumba el análisis", {
  # Ruta imposible: la función debe devolver FALSE, no lanzar.
  expect_false(append_audit_log("x", list(), outputs_dir = "/proc/no/existe/jamas"))
})

test_that("sin resultados no se persiste nada", {
  base <- withr::local_tempdir()
  expect_null(persist_deg_analysis(NULL, base_dir = base, outputs_dir = base))
  expect_null(persist_deg_analysis(list(results = NULL), base_dir = base,
                                   outputs_dir = base))
})
