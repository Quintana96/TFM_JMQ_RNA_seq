#' test-script-parametrizado.R
#' El script exportado tiene que PARSEAR para todas las combinaciones que la
#' interfaz permite, no solo para la que se probo al escribirlo.
#'
#' Existe porque la bateria anterior solo cubria DESeq2 con prefiltrado NULL, sin
#' batch, sin variables sustitutas y sin modo de outliers. Por ese hueco se colo
#' un fallo real: con Wilcoxon, run_deg() devolvia como `design` una etiqueta en
#' PROSA ("sin modelo (test de dos muestras...)") y el generador la interpolaba
#' dentro de model.matrix(), produciendo un script que no parseaba. Como el
#' prefiltrado por defecto es automatico, le ocurria a cualquiera que exportase
#' un analisis Wilcoxon.

ajusta <- function(method, prefilter_mode = "auto", batch = FALSE,
                   design_formula = NULL, outliers = "na", n_per_group = 4) {
  counts <- make_test_counts(n_genes = 150, n_per_group = n_per_group)
  meta <- make_test_meta(counts, batch = batch)
  cm <- if (identical(prefilter_mode, "ninguno")) counts else
    prefilter_counts(counts, mode = prefilter_mode, group = meta$condition,
                     min_count = 5, min_samples = 2)
  args <- list(cm, meta, method = method, ref_level = "ctrl", contrast_num = "trt",
               fdr = 0.05, shrink = FALSE, design_formula = design_formula,
               batch = if (batch) "lote" else NULL)
  if (identical(method, "DESeq2")) args$outliers <- outliers
  res <- suppressMessages(suppressWarnings(do.call(run_deg, args)))
  if (is.null(res$table)) return(NULL)
  list(results = res$table, meta = meta, method = method, design = res$design,
       design_code = res$design_code, coef = res$coef, contrast = res$contrast,
       fdr = 0.05, lfc_threshold = 0, shrink = "ninguno", padj_method = "BH",
       prefilter = attr(cm, "prefilter"), ref_level = "ctrl", contrast_num = "trt",
       batch = if (batch) "lote" else NULL, outliers = outliers,
       run_at = Sys.time(), seeds = list())
}

motores <- c("DESeq2", "edgeR", "limma-voom", "Wilcoxon",
             if (isTRUE(HAS_DEARSEQ)) "dearseq")

for (m in motores) {
  for (pf in c("auto", "manual", "ninguno")) {
    test_that(sprintf("el script parsea: %s con prefiltrado %s", m, pf), {
      skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
      rv <- ajusta(m, prefilter_mode = pf)
      skip_if(is.null(rv), paste("el motor", m, "no produjo tabla"))
      s <- build_deg_r_script(rv)
      expect_true(is.character(s) && nzchar(s))
      expect_silent(parse(text = s))
    })
  }
}

test_that("el script parsea con batch y con formula libre", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  rv <- ajusta("DESeq2", batch = TRUE)
  expect_silent(parse(text = build_deg_r_script(rv)))

  rv2 <- ajusta("DESeq2", batch = TRUE, design_formula = "~ lote + condition")
  expect_silent(parse(text = build_deg_r_script(rv2)))
})

test_that("el script parsea con los tres modos de outliers", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  for (o in c("na", "refit", "keep")) {
    rv <- ajusta("DESeq2", outliers = o)
    s <- build_deg_r_script(rv)
    expect_silent(parse(text = s))
    # Y el modo elegido tiene que aparecer: si no, dos analisis con resultados
    # distintos producirian el mismo script.
    if (identical(o, "keep")) expect_true(grepl("cooksCutoff = FALSE", s, fixed = TRUE))
    if (identical(o, "refit")) expect_true(grepl("minReplicatesForReplace", s, fixed = TRUE))
  }
})

test_that("el script de Swish no referencia un objeto counts inexistente", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  rv <- ajusta("DESeq2")
  rv$method <- "Swish"; rv$quant_tool <- "salmon"
  rv$seeds <- list(swish = 1L, swish_nperms = 30L)
  rv$prefilter <- NULL
  s <- build_deg_r_script(rv)
  expect_silent(parse(text = s))
  # Swish parte de las cuantificaciones, no de la matriz: no puede alinear el
  # samplesheet contra colnames(counts).
  expect_false(grepl("colnames(counts)", s, fixed = TRUE))
})

test_that("Wilcoxon normaliza con TMM tambien en el script", {
  skip_if_not(requireNamespace("edgeR", quietly = TRUE), "edgeR no instalado")
  rv <- ajusta("Wilcoxon")
  skip_if(is.null(rv), "Wilcoxon no produjo tabla")
  s <- build_deg_r_script(rv)
  # El motor usa normalized_cpm() (TMM). Emitir CPM por tamano de libreria
  # reintroduciria el sesgo de composicion que el propio codigo documenta.
  expect_true(grepl("normLibSizes", s, fixed = TRUE))
})

test_that("el diseno legible y el ejecutable no se confunden", {
  rv <- ajusta("Wilcoxon")
  skip_if(is.null(rv), "Wilcoxon no produjo tabla")
  # La etiqueta es prosa; el codigo, NULL. Interpolar la primera donde va el
  # segundo es lo que rompia el script.
  expect_match(rv$design, "sin modelo", fixed = TRUE)
  expect_null(rv$design_code)
})
