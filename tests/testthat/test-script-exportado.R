#' test-script-exportado.R
#' El script R descargable es el entregable de reproducibilidad de la app: si no
#' se ejecuta, o se ejecuta pero reproduce OTRO analisis, no sirve para nada.
#'
#' El caso critico es el contraste cuyo denominador no es el primer nivel en
#' orden alfabetico ("ctrl vs trt": el denominador es "trt"). El script tomaba el
#' primer nivel alfabetico como referencia, releveleaba mal y pedia despues un
#' coeficiente inexistente, de modo que fallaba al ejecutarse.

fit_rv <- function(counts, meta, num, den, method = "DESeq2") {
  res <- run_deg(counts, meta, method = method, ref_level = den,
                 contrast_num = num, fdr = 0.05, shrink = FALSE)
  list(results = res$table, meta = meta, method = method,
       design = res$design, coef = res$coef, contrast = res$contrast,
       fdr = 0.05, lfc_threshold = 0, shrink = "ninguno",
       padj_method = "BH", prefilter = NULL,
       ref_level = den, contrast_num = num, seeds = list())
}

test_that("el script usa el denominador real del contraste, no el alfabetico", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts)

  # Denominador "trt": NO es el primer nivel alfabetico (lo es "ctrl").
  rv <- fit_rv(counts, meta, num = "ctrl", den = "trt")
  script <- build_deg_r_script(rv)

  expect_true(grepl('relevel(factor(meta$condition), ref = "trt")', script,
                    fixed = TRUE))
  expect_false(grepl('ref = "ctrl"', script, fixed = TRUE))
  # Y el coeficiente pedido debe existir en ese diseno relevelado.
  expect_true(grepl(rv$coef, script, fixed = TRUE))
})

test_that("el script exportado se ejecuta y reproduce las llamadas de significacion", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts)
  rv <- fit_rv(counts, meta, num = "ctrl", den = "trt")   # el caso que fallaba

  dir <- withr::local_tempdir()
  utils::write.table(counts, file.path(dir, "counts.tsv"), sep = "\t",
                     quote = FALSE, col.names = NA)
  utils::write.table(meta, file.path(dir, "meta.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  script_path <- file.path(dir, "analisis.R")
  writeLines(build_deg_r_script(rv), script_path)

  # Se ejecuta de verdad, en el directorio con las entradas que el script declara.
  out <- withr::with_dir(dir, {
    tryCatch({
      env <- new.env(parent = globalenv())
      suppressMessages(sys.source(script_path, envir = env))
      list(ok = TRUE, res = get("res", envir = env))
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  })

  expect_true(out$ok, info = if (!out$ok) out$msg else "")

  # Mismas llamadas de significacion que el ajuste de la app.
  app_sig <- rv$results$gene[!is.na(rv$results$padj) & rv$results$padj <= 0.05]
  scr_sig <- rownames(out$res)[!is.na(out$res$padj) & out$res$padj <= 0.05]
  expect_setequal(sort(app_sig), sort(scr_sig))
})

test_that("el script declara la semilla y no deja motores sin rama", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts)
  rv <- fit_rv(counts, meta, num = "trt", den = "ctrl")

  expect_match(build_deg_r_script(rv), "set.seed(", fixed = TRUE)

  # Swish es lanzable desde la interfaz, asi que no puede caer en el generico
  # "# Motor no reconocido".
  rv_swish <- rv
  rv_swish$method <- "Swish"
  rv_swish$quant_tool <- "salmon"
  rv_swish$seeds <- list(swish = 1L, swish_nperms = 30L)
  s <- build_deg_r_script(rv_swish)
  expect_false(grepl("Motor no reconocido", s, fixed = TRUE))
  expect_true(grepl("swish(", s, fixed = TRUE))
  expect_true(grepl("txOut = TRUE", s, fixed = TRUE))
})

test_that("con IHW el script atacha S4Vectors", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts)
  rv <- fit_rv(counts, meta, num = "trt", den = "ctrl")
  rv$padj_method <- "IHW"

  s <- build_deg_r_script(rv)
  # Sin S4Vectors atachado, IHW falla con "no se pudo encontrar la funcion mcols".
  expect_true(grepl("library(S4Vectors)", s, fixed = TRUE))
  expect_true(grepl("filterFun = IHW::ihw", s, fixed = TRUE))
})
