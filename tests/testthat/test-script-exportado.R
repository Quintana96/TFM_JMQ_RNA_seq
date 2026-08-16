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

# ── Disenos con variables sustitutas ────────────────────────────────────────
#
# El caso que fallaba: sva + prefiltrado automatico, que es el modo por defecto.
# El bloque de prefiltrado se emite ANTES del que crea las SV, pero interpolaba
# `design_code` (que ya las incluye), de modo que el script moria en la primera
# linea util con "object 'SV1' not found". Y el bloque de sva reestimaba las
# variables con `~ condition` porque leia `rv$design_base`, un campo que ningun
# sitio de la aplicacion rellenaba: aunque el orden hubiera sido correcto, el
# script no reproducia un ajuste con batch o con formula libre.

#' Reproduce lo que hace el observer de la pestana 4 con sva activado:
#' prefiltra, estima las SV sobre la matriz YA prefiltrada y ajusta con el
#' diseno ampliado, registrando `design_base` y `design_code` por separado.
fit_rv_sva <- function(counts, meta, num, den, batch = "lote", n_sv = 1L) {
  base_formula <- paste0("~ ", batch, " + condition")
  d <- build_design(meta, den, batch)
  cm_f <- prefilter_counts(counts, mode = "auto",
                           design = stats::model.matrix(d$formula, data = d$meta),
                           group = meta$condition)
  sv <- estimate_surrogate_vars(cm_f, meta, stats::as.formula(base_formula),
                                n_sv = n_sv, seed = ANALYSIS_SEED)
  if (is.null(sv$sv)) return(NULL)
  meta_sv <- cbind(meta, as.data.frame(sv$sv))
  full_formula <- paste(base_formula, "+",
                        paste(colnames(sv$sv), collapse = " + "))
  res <- run_deg(cm_f, meta_sv, method = "DESeq2", ref_level = den,
                 contrast_num = num, batch = batch, fdr = 0.05, shrink = FALSE,
                 design_formula = full_formula)
  list(counts = cm_f, meta = meta_sv, results = res$table, method = "DESeq2",
       design = res$design, design_code = res$design_code,
       design_base = base_formula, coef = res$coef, contrast = res$contrast,
       fdr = 0.05, lfc_threshold = 0, shrink = "ninguno", padj_method = "BH",
       prefilter = attr(cm_f, "prefilter"), batch = batch,
       ref_level = den, contrast_num = num,
       seeds = list(sva = ANALYSIS_SEED, n_sv = sv$n_sv))
}

test_that("con sva el prefiltrado no referencia variables que aun no existen", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  skip_if_not(requireNamespace("sva", quietly = TRUE), "sva no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts, batch = TRUE)
  rv <- fit_rv_sva(counts, meta, num = "trt", den = "ctrl")
  skip_if(is.null(rv), "sva no ha estimado variables sustitutas con estos datos")

  script <- build_deg_r_script(rv)
  lineas <- strsplit(script, "\n", fixed = TRUE)[[1]]

  # El diseno del prefiltrado es el BASE: sin SV, porque en la app las SV se
  # estiman despues, sobre la matriz ya prefiltrada.
  pref <- grep("^design <- model.matrix", lineas, value = TRUE)
  expect_length(pref, 1L)
  expect_false(grepl("SV", pref, fixed = TRUE))
  expect_true(grepl("lote", pref, fixed = TRUE))

  # Y ninguna linea puede mencionar SV1 antes de la que lo crea.
  crea_sv <- grep("meta <- cbind(meta, as.data.frame(sv))", lineas, fixed = TRUE)
  usa_sv  <- grep("SV1", lineas, fixed = TRUE)
  expect_length(crea_sv, 1L)
  expect_true(all(usa_sv[usa_sv < crea_sv] %in% grep("^#", lineas)))
})

test_that("con sva el modelo de svaseq es el diseno base real, no ~ condition", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  skip_if_not(requireNamespace("sva", quietly = TRUE), "sva no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts, batch = TRUE)
  rv <- fit_rv_sva(counts, meta, num = "trt", den = "ctrl")
  skip_if(is.null(rv), "sva no ha estimado variables sustitutas con estos datos")

  script <- build_deg_r_script(rv)
  mod <- grep("^mod  <- model.matrix", strsplit(script, "\n", fixed = TRUE)[[1]],
              value = TRUE)
  expect_length(mod, 1L)
  # El ajuste llevaba batch: reestimar las SV con `~ condition` daria otras
  # variables y otro resultado.
  expect_true(grepl("~ lote + condition", mod, fixed = TRUE))
  expect_false(grepl("model.matrix(~ condition,", mod, fixed = TRUE))
})

test_that("el script con sva parsea y se ejecuta de principio a fin", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  skip_if_not(requireNamespace("sva", quietly = TRUE), "sva no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts, batch = TRUE)
  rv <- fit_rv_sva(counts, meta, num = "trt", den = "ctrl")
  skip_if(is.null(rv), "sva no ha estimado variables sustitutas con estos datos")

  script <- build_deg_r_script(rv)
  expect_silent(parse(text = script))

  dir <- withr::local_tempdir()
  # Se escriben las entradas que el script DECLARA: la matriz sin prefiltrar y
  # el samplesheet SIN las columnas SV, que es justo lo que el script reestima.
  utils::write.table(counts, file.path(dir, "counts.tsv"), sep = "\t",
                     quote = FALSE, col.names = NA)
  utils::write.table(meta, file.path(dir, "meta.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  script_path <- file.path(dir, "analisis.R")
  writeLines(script, script_path)

  out <- withr::with_dir(dir, {
    tryCatch({
      env <- new.env(parent = globalenv())
      suppressMessages(suppressWarnings(sys.source(script_path, envir = env)))
      list(ok = TRUE, res = get("res", envir = env))
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  })

  expect_true(out$ok, info = if (!out$ok) out$msg else "")
  # Y reproduce el ajuste: mismas llamadas de significacion que la app.
  app_sig <- rv$results$gene[!is.na(rv$results$padj) & rv$results$padj <= 0.05]
  scr_sig <- rownames(out$res)[!is.na(out$res$padj) & out$res$padj <= 0.05]
  expect_setequal(sort(app_sig), sort(scr_sig))
})

test_that("sin sva el prefiltrado sigue usando el diseno con batch", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts, batch = TRUE)
  rv <- fit_rv(counts, meta, num = "trt", den = "ctrl")
  rv$batch <- "lote"
  rv$design <- rv$design_code <- "~ lote + condition"
  rv$prefilter <- list(mode = "filterByExpr", n_before = 400, n_after = 380)

  pref <- grep("^design <- model.matrix",
               strsplit(build_deg_r_script(rv), "\n", fixed = TRUE)[[1]],
               value = TRUE)
  # El arreglo no puede haber degradado el caso normal a `~ condition`.
  expect_true(grepl("~ lote + condition", pref, fixed = TRUE))
})

test_that("con Wilcoxon el prefiltrado usa group= y no una formula", {
  skip_if_not(requireNamespace("edgeR", quietly = TRUE), "edgeR no instalado")
  counts <- make_test_counts()
  meta <- make_test_meta(counts)
  rv <- fit_rv(counts, meta, num = "trt", den = "ctrl", method = "Wilcoxon")
  rv$design_code <- NULL
  rv$prefilter <- list(mode = "filterByExpr", n_before = 400, n_after = 380)

  script <- build_deg_r_script(rv)
  # Wilcoxon no ajusta modelo: interpolar su etiqueta ("sin modelo...") dentro de
  # model.matrix() producia un fichero que ni siquiera parsea. La etiqueta si
  # aparece en el comentario de cabecera, que es donde tiene que estar.
  expect_silent(parse(text = script))
  expect_true(grepl("filterByExpr(y, group = meta$condition)", script, fixed = TRUE))
  codigo <- grep("^\\s*#", strsplit(script, "\n", fixed = TRUE)[[1]],
                 value = TRUE, invert = TRUE)
  expect_false(any(grepl("sin modelo", codigo, fixed = TRUE)))
  expect_false(any(grepl("model.matrix", codigo, fixed = TRUE)))
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
