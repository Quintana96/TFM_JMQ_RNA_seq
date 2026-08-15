#' utils_report.R
#' Informe reproducible y exportacion del script R equivalente.
#'
#' Por que existe (docs/REVISION_ESTADISTICA.md, B11): una aplicacion grafica
#' agrava el problema de trazabilidad que Pall et al. asocian al sesgo del campo,
#' porque al no haber script no queda registro de que se hizo. Dos salidas de
#' bajo coste y alto valor: un informe con todos los parametros y los
#' diagnosticos, y el script R equivalente al analisis ejecutado — que convierte
#' la app en herramienta de aprendizaje y permite auditar el resultado fuera de
#' ella.
#'
#' El informe se genera como HTML autocontenido escrito a mano, NO via
#' rmarkdown::render(): pandoc no esta disponible en todos los entornos donde
#' corre la app (no lo esta en esta maquina), y una descarga que falla la mitad
#' de las veces no sirve. El script R exportado si es texto plano, sin
#' dependencias.

#' Escapa texto para insertarlo en HTML.
html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

#' Tabla HTML de dos columnas a partir de una lista clave-valor.
html_kv_table <- function(x) {
  if (!length(x)) return("")
  rows <- vapply(names(x), function(k) {
    paste0("<tr><th>", html_escape(k), "</th><td>",
           html_escape(paste(x[[k]], collapse = ", ")), "</td></tr>")
  }, character(1))
  paste0("<table class=\"kv\">", paste(rows, collapse = ""), "</table>")
}

#' Parametros del analisis DEG en forma de lista ordenada, para el informe y el
#' script. Es la fuente unica de verdad de "que se ejecuto".
deg_run_parameters <- function(rv) {
  if (is.null(rv) || is.null(rv$results)) return(NULL)
  pf <- rv$prefilter
  list(
    "Fecha del analisis"       = format(rv$run_at %||% Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "Motor"                    = rv$method %||% "—",
    "Diseno"                   = rv$design %||% "~ condition",
    "Contraste"                = rv$contrast %||% "—",
    "Coeficiente testeado"     = rv$coef %||% "—",
    "FDR objetivo (alpha)"     = rv$fdr %||% 0.05,
    "Correccion multiple"      = rv$padj_method %||% "BH",
    "Umbral |log2FC| del test" = rv$lfc_threshold %||% 0,
    "Encogido de log2FC"       = rv$shrink %||% "ninguno",
    "Prefiltrado"              = if (is.null(pf)) "—" else
      paste0(pf$mode, ": ", pf$n_before, " -> ", pf$n_after, " genes"),
    "Muestras"                 = if (!is.null(rv$meta)) nrow(rv$meta) else NA,
    "Genes analizados"         = if (!is.null(rv$results)) nrow(rv$results) else NA,
    "Correccion solo grafica"  = rv$viz_note %||% "ninguna"
  )
}

#' Informe HTML autocontenido del analisis.
#'
#' @param rv `state$deg_rv`
#' @param diagnostics list(pi0, verdict, na_breakdown, cooks_dominant, replicability)
build_deg_report_html <- function(rv, diagnostics = NULL) {
  params <- deg_run_parameters(rv)
  if (is.null(params)) return(NULL)

  sig_n <- sum(!is.na(rv$results$padj) & rv$results$padj <= (rv$fdr %||% 0.05))
  up <- sum(!is.na(rv$results$padj) & rv$results$padj <= (rv$fdr %||% 0.05) &
              !is.na(rv$results$log2FC) & rv$results$log2FC > 0)

  diag_html <- ""
  if (!is.null(diagnostics)) {
    bits <- list()
    if (!is.null(diagnostics$verdict)) {
      bits[["Forma del histograma de p-valores"]] <-
        paste0(diagnostics$verdict$label, " — ", diagnostics$verdict$detail)
    }
    if (!is.null(diagnostics$pi0) && !is.na(diagnostics$pi0)) {
      bits[["pi0 (proporcion de nulas ciertas)"]] <- round(diagnostics$pi0, 4)
    }
    if (!is.null(diagnostics$na_breakdown)) {
      bits[["Genes sin p-valor ajustado"]] <-
        padj_na_breakdown_text(diagnostics$na_breakdown)
    }
    if (!is.null(diagnostics$cooks_dominant) && !is.na(diagnostics$cooks_dominant)) {
      bits[["Muestra que concentra los outliers de Cook"]] <- diagnostics$cooks_dominant
    }
    # El objeto de replicabilidad puede venir de un intento fallido (sin summary),
    # asi que se comprueba antes de leerlo en lugar de asumir que esta completo.
    r <- diagnostics$replicability
    if (!is.null(r)) {
      bits[["Replicabilidad (bootstrap)"]] <- if (!is.null(r$summary) &&
                                                  is.numeric(r$summary$mediana)) {
        paste0(r$interpretation$label, " — Spearman mediana ",
               round(r$summary$mediana[1], 3), " sobre ", r$n_ok, " remuestreos")
      } else {
        paste0("No estimada: ", r$error %||% "motivo desconocido")
      }
    }
    if (length(bits)) {
      diag_html <- paste0("<h2>Diagnosticos</h2>", html_kv_table(bits))
    }
  }

  top <- rv$results[!is.na(rv$results$padj), , drop = FALSE]
  top <- top[order(top$padj), , drop = FALSE]
  top <- utils::head(top, 25)
  cols <- intersect(c("gene", "baseMean", "log2FC", "log2FC_shrunk", "pvalue", "padj"),
                    names(top))
  top_rows <- if (nrow(top)) paste(vapply(seq_len(nrow(top)), function(i) {
    paste0("<tr>", paste(vapply(cols, function(c) {
      v <- top[[c]][i]
      paste0("<td>", html_escape(if (is.numeric(v)) signif(v, 4) else v), "</td>")
    }, character(1)), collapse = ""), "</tr>")
  }, character(1)), collapse = "") else ""

  si <- utils::capture.output(utils::sessionInfo())
  pkgs <- c("DESeq2", "edgeR", "limma", "apeglm", "IHW", "sva", "clusterProfiler",
            "fgsea", "tximport", "dearseq")
  pkg_versions <- stats::setNames(lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p))
    else "no instalado"
  }), pkgs)

  paste0(
'<!DOCTYPE html><html lang="es"><head><meta charset="utf-8">',
'<title>Informe de expresion diferencial</title><style>',
'body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:900px;',
'margin:2rem auto;padding:0 1rem;color:#20332A;line-height:1.5}',
'h1{border-bottom:3px solid #7BBF9A;padding-bottom:.3rem}',
'h2{margin-top:2rem;color:#244B34}',
'table{border-collapse:collapse;width:100%;margin:.5rem 0;font-size:.92rem}',
'th,td{border:1px solid #D9E2DC;padding:.35rem .5rem;text-align:left}',
'table.kv th{width:34%;background:#F2F7F4}',
'thead th{background:#E8F1EC}',
'.metric{display:inline-block;background:#F2F7F4;border:1px solid #D9E2DC;',
'border-radius:6px;padding:.5rem .9rem;margin:.25rem .4rem .25rem 0}',
'.metric b{display:block;font-size:1.3rem;color:#244B34}',
'pre{background:#F7F9F8;border:1px solid #E3EAE6;padding:.6rem;overflow-x:auto;',
'font-size:.8rem}',
'footer{margin-top:2.5rem;font-size:.82rem;color:#60756A;border-top:1px solid #D9E2DC;',
'padding-top:.6rem}</style></head><body>',
'<h1>Informe de expresion diferencial</h1>',
'<p><span class="metric"><b>', fmt_int(nrow(rv$results)), '</b>genes analizados</span>',
'<span class="metric"><b>', fmt_int(sig_n), '</b>significativos a FDR &le; ',
html_escape(rv$fdr %||% 0.05), '</span>',
'<span class="metric"><b>', fmt_int(up), ' / ', fmt_int(sig_n - up),
'</b>up / down</span></p>',
'<h2>Parametros del analisis</h2>', html_kv_table(params),
diag_html,
'<h2>Metadatos de las muestras</h2>',
if (!is.null(rv$meta)) paste0(
  '<table><thead><tr>',
  paste0('<th>', vapply(names(rv$meta), html_escape, character(1)), '</th>', collapse = ''),
  '</tr></thead><tbody>',
  paste(vapply(seq_len(nrow(rv$meta)), function(i) paste0('<tr>',
    paste0('<td>', vapply(rv$meta[i, ], function(v) html_escape(as.character(v)), character(1)),
           '</td>', collapse = ''), '</tr>'), character(1)), collapse = ''),
  '</tbody></table>') else '<p>Sin metadatos.</p>',
'<h2>Top 25 genes por p-valor ajustado</h2>',
'<table><thead><tr>',
paste0('<th>', vapply(cols, html_escape, character(1)), '</th>', collapse = ''),
'</tr></thead><tbody>', top_rows, '</tbody></table>',
'<h2>Versiones de los paquetes</h2>', html_kv_table(pkg_versions),
'<h2>sessionInfo()</h2><pre>', html_escape(paste(si, collapse = "\n")), '</pre>',
'<footer>Generado por RNA-seq Workflow Runner el ',
html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
'. Este informe recoge los parametros con los que se ejecuto el analisis; ',
'el script R equivalente se puede descargar aparte para reproducirlo fuera de la app.',
'</footer></body></html>')
}

#' Script R equivalente al analisis ejecutado.
#'
#' No es una transcripcion de la app: es el codigo minimo que reproduce el mismo
#' resultado usando directamente los paquetes de Bioconductor, para que se pueda
#' auditar y correr sin la app.
build_deg_r_script <- function(rv) {
  if (is.null(rv) || is.null(rv$results)) return(NULL)
  method <- rv$method %||% "DESeq2"
  design <- rv$design %||% "~ condition"
  fdr <- rv$fdr %||% 0.05
  lfc <- rv$lfc_threshold %||% 0
  coef <- rv$coef %||% ""
  # Denominador REAL del contraste ajustado. Antes se tomaba el primer nivel en
  # orden ALFABETICO, que solo coincide con el denominador elegido por
  # casualidad: en cuanto el usuario contrastaba, por ejemplo, "ctrl vs trt", el
  # script releveleaba a "ctrl" y luego pedia un coeficiente que no existe, de
  # modo que fallaba al ejecutarse. El contraste es parte del ajuste, no algo
  # que se pueda deducir de los datos a posteriori.
  ref <- rv$ref_level %||%
    # Respaldo para estados guardados antes de que se registrara ref_level:
    # el contraste se muestra como "num vs den".
    (if (!is.null(rv$contrast) && grepl(" vs ", rv$contrast, fixed = TRUE))
       sub("^.* vs ", "", rv$contrast) else NULL) %||%
    (if (!is.null(rv$meta) && "condition" %in% names(rv$meta))
       levels(as.factor(as.character(rv$meta$condition)))[1] else "")
  num <- rv$contrast_num %||%
    (if (!is.null(rv$contrast) && grepl(" vs ", rv$contrast, fixed = TRUE))
       sub(" vs .*$", "", rv$contrast) else "")
  seeds <- rv$seeds %||% list()
  pf <- rv$prefilter

  header <- c(
    "# Script equivalente al analisis ejecutado en RNA-seq Workflow Runner",
    paste0("# Generado el ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("# Contraste: ", rv$contrast %||% "—",
           "   |   Diseno: ", design,
           "   |   Motor: ", method),
    "#",
    "# Entradas que hay que proporcionar:",
    if (identical(method, "Swish"))
      "#   quant_dir   directorio 03_alignments de la ejecucion (salmon/kallisto)"
    else "#   counts.tsv  matriz de conteos (genes x muestras)",
    "#   meta.tsv    samplesheet con sample_id, condition y las covariables del diseno",
    "",
    "# Semilla de todo lo estocastico del analisis. Sin fijarla, los resultados",
    "# que dependen de permutaciones no son reproducibles.",
    paste0("set.seed(", seeds$sva %||% seeds$swish %||% 1L, ")"),
    "",
    if (identical(method, "Swish")) character(0) else
      'counts <- as.matrix(read.delim("counts.tsv", row.names = 1, check.names = FALSE))',
    'meta   <- read.delim("meta.tsv", stringsAsFactors = FALSE)',
    'meta   <- meta[match(colnames(counts), meta$sample_id), , drop = FALSE]',
    paste0('# El denominador del contraste es el nivel de referencia del factor.'),
    paste0('meta$condition <- relevel(factor(meta$condition), ref = "', ref, '")'),
    ""
  )

  # Variables sustitutas: la formula del diseno las menciona (SV1, SV2...) pero
  # no estan en meta.tsv, asi que hay que volver a estimarlas o el script no
  # correria. Se reproducen con la misma semilla y el mismo numero.
  sva_block <- if (!is.null(seeds$n_sv) && seeds$n_sv > 0) c(
    "# El diseno incluye variables sustitutas estimadas con sva. Para reproducir",
    "# el ajuste hay que volver a estimarlas: no viajan en el samplesheet.",
    "mod  <- model.matrix(~ condition, data = meta)",
    "mod0 <- model.matrix(~ 1, data = meta)",
    "cm_sv <- as.matrix(counts)[rowMeans(as.matrix(counts)) > 1, , drop = FALSE]",
    paste0("sv <- sva::svaseq(cm_sv, mod, mod0, n.sv = ", seeds$n_sv, ")$sv"),
    "colnames(sv) <- paste0(\"SV\", seq_len(ncol(sv)))",
    "meta <- cbind(meta, as.data.frame(sv))",
    ""
  ) else character(0)

  # Swish no parte de la matriz de conteos, asi que el prefiltrado no aplica:
  # emitirlo produciria un script que referencia un objeto `counts` inexistente.
  prefilter <- if (identical(method, "Swish")) character(0) else
    if (!is.null(pf) && identical(pf$mode, "filterByExpr")) c(
    "# Prefiltrado: filterByExpr usa el tamano del grupo mas pequeno",
    paste0("design <- model.matrix(", design, ", data = meta)"),
    "y <- edgeR::DGEList(counts = round(counts), group = meta$condition)",
    "keep <- edgeR::filterByExpr(y, design = design)",
    "counts <- counts[keep, , drop = FALSE]",
    ""
  ) else if (!is.null(pf)) c(
    paste0("# Prefiltrado manual: >= ", pf$min_count, " conteos en >= ",
           pf$min_samples, " muestras"),
    paste0("keep <- rowSums(counts >= ", pf$min_count, ") >= ", pf$min_samples),
    "counts <- counts[keep, , drop = FALSE]",
    ""
  ) else character(0)

  body <- switch(method,
    "DESeq2" = c(
      "library(DESeq2)",
      # IHW necesita S4Vectors ATACHADO, no solo cargado: sin esto falla con
      # "no se pudo encontrar la funcion mcols". Es el mismo tropiezo que tuvo
      # la propia app, asi que el script no puede omitirlo.
      if (identical(rv$padj_method, "IHW"))
        "library(IHW); library(S4Vectors)  # IHW necesita S4Vectors atachado"
      else character(0),
      paste0("dds <- DESeqDataSetFromMatrix(round(counts), meta, ", design, ")"),
      "dds <- DESeq(dds)",
      paste0("# alpha = FDR objetivo: calibra el filtrado independiente"),
      paste0("res <- results(dds, name = \"", coef, "\", alpha = ", fdr,
             if (lfc > 0) paste0(", lfcThreshold = ", lfc,
                                 ", altHypothesis = \"greaterAbs\"") else "",
             if (identical(rv$padj_method, "IHW")) ", filterFun = IHW::ihw" else "",
             ")"),
      if (!identical(rv$shrink %||% "ninguno", "ninguno")) c(
        "# El encogido cambia el log2FC pero no los p-valores",
        paste0("shr <- lfcShrink(dds, coef = \"", coef, "\", type = \"",
               rv$shrink, "\")")) else character(0),
      "res <- as.data.frame(res)"
    ),
    "edgeR" = c(
      "library(edgeR)",
      paste0("design <- model.matrix(", design, ", data = meta)"),
      "y <- DGEList(counts = round(counts), group = meta$condition)",
      "y <- normLibSizes(y)                 # edgeR v4; antes calcNormFactors",
      "fit <- glmQLFit(y, design)           # v4 estima las dispersiones aqui",
      if (lfc > 0)
        paste0("test <- glmTreat(fit, coef = \"", coef, "\", lfc = ", lfc, ")")
      else paste0("test <- glmQLFTest(fit, coef = \"", coef, "\")"),
      "res <- topTags(test, n = Inf, sort.by = \"none\")$table"
    ),
    "limma-voom" = c(
      "library(limma); library(edgeR)",
      paste0("design <- model.matrix(", design, ", data = meta)"),
      "y <- normLibSizes(DGEList(counts = round(counts), group = meta$condition))",
      "v <- voomWithQualityWeights(y, design)",
      "fit <- lmFit(v, design)",
      if (lfc > 0) c(
        paste0("fit <- treat(fit, lfc = ", lfc, ", robust = TRUE)"),
        paste0("res <- topTreat(fit, coef = \"", coef, "\", number = Inf, sort.by = \"none\")"))
      else c(
        "fit <- eBayes(fit, robust = TRUE)",
        paste0("res <- topTable(fit, coef = \"", coef, "\", number = Inf, sort.by = \"none\")"))
    ),
    "Wilcoxon" = c(
      "# Wilcoxon rank-sum sobre CPM: sin modelo, no ajusta covariables",
      "libs <- colSums(counts); cpm <- t(t(counts) / libs) * 1e6",
      "g <- as.character(meta$condition)",
      paste0("num <- \"", sub(" vs .*$", "", rv$contrast %||% ""), "\"; ",
             "den <- \"", sub("^.* vs ", "", rv$contrast %||% ""), "\""),
      "pv <- apply(cpm, 1, function(x) wilcox.test(x[g == num], x[g == den],",
      "                                            exact = FALSE)$p.value)",
      "res <- data.frame(pvalue = pv, padj = p.adjust(pv, method = \"BH\"))"
    ),
    "dearseq" = c(
      "library(dearseq)",
      "v2t <- model.matrix(~ condition, data = meta)[, -1, drop = FALSE]",
      "res <- dear_seq(exprmat = round(counts), variables2test = v2t,",
      "                which_test = \"asymptotic\", preprocessed = FALSE)$pvals"
    ),
    "Swish" = c(
      "library(fishpond); library(tximport); library(SummarizedExperiment)",
      "# Swish no parte de la matriz de conteos: la incertidumbre de la",
      "# cuantificacion vive en las replicas inferenciales de salmon/kallisto.",
      "# El analisis queda a nivel de TRANSCRITO (txOut = TRUE).",
      paste0('qfiles <- list.files("quant_dir", pattern = "',
             if (identical(rv$quant_tool %||% "", "kallisto")) "abundance.h5"
             else "quant.sf",
             '", recursive = TRUE, full.names = TRUE)'),
      "names(qfiles) <- basename(dirname(qfiles))",
      "qfiles <- qfiles[match(meta$sample_id, names(qfiles))]",
      paste0('txi <- tximport(qfiles, type = "', rv$quant_tool %||% "salmon",
             '", txOut = TRUE, dropInfReps = FALSE)'),
      "se <- SummarizedExperiment(",
      "  assays = c(list(counts = txi$counts, abundance = txi$abundance,",
      "                  length = txi$length), txi$infReps),",
      "  colData = S4Vectors::DataFrame(meta))",
      paste0('colData(se)$condition <- factor(as.character(meta$condition),',
             ' levels = c("', ref, '", "', num, '"))'),
      "y <- scaleInfReps(se)",
      "y <- labelKeep(y)",
      "y <- y[mcols(y)$keep, ]",
      paste0("y <- swish(y, x = \"condition\", nperms = ",
             seeds$swish_nperms %||% 30L, ")"),
      "res <- as.data.frame(mcols(y))",
      "# Swish reporta qvalue, no padj; se renombra para el resumen de abajo.",
      "res$padj <- res$qvalue"
    ),
    "# Motor no reconocido"
  )

  footer <- c(
    "",
    paste0("# Genes significativos a FDR <= ", fdr),
    paste0("sum(res$padj <= ", fdr, ", na.rm = TRUE)"),
    "",
    "sessionInfo()"
  )
  paste(c(header, sva_block, prefilter, body, footer), collapse = "\n")
}
