#' utils_swish.R
#' Swish (fishpond): analisis diferencial que propaga la incertidumbre de la
#' cuantificacion.
#'
#' Por que existe (docs/REVISION_ESTADISTICA.md, B9): la app soporta salmon y
#' kallisto pero trataba sus `est_counts` como conteos exactos. No lo son: son
#' estimaciones con incertidumbre, especialmente en genes cuyos transcritos
#' comparten secuencia. Swish (Zhu et al., NAR 2019) promedia el estadistico
#' sobre 20-30 REPLICAS INFERENCIALES (muestras de Gibbs en salmon, bootstrap en
#' kallisto) y controla mejor la FDR en genes con alta incertidumbre.
#'
#' Esto senalaba una asimetria real: elegir kallisto en vez de bowtie2 no
#' cambiaba nada estadisticamente, cuando deberia habilitar metodos mejores.
#' Ahora si lo hace, a cambio de pedir las replicas al cuantificar
#' (`--INFERENTIAL_REPS` en workflow.sh, por defecto 20).

#' TRUE si el directorio de cuantificacion tiene replicas inferenciales.
#'
#' salmon las deja en `aux_info/bootstrap/`; kallisto las guarda dentro del
#' `abundance.h5`. Sin ellas Swish no aporta nada sobre los demas motores, asi
#' que conviene detectarlo ANTES de ofrecerlo.
has_inferential_replicates <- function(output_dir, tool) {
  if (is.null(output_dir) || !nzchar(output_dir %||% "")) return(FALSE)
  aln <- file.path(output_dir, "03_alignments", tool)
  if (!dir.exists(aln)) return(FALSE)
  sdirs <- list.dirs(aln, recursive = FALSE, full.names = TRUE)
  if (!length(sdirs)) return(FALSE)
  if (identical(tool, "salmon")) {
    return(any(vapply(sdirs, function(d) {
      b <- file.path(d, "aux_info", "bootstrap")
      dir.exists(b) && length(list.files(b)) > 0
    }, logical(1))))
  }
  if (identical(tool, "kallisto")) {
    return(any(vapply(sdirs, function(d) {
      f <- file.path(d, "abundance.h5")
      if (!file.exists(f) || !requireNamespace("rhdf5", quietly = TRUE)) return(FALSE)
      ls <- tryCatch(rhdf5::h5ls(f)$name, error = function(e) character(0))
      "bootstrap" %in% ls
    }, logical(1))))
  }
  FALSE
}

#' Motor Swish sobre las replicas inferenciales de una ejecucion.
#'
#' A diferencia de los demas motores, este NO parte de la matriz de conteos: usa
#' los ficheros de cuantificacion, porque la incertidumbre esta ahi y no en la
#' matriz ya resumida. Por eso recibe `output_dir` y no `counts`.
#'
#' @param meta metadatos con sample_id y condition
#' @param output_dir directorio de la ejecucion del workflow
#' @param tool "salmon" o "kallisto"
#' @param annotation_file GFF/GTF para resumir a gen (opcional: sin el, el
#'   analisis queda a nivel de transcrito)
#' @param seed semilla para `fishpond::swish`, que permuta etiquetas de muestra.
#'   Sin ella el mismo input devuelve p-valores y q-valores distintos en cada
#'   ejecucion, de modo que el analisis no seria reproducible.
run_deg_swish <- function(meta, output_dir, tool, annotation_file = NULL,
                          ref_level = NULL, contrast_num = NULL, batch = NULL,
                          fdr = 0.05, n_perms = SWISH_NPERMS, seed = 1L) {
  info <- list(contrast = NA_character_, coef = NA_character_,
               n_levels = NA_integer_, shrink = "ninguno", padj_method = "qvalue",
               disp_data = NULL, cooks = NULL, coef_available = character(0),
               n_perms = n_perms, seed = seed)
  fail <- function(msg) c(list(table = NULL, error = msg), info)

  if (!requireNamespace("fishpond", quietly = TRUE)) {
    return(fail("fishpond no esta instalado."))
  }
  if (!requireNamespace("tximeta", quietly = TRUE) &&
      !requireNamespace("tximport", quietly = TRUE)) {
    return(fail("Hace falta tximport o tximeta para leer las replicas inferenciales."))
  }
  qfiles <- quant_files_for_run(output_dir, tool)
  if (is.null(qfiles)) return(fail("No se han encontrado ficheros de cuantificacion."))
  if (!has_inferential_replicates(output_dir, tool)) {
    return(fail(paste0(
      "Esta ejecucion no tiene replicas inferenciales, que son lo que Swish ",
      "necesita. Relanza el workflow con --INFERENTIAL_REPS 20 (salmon: ",
      "--numGibbsSamples, kallisto: -b).")))
  }

  d <- build_design(meta, ref_level, NULL)
  num <- if (!is.null(contrast_num) && nzchar(contrast_num %||% "")) contrast_num
         else utils::tail(d$levels, 1)
  den <- d$ref
  info$n_levels <- length(d$levels)

  out <- tryCatch({
    m <- d$meta
    keep <- as.character(m$condition) %in% c(num, den)
    m <- m[keep, , drop = FALSE]
    qf <- qfiles[match(m$sample_id, names(qfiles))]
    if (any(is.na(qf))) {
      stop("Hay muestras del samplesheet sin fichero de cuantificacion.")
    }
    # `txOut = TRUE` mantiene el nivel de transcrito: Swish trabaja sobre las
    # replicas inferenciales, que son por transcrito.
    txi <- tximport::tximport(qf, type = tool, txOut = TRUE,
                              dropInfReps = FALSE)
    if (is.null(txi$infReps) || !length(txi$infReps)) {
      stop("tximport no ha devuelto replicas inferenciales.")
    }
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = c(list(counts = txi$counts, abundance = txi$abundance,
                      length = txi$length), txi$infReps),
      colData = S4Vectors::DataFrame(m)
    )
    # Swish pide la variable de interes como factor de dos niveles.
    SummarizedExperiment::colData(se)$condition <-
      factor(as.character(m$condition), levels = c(den, num))

    # `scaleInfReps` tambien remuestrea, asi que la semilla envuelve todo el
    # bloque estocastico y no solo a swish. `with_seed` restaura el estado del
    # RNG al salir, para no alterar el resto de la sesion.
    y <- withr::with_seed(seed, {
      y0 <- fishpond::scaleInfReps(se)
      y0 <- fishpond::labelKeep(y0)
      y0 <- y0[SummarizedExperiment::mcols(y0)$keep, ]
      fishpond::swish(y0, x = "condition",
                      cov = if (!is.null(batch) && nzchar(batch %||% "") &&
                                batch %in% names(m)) batch else NULL,
                      nperms = n_perms)
    })
    mc <- as.data.frame(SummarizedExperiment::mcols(y))
    # Acceso por [[ ]] y con longitud explicita: swish no garantiza los nombres
    # de sus columnas entre versiones, y `mc$x` sobre una columna ausente
    # devolveria NULL (o peor, un partial match) y reciclaria en el data.frame.
    col <- function(nm) {
      v <- mc[[nm]]
      if (is.null(v)) rep(NA_real_, nrow(mc)) else as.numeric(v)
    }
    # baseMean se toma directamente de los conteos escalados: es la definicion
    # que usan los demas motores, y no depende de como llame swish a su media.
    base_mean <- tryCatch(
      as.numeric(rowMeans(SummarizedExperiment::assay(y, "counts"), na.rm = TRUE)),
      error = function(e) rep(NA_real_, nrow(mc)))
    tab <- data.frame(
      gene = rownames(y),
      baseMean = base_mean,
      log2FC = col("log2FC"),
      log2FC_shrunk = NA_real_,
      lfcSE = NA_real_,
      stat = col("stat"),
      pvalue = col("pvalue"),
      padj = col("qvalue"),
      stringsAsFactors = FALSE
    )
    list(table = tab, coef = paste0("condition", num))
  }, error = function(e) e)

  if (inherits(out, "error")) return(fail(conditionMessage(out)))
  info$coef     <- out$coef
  info$contrast <- paste(num, "vs", den)
  c(list(table = out$table, error = NULL), info)
}
