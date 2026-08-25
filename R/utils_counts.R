#' utils_counts.R
#' Carga de matrices de conteos generadas por el workflow (bowtie2, salmon, kallisto).

#' Combina quant.sf / abundance.tsv en una matriz de conteos (genes x muestras)
load_quant_counts <- function(qfiles, tool) {
  mats <- lapply(qfiles, function(f) {
    x <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(x)) return(NULL)

    if (identical(tool, "salmon") && all(c("Name", "NumReads") %in% names(x))) {
      vals <- x$NumReads
      names(vals) <- x$Name
      return(vals)
    }
    if (identical(tool, "kallisto") && all(c("target_id", "est_counts") %in% names(x))) {
      vals <- x$est_counts
      names(vals) <- x$target_id
      return(vals)
    }
    NULL
  })
  valid <- !vapply(mats, is.null, logical(1))
  if (!any(valid)) return(NULL)

  mats <- mats[valid]
  features <- Reduce(union, lapply(mats, names))
  out <- matrix(0, nrow = length(features), ncol = length(mats),
                dimnames = list(features, names(qfiles)[valid]))
  for (nm in names(mats)) {
    vals <- mats[[nm]]
    out[names(vals), nm] <- vals
  }
  round(out)
}

#' Carga la matriz count_matrix.tsv generada por featureCounts (bowtie2)
load_count_matrix_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) return(NULL)

  header_idx <- grep("^Geneid\\t|^geneid\\t|^Gene\\t|^gene\\t", lines, ignore.case = TRUE)[1]
  if (is.na(header_idx)) header_idx <- 1L

  df <- tryCatch(
    read.delim(
      text = paste(lines[header_idx:length(lines)], collapse = "\n"),
      header = TRUE, row.names = 1, check.names = FALSE, comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)

  names(df) <- sub("\\.bam$", "", basename(names(df)))
  df
}

#' Localiza los ficheros de cuantificación de salmon/kallisto de una ejecución.
quant_files_for_run <- function(output_dir, tool) {
  aln_dir <- file.path(output_dir, "03_alignments", tool)
  if (!dir.exists(aln_dir)) return(NULL)
  sdirs <- list.dirs(aln_dir, recursive = FALSE, full.names = TRUE)
  if (!length(sdirs)) return(NULL)
  qfiles <- if (identical(tool, "salmon")) file.path(sdirs, "quant.sf")
            else file.path(sdirs, "abundance.h5")
  if (!all(file.exists(qfiles)) && !identical(tool, "salmon")) {
    qfiles <- file.path(sdirs, "abundance.tsv")
  }
  valid <- file.exists(qfiles)
  if (!any(valid)) return(NULL)
  qfiles <- qfiles[valid]
  names(qfiles) <- basename(sdirs[valid])
  qfiles
}

#' Identificadores de transcrito presentes en un fichero de cuantificación.
#' Sirve para comprobar contra que anotación casan ANTES de llamar a tximport.
quant_tx_ids <- function(qfile, tool) {
  if (!length(qfile) || !file.exists(qfile[1])) return(character(0))
  f <- qfile[1]
  if (grepl("[.]h5$", f)) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) return(character(0))
    return(tryCatch(as.character(rhdf5::h5read(f, "aux/ids")),
                    error = function(e) character(0)))
  }
  x <- tryCatch(read.delim(f, check.names = FALSE, nrows = -1),
                error = function(e) NULL)
  if (is.null(x)) return(character(0))
  if (identical(tool, "salmon") && "Name" %in% names(x)) return(as.character(x$Name))
  if ("target_id" %in% names(x)) return(as.character(x$target_id))
  character(0)
}

#' Carga la matriz de conteos generada por el workflow para una tool dada.
#'
#' Para salmon/kallisto se intenta la ruta correcta (tximport con `tx2gene`) y
#' solo se cae a la suma cruda de `est_counts` si no hay anotación utilizable.
#'
#' Por qué importa (docs/REVISION_ESTADISTICA.md, A8): la versión anterior
#' llamaba a `tximport()` SIN `tx2gene` y sin `txOut = TRUE`, lo que falla
#' siempre ("tximport failed at summarizing to the gene-level"), y el fallo
#' quedaba escondido en un `tryCatch` que devolvía NULL. Resultado: la rama de
#' tximport nunca se ejecutaba con exito y todos los análisis con salmon o
#' kallisto usaban `round(est_counts)`, perdiendo los offsets de longitud media
#' de transcrito por muestra — los que corrigen los cambios de longitud efectiva
#' de gen por uso diferencial de isoformas.
#'
#' Ahora la degradación es explícita: el resultado lleva un atributo
#' "counts_source" con lo que se ha hecho y por qué, para que la interfaz lo
#' muestre en lugar de fingir que se uso la ruta buena.
#'
#' @param annotation_file GFF/GTF con el que construir `tx2gene`. Sin el, no hay
#'   forma de resumir transcritos a genes.
load_counts_from_workflow <- function(output_dir, tool, annotation_file = NULL) {
  # Las dos rutas de alineamiento acaban en featureCounts, así que la matriz se
  # lee igual: lo que cambia es quién colocó las lecturas, no cómo se cuentan.
  if (tool %in% c("bowtie2", "subjunc")) {
    m <- load_count_matrix_tsv(file.path(output_dir, "04_counts", "count_matrix.tsv"))
    if (!is.null(m)) attr(m, "counts_source") <- list(
      method = "featureCounts",
      detail = paste0("Conteos por gen de featureCounts sobre alineamientos de ", tool, "."),
      ok = TRUE)
    return(m)
  }
  if (!tool %in% c("salmon", "kallisto")) return(NULL)

  qfiles <- quant_files_for_run(output_dir, tool)
  if (is.null(qfiles)) return(NULL)

  fallback <- function(reason) {
    m <- load_quant_counts(qfiles, tool)
    if (!is.null(m)) attr(m, "counts_source") <- list(
      method = "est_counts redondeados", ok = FALSE,
      detail = paste0(reason, " Se han sumado los conteos estimados por ",
                      "transcrito sin resumir a gen ni aplicar los offsets de ",
                      "longitud: los resultados pueden estar sesgados en genes ",
                      "con varias isoformas."))
    m
  }

  if (!isTRUE(HAS_TXIMPORT)) return(fallback("tximport no está instalado."))

  t2g_all <- build_tx2gene_from_annotation(annotation_file)
  if (is.null(t2g_all)) {
    return(fallback(paste0(
      "No se ha podido construir el mapa transcrito-gen",
      if (is.null(annotation_file) || !nzchar(annotation_file %||% ""))
        " (no se indico fichero de anotación)." else " desde la anotación indicada.")))
  }
  ids <- quant_tx_ids(qfiles, tool)
  if (!length(ids)) return(fallback("No se han podido leer los identificadores de transcrito."))

  pick <- pick_tx2gene_for_quant(t2g_all, ids)
  # Umbral: por debajo de la mitad de los transcritos emparejados, la anotación
  # no corresponde al transcriptoma usado y resumir a gen daría basura.
  if (is.null(pick$tx2gene) || !is.finite(pick$rate) || pick$rate < 0.5) {
    return(fallback(paste0(
      "La anotación solo casa con ", pick$n_matched, " de ", pick$n_tx,
      " transcritos (", round(100 * (pick$rate %||% 0), 1),
      " %): probablemente no corresponde al transcriptoma con el que se cuantifico.")))
  }

  # ignoreTxVersion = FALSE a propósito: las variantes de identificador ya se han
  # resuelto en pick_tx2gene_for_quant(), que elige la convencion que realmente
  # casa. Dejarlo en TRUE hace que tximport corte por el primer punto y descarte
  # en silencio los transcritos cuyo ID lleva un accession con versión
  # ("cds-AAC73112.1"), colapsando la matriz a un puñado de genes.
  txi <- tryCatch(
    tximport::tximport(qfiles, type = tool, tx2gene = pick$tx2gene,
                       countsFromAbundance = "lengthScaledTPM",
                       ignoreTxVersion = FALSE),
    error = function(e) e
  )
  if (inherits(txi, "error")) {
    return(fallback(paste0("tximport fallo: ", conditionMessage(txi), ".")))
  }
  m <- round(txi$counts)
  attr(m, "counts_source") <- list(
    method = "tximport (lengthScaledTPM)", ok = TRUE,
    detail = paste0("Transcritos resumidos a gen con la anotación indicada: ",
                    pick$n_matched, " de ", pick$n_tx, " transcritos emparejados (",
                    round(100 * pick$rate, 1), " %), convencion de identificador '",
                    pick$alias_type, "'."))
  m
}
