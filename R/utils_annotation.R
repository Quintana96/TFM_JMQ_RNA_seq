#' utils_annotation.R
#' Funciones puras (sin Shiny) para extraer del GFF/GTF de anotacion lo que el
#' analisis necesita: el mapa transcrito -> gen para tximport y los genes de
#' rRNA para la metrica de QC.
#'
#' Criterio de identificadores: el workflow llama a
#' `featureCounts -t gene -g locus_tag`, asi que la ruta de alineamiento produce
#' IDs de gen = locus_tag. Para que las dos rutas (alineamiento y
#' pseudoalineamiento) den los MISMOS identificadores, el gen canonico de
#' `tx2gene` es tambien el locus_tag.

#' Cache de la anotacion parseada.
#'
#' Parsear el GFF es CARO: 7 segundos para los 9.523 features de E. coli, y un
#' GTF humano es dos ordenes de magnitud mayor. Cuatro funciones distintas lo
#' necesitan (`tx2gene`, longitudes de gen, deteccion de splicing, genes de
#' rRNA), invocadas desde cuatro puntos de la aplicacion, asi que sin cache el
#' mismo fichero se parseaba cuatro veces por sesion como minimo, y de nuevo en
#' cada invalidacion reactiva.
#'
#' Se guarda una sola entrada: el patron de uso es una anotacion por ejecucion, y
#' el data.frame de un GTF humano ocupa cientos de MB. La clave incluye mtime y
#' tamano para que un fichero modificado se vuelva a leer.
.annotation_cache <- new.env(parent = emptyenv())

#' Lee un GFF/GTF y devuelve sus atributos como data.frame.
#' Devuelve NULL (sin error) si no se puede leer, para que quien llame decida.
read_annotation_features <- function(path, use_cache = TRUE) {
  if (is.null(path) || !length(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  if (!requireNamespace("rtracklayer", quietly = TRUE)) return(NULL)

  info <- file.info(path)
  key <- paste(normalizePath(path, mustWork = FALSE), info$mtime, info$size, sep = "|")
  if (isTRUE(use_cache) && identical(.annotation_cache$key, key)) {
    return(.annotation_cache$value)
  }

  df <- tryCatch({
    gr <- rtracklayer::import(path)
    out <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
    out$.width <- BiocGenerics::width(gr)
    out
  }, error = function(e) NULL)

  if (isTRUE(use_cache) && !is.null(df)) {
    .annotation_cache$key <- key
    .annotation_cache$value <- df
  }
  df
}

#' Primera columna no vacia de entre varias candidatas.
first_present <- function(df, cands) {
  for (c in cands) {
    if (c %in% names(df)) {
      v <- as.character(df[[c]])
      if (any(!is.na(v) & nzchar(v))) return(v)
    }
  }
  rep(NA_character_, nrow(df))
}

#' Construye el mapa transcrito -> gen a partir de la anotacion.
#'
#' El problema real no es leer el GFF, es que los nombres de los transcritos que
#' produce salmon/kallisto son los encabezados del FASTA de transcriptoma, y esos
#' dependen de con que herramienta se extrajo (gffread, NCBI cds_from_genomic,
#' Ensembl Bacteria...). Por eso NO se adivina una convencion: se emite una fila
#' por cada ALIAS plausible del transcrito, y luego
#' `pick_tx2gene_for_quant()` elige el subconjunto que de verdad casa con los
#' identificadores del fichero de cuantificacion.
#'
#' @return data.frame(TXNAME, GENEID, alias_type) o NULL
build_tx2gene_from_annotation <- function(path, tx_types = NULL) {
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df) || !"type" %in% names(df)) return(NULL)
  types <- as.character(df$type)
  # Features que representan un transcrito. En procariotas el analisis se hace
  # tipicamente a nivel de CDS; se incluyen tambien los RNA no codificantes.
  tx_types <- tx_types %||% c("CDS", "mRNA", "transcript", "ncRNA", "rRNA",
                              "tRNA", "tmRNA", "antisense_RNA", "RNase_P_RNA",
                              "SRP_RNA", "exon")
  sub <- df[types %in% tx_types, , drop = FALSE]
  if (!nrow(sub)) return(NULL)

  gene_id <- first_present(sub, c("locus_tag", "gene_id", "gene", "Name"))
  # Aliases del transcrito, en orden de especificidad.
  alias_cols <- c("ID", "transcript_id", "Name", "protein_id", "orig_protein_id",
                  "Parent", "locus_tag")
  out <- list()
  for (col in alias_cols) {
    if (!col %in% names(sub)) next
    v <- sub[[col]]
    # `Parent` puede venir como CharacterList: se toma el primer padre.
    if (is.list(v)) v <- vapply(v, function(x) if (length(x)) as.character(x)[1] else NA_character_,
                                character(1))
    v <- as.character(v)
    ok <- !is.na(v) & nzchar(v) & !is.na(gene_id) & nzchar(gene_id)
    if (!any(ok)) next
    out[[col]] <- data.frame(TXNAME = v[ok], GENEID = gene_id[ok],
                             alias_type = col, stringsAsFactors = FALSE)
    # Variantes que producen los distintos pipelines: sin el prefijo "cds-"/
    # "rna-", y sin el sufijo de version.
    #
    # Las variantes se generan AQUI en lugar de delegar en
    # `tximport(ignoreTxVersion = TRUE)` porque ese flag corta por el primer
    # punto y no distingue una version de Ensembl de un accession de GenBank:
    # con IDs como "cds-AAC73112.1" (donde el ".1" es parte del accession)
    # descarta el 99 % de los transcritos en silencio. Generando los alias
    # explicitamente, `pick_tx2gene_for_quant()` elige el que de verdad casa y
    # tximport se llama con ignoreTxVersion = FALSE.
    variants <- list(
      sin_prefijo = sub("^(cds|rna|gene|id)-", "", v[ok]),
      sin_version = sub("[.][0-9]+$", "", v[ok]),
      sin_prefijo_ni_version = sub("[.][0-9]+$", "", sub("^(cds|rna|gene|id)-", "", v[ok]))
    )
    for (vn in names(variants)) {
      sv <- variants[[vn]]
      if (!any(sv != v[ok])) next
      out[[paste0(col, "_", vn)]] <- data.frame(
        TXNAME = sv, GENEID = gene_id[ok],
        alias_type = paste0(col, "_", vn), stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  res[!duplicated(res[, c("TXNAME", "GENEID", "alias_type")]), , drop = FALSE]
}

#' Elige el subconjunto de `tx2gene` que casa con los transcritos reales.
#'
#' Devuelve el mapa filtrado y la tasa de emparejamiento, para que la interfaz
#' pueda decir "8.412 de 8.500 transcritos mapeados (99 %)" o avisar de que la
#' anotacion no corresponde al transcriptoma usado.
pick_tx2gene_for_quant <- function(tx2gene, quant_tx_ids) {
  if (is.null(tx2gene) || !nrow(tx2gene) || !length(quant_tx_ids)) {
    return(list(tx2gene = NULL, rate = NA_real_, n_matched = 0L,
                n_tx = length(quant_tx_ids), alias_type = NA_character_))
  }
  ids <- unique(as.character(quant_tx_ids))
  # Se prueba cada convencion de alias por separado y gana la de mayor cobertura.
  best <- NULL
  for (at in unique(tx2gene$alias_type)) {
    sub <- tx2gene[tx2gene$alias_type == at, , drop = FALSE]
    n <- length(intersect(ids, sub$TXNAME))
    if (is.null(best) || n > best$n) best <- list(at = at, n = n, sub = sub)
  }
  if (is.null(best) || best$n == 0L) {
    return(list(tx2gene = NULL, rate = 0, n_matched = 0L, n_tx = length(ids),
                alias_type = NA_character_))
  }
  m <- best$sub[, c("TXNAME", "GENEID")]
  m <- m[!duplicated(m$TXNAME), , drop = FALSE]
  list(tx2gene = m, rate = best$n / length(ids), n_matched = best$n,
       n_tx = length(ids), alias_type = best$at)
}

#' Longitud de cada gen segun la anotacion, indexada por el identificador que usa
#' la matriz de conteos (locus_tag).
#'
#' Para genes con varios exones se suma la longitud de los exones, que es la
#' longitud efectiva que ve el contaje; para el resto se usa la del propio
#' feature. La necesita el diagnostico de sesgo de longitud (B3c).
gene_lengths_from_annotation <- function(path) {
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df) || !".width" %in% names(df)) return(NULL)
  types <- if ("type" %in% names(df)) as.character(df$type) else rep(NA_character_, nrow(df))
  gid <- first_present(df, c("locus_tag", "gene_id", "gene", "Name"))
  ok <- !is.na(gid) & nzchar(gid)
  if (!any(ok)) return(NULL)

  # Si hay exones, la longitud del gen es la suma de sus exones.
  ex <- ok & types %in% c("exon", "CDS")
  if (any(ex)) {
    len <- tapply(df$.width[ex], gid[ex], sum, na.rm = TRUE)
  } else {
    len <- tapply(df$.width[ok], gid[ok], max, na.rm = TRUE)
  }
  # Completa con la longitud del feature `gene` los que no tuvieran exones.
  g <- ok & types %in% c("gene", "pseudogene")
  if (any(g)) {
    gl <- tapply(df$.width[g], gid[g], max, na.rm = TRUE)
    missing <- setdiff(names(gl), names(len))
    if (length(missing)) len <- c(len, gl[missing])
  }
  out <- as.numeric(len)
  names(out) <- names(len)
  out[is.finite(out) & out > 0]
}

#' Detecta si la anotacion describe genes con splicing (varios exones por gen).
#'
#' Por que importa (docs/REVISION_ESTADISTICA.md, B10): `bowtie2` NO es
#' splice-aware. Sobre un procariota es correcto porque no hay intrones, pero
#' aplicar la misma ruta a un eucariota subestima sistematicamente los conteos en
#' las uniones exon-exon, y nada en la app lo advertia. Se detecta de la propia
#' anotacion en lugar de preguntar al usuario por el organismo: si los genes
#' tienen varios exones, la ruta de alineamiento con bowtie2 no es adecuada.
#'
#' @return list(spliced, n_genes_multiexon, n_genes_with_exons, frac, error)
detect_spliced_annotation <- function(path, min_frac = 0.1) {
  out <- list(spliced = NA, n_genes_multiexon = 0L, n_genes_with_exons = 0L,
              frac = NA_real_, error = NULL)
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df) || !"type" %in% names(df)) {
    out$error <- "No se ha podido leer la anotacion."
    return(out)
  }
  types <- as.character(df$type)
  ex <- types == "exon"
  if (!any(ex)) {
    # Sin features `exon` no se puede afirmar nada; en muchos GFF procariotas no
    # existen, lo que ya es un indicio de ausencia de splicing.
    out$spliced <- FALSE
    return(out)
  }
  parent <- df$Parent
  if (is.list(parent)) {
    parent <- vapply(parent, function(x) if (length(x)) as.character(x)[1] else NA_character_,
                     character(1))
  }
  parent <- as.character(parent)[ex]
  key <- if (all(is.na(parent))) first_present(df, c("locus_tag", "gene_id"))[ex] else parent
  key <- key[!is.na(key) & nzchar(key)]
  if (!length(key)) { out$spliced <- FALSE; return(out) }
  per_gene <- table(key)
  out$n_genes_with_exons <- length(per_gene)
  out$n_genes_multiexon  <- sum(per_gene > 1)
  out$frac <- out$n_genes_multiexon / out$n_genes_with_exons
  out$spliced <- isTRUE(out$frac >= min_frac)
  out
}

#' Identificadores de genes de rRNA segun la anotacion.
#'
#' Se recogen por dos vias (`type == "rRNA"` y `gene_biotype == "rRNA"`) y se
#' devuelven todos los alias, porque la matriz de conteos puede usar locus_tag o
#' nombre de gen segun como se generase.
rrna_ids_from_annotation <- function(path) {
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df)) return(character(0))
  is_rrna <- rep(FALSE, nrow(df))
  if ("type" %in% names(df)) {
    is_rrna <- is_rrna | as.character(df$type) %in% c("rRNA", "rRNA_gene")
  }
  for (col in c("gene_biotype", "biotype", "gbkey")) {
    if (col %in% names(df)) {
      is_rrna <- is_rrna | (!is.na(df[[col]]) & as.character(df[[col]]) == "rRNA")
    }
  }
  if (!any(is_rrna)) return(character(0))
  sub <- df[is_rrna, , drop = FALSE]
  ids <- character(0)
  for (col in c("locus_tag", "gene", "Name", "ID", "gene_id")) {
    if (col %in% names(sub)) {
      v <- sub[[col]]
      if (is.list(v)) v <- unlist(v, use.names = FALSE)
      v <- as.character(v)
      ids <- c(ids, v[!is.na(v) & nzchar(v)])
    }
  }
  unique(c(ids, sub("^(gene|rna)-", "", ids)))
}

# ── Traduccion de identificadores entre atributos de la anotacion ───────────
#
# El pipeline llama a `featureCounts -g locus_tag`, de modo que la matriz de
# conteos viene con locus tags (BW25113_RS00005). Ningun OrgDb los conoce: el
# enriquecimiento sobre la matriz tal cual mapea el 0 %. Y no se arregla
# pidiendole a featureCounts que agrupe por `gene`, porque ese atributo no esta
# en todos los registros del GTF y la herramienta aborta.
#
# La anotacion SI relaciona los dos: es la unica fuente que sabe que
# BW25113_RS00005 es thrL. Por eso la traduccion se hace con el GTF y no con el
# OrgDb, y por eso vive aqui y no en utils_enrich.R.
#
# Se traduce en el enriquecimiento, no al cargar la matriz: la tabla DEG, los
# informes y el script exportado conservan asi los identificadores que el
# pipeline produjo realmente, que es lo que hace el resultado reproducible.

#' Atributos del GFF/GTF que pueden identificar un gen, de mas a menos
#' especifico. `gene_id` primero porque es el que featureCounts escribe.
ANNOTATION_ID_ATTRS <- c("gene_id", "locus_tag", "old_locus_tag", "gene",
                         "protein_id", "transcript_id")

#' Atributos utiles como destino de una traduccion, con su etiqueta.
ANNOTATION_TARGET_ATTRS <- c(
  "gene"          = "Simbolo del gen (gene)",
  "locus_tag"     = "Locus tag (locus_tag)",
  "old_locus_tag" = "Locus tag antiguo (old_locus_tag)",
  "protein_id"    = "Identificador de proteina (protein_id)"
)

#' Que atributos de la anotacion estan realmente poblados.
#' Devuelve los nombres de columna con al menos un valor no vacio.
annotation_available_attrs <- function(path, candidatos = NULL) {
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df)) return(character(0))
  cand <- candidatos %||% union(ANNOTATION_ID_ATTRS, names(ANNOTATION_TARGET_ATTRS))
  cand <- intersect(cand, names(df))
  cand[vapply(cand, function(n) {
    v <- as.character(df[[n]])
    any(!is.na(v) & nzchar(v))
  }, logical(1))]
}

#' Deduce a que atributo de la anotacion corresponden unos identificadores.
#'
#' Compara la lista contra cada atributo candidato y se queda con el de mayor
#' cobertura. Se prefiere deducirlo a preguntarlo: el usuario no tiene por que
#' saber con que `-g` se lanzo featureCounts, y equivocarse produce una
#' traduccion silenciosamente vacia.
#'
#' @return list(attr, rate, tabla) o NULL si no se pudo leer la anotacion.
detect_annotation_keytype <- function(ids, path) {
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df)) return(NULL)
  ids <- unique(as.character(ids %||% character(0)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(NULL)

  cand <- intersect(ANNOTATION_ID_ATTRS, names(df))
  if (!length(cand)) return(NULL)
  tasas <- vapply(cand, function(n) {
    v <- as.character(df[[n]])
    mean(ids %in% v[!is.na(v) & nzchar(v)])
  }, numeric(1))
  mejor <- which.max(tasas)
  list(attr = cand[mejor], rate = unname(tasas[mejor]),
       tabla = data.frame(atributo = cand, cobertura = unname(tasas),
                          stringsAsFactors = FALSE))
}

#' Mapa de traduccion entre dos atributos de la anotacion.
#'
#' Cuando un mismo valor de origen apunta a varios de destino se queda con el
#' primero: en un GTF eso ocurre porque el gen tiene varias filas (gene, CDS,
#' exon) que repiten el mismo par, no porque haya ambiguedad real.
annotation_id_map <- function(path, from, to) {
  df <- read_annotation_features(path)
  if (is.null(df) || !nrow(df)) return(character(0))
  if (!all(c(from, to) %in% names(df))) return(character(0))
  o <- as.character(df[[from]]); d <- as.character(df[[to]])
  ok <- !is.na(o) & nzchar(o) & !is.na(d) & nzchar(d)
  if (!any(ok)) return(character(0))
  o <- o[ok]; d <- d[ok]
  dup <- duplicated(o)
  stats::setNames(d[!dup], o[!dup])
}

#' Traduce identificadores con la anotacion.
#'
#' La forma del valor devuelto imita la de `translate_to_entrez()` para que el
#' enriquecimiento pueda encadenar las dos traducciones sin casos especiales.
#'
#' @param from atributo de origen; si es NULL se deduce con detect_annotation_keytype().
#' @return list(ids, back, mapping, error, from). `back` va del identificador
#'   traducido al original, que es lo que permite devolver las tablas a los
#'   identificadores que el usuario reconoce.
translate_ids_with_annotation <- function(ids, path, from = NULL, to = "gene") {
  ids <- unique(as.character(ids %||% character(0)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  vacio <- function(msg, from_used = from) list(
    ids = character(0), back = character(0), error = msg, from = from_used,
    mapping = list(n_input = length(ids), n_mapped = 0L, rate = 0,
                   keytype = from_used %||% NA_character_, source = to))

  if (!length(ids)) return(vacio("Lista de genes vacia."))
  df <- read_annotation_features(path)
  if (is.null(df)) return(vacio("No se pudo leer la anotacion."))
  if (!to %in% names(df)) {
    return(vacio(paste0("La anotacion no trae el atributo '", to, "'.")))
  }
  if (is.null(from)) {
    det <- detect_annotation_keytype(ids, path)
    if (is.null(det) || det$rate == 0) {
      return(vacio(paste0("Ningun atributo de la anotacion reconoce estos ",
                          "identificadores.")))
    }
    from <- det$attr
  }
  if (identical(from, to)) {
    # Nada que traducir: se devuelve la identidad para que quien llame no tenga
    # que distinguir este caso.
    return(list(ids = ids, back = stats::setNames(ids, ids), error = NULL,
                from = from,
                mapping = list(n_input = length(ids), n_mapped = length(ids),
                               rate = 1, keytype = from, source = to)))
  }
  mapa <- annotation_id_map(path, from, to)
  if (!length(mapa)) return(vacio(paste0("No hay correspondencia entre '", from,
                                         "' y '", to, "' en la anotacion."), from))
  tr <- unname(mapa[ids])
  ok <- !is.na(tr) & nzchar(tr)
  if (!any(ok)) return(vacio(paste0("Ningun identificador se pudo traducir de '",
                                    from, "' a '", to, "'."), from))
  originales <- ids[ok]; traducidos <- tr[ok]
  # Varios identificadores de origen pueden caer en el mismo destino (en E. coli
  # son ~50 simbolos repartidos entre parologos y fragmentos anotados por
  # separado). Se conserva el primero y se cuenta cuantos se perdieron, porque
  # un enriquecimiento con genes duplicados infla el recuento de los terminos
  # que los contienen.
  dup <- duplicated(traducidos)
  list(
    ids = traducidos[!dup],
    back = stats::setNames(originales[!dup], traducidos[!dup]),
    error = NULL, from = from,
    mapping = list(n_input = length(ids), n_mapped = sum(!dup),
                   rate = sum(!dup) / length(ids), keytype = from, source = to,
                   n_colapsados = sum(dup))
  )
}

#' Traduce un vector con nombres (el ranking de GSEA) conservando los valores.
#'
#' Cuando dos identificadores de origen caen en el mismo destino se conserva el
#' de mayor valor absoluto. Es el criterio "max" que usa el GSEA original al
#' colapsar sondas a genes: quedarse con el primero seria arbitrario, y
#' promediar diluye la senal que el test busca.
translate_ranking_with_annotation <- function(ranked, path, from = NULL, to = "gene") {
  tr <- translate_ids_with_annotation(names(ranked), path, from = from, to = to)
  if (!length(tr$ids)) return(list(ranked = NULL, error = tr$error,
                                   mapping = tr$mapping, from = tr$from))
  mapa <- annotation_id_map(path, tr$from, to)
  destino <- unname(mapa[names(ranked)])
  ok <- !is.na(destino) & nzchar(destino)
  v <- ranked[ok]; destino <- destino[ok]
  orden <- order(abs(v), decreasing = TRUE)
  v <- v[orden]; destino <- destino[orden]
  dup <- duplicated(destino)
  out <- stats::setNames(unname(v[!dup]), destino[!dup])
  out <- sort(out, decreasing = TRUE)
  list(ranked = out, error = NULL, from = tr$from,
       mapping = list(n_input = length(ranked), n_mapped = length(out),
                      rate = length(out) / length(ranked),
                      keytype = tr$from, source = to,
                      n_colapsados = sum(dup)))
}
