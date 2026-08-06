#' utils_enrich.R
#' Funciones puras para enriquecimiento funcional (GO / KEGG) via clusterProfiler.
#' No dependen de Shiny. Tanto si los paquetes no estan instalados como si la
#' query falla, devuelven list(table = NULL, error = mensaje, mapping = ...).
#'
#' Dos requisitos que Wijesooriya et al. (PLoS Comput Biol 2022) encontraron
#' incumplidos en el 95 % de los analisis de sobre-representacion publicados, y
#' que aqui se respetan explicitamente (ver docs/REVISION_ESTADISTICA.md, B3):
#'
#'   1. LISTA DE FONDO. Tanto GO como KEGG reciben `universe` = los genes
#'      efectivamente testeados. Usar todo el genoma como fondo infla los
#'      enriquecimientos.
#'   2. TASA DE MAPEO VISIBLE. `mapping` viaja siempre en el resultado, incluso
#'      cuando no hay terminos, para que la interfaz pueda distinguir "no hay
#'      enriquecimiento" de "solo mapeo el 12 % de los IDs".

#' Keytypes disponibles en un OrgDb (para poblar el selector de la interfaz).
#' Devuelve character(0) si no se puede consultar.
orgdb_keytypes <- function(OrgDb) {
  if (is.character(OrgDb)) {
    if (!nzchar(OrgDb) || !requireNamespace(OrgDb, quietly = TRUE)) return(character(0))
    OrgDb <- getFromNamespace(OrgDb, OrgDb)
  }
  if (is.null(OrgDb) || !requireNamespace("AnnotationDbi", quietly = TRUE)) return(character(0))
  tryCatch(AnnotationDbi::keytypes(OrgDb), error = function(e) character(0))
}

#' Tasa de mapeo de una lista de genes contra un OrgDb y un keyType.
#'
#' Un enriquecimiento con mapeo bajo no es interpretable. Los IDs que salen de
#' featureCounts sobre un GFF de E. coli son locus tags (b0001) o IDs de Ensembl
#' Bacteria, no simbolos, asi que con keyType = "SYMBOL" el mapeo puede ser
#' casi nulo sin que nada falle visiblemente.
gene_mapping_rate <- function(genes, OrgDb, keyType = "SYMBOL") {
  genes <- unique(as.character(genes[!is.na(genes)]))
  out <- list(n_input = length(genes), n_mapped = NA_integer_, rate = NA_real_,
              keytype = keyType, source = "OrgDb")
  if (!length(genes) || is.null(OrgDb)) return(out)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(out)
  ks <- tryCatch(AnnotationDbi::keys(OrgDb, keytype = keyType), error = function(e) NULL)
  if (is.null(ks)) return(out)
  out$n_mapped <- length(intersect(genes, ks))
  out$rate <- out$n_mapped / length(genes)
  out
}

#' Tasa de mapeo derivada del denominador de GeneRatio ("k/n"): `n` es el numero
#' de genes de entrada que clusterProfiler ha podido anotar. Es la via para
#' KEGG, donde no hay una base local contra la que comparar.
mapping_from_generatio <- function(enrich_df, n_input, keyType = NA_character_) {
  out <- list(n_input = n_input, n_mapped = NA_integer_, rate = NA_real_,
              keytype = keyType, source = "GeneRatio")
  if (is.null(enrich_df) || !nrow(enrich_df) || !"GeneRatio" %in% names(enrich_df)) return(out)
  den <- suppressWarnings(as.integer(sub(".*/", "", as.character(enrich_df$GeneRatio[1]))))
  if (is.na(den)) return(out)
  out$n_mapped <- den
  out$rate <- if (n_input > 0) den / n_input else NA_real_
  out
}

#' Texto corto con la tasa de mapeo, listo para mostrar bajo el grafico.
mapping_rate_text <- function(mapping) {
  if (is.null(mapping) || is.null(mapping$n_input)) return(NULL)
  if (is.na(mapping$n_mapped %||% NA)) {
    return(paste0(mapping$n_input, " genes de entrada; tasa de mapeo no determinable."))
  }
  paste0(
    fmt_int(mapping$n_mapped), " de ", fmt_int(mapping$n_input), " genes mapeados (",
    round(100 * mapping$rate, 1), " %)",
    if (!is.na(mapping$keytype)) paste0(" con keyType = ", mapping$keytype) else ""
  )
}

#' Enriquecimiento GO. Si OrgDb es NULL, devuelve mensaje pidiendo OrgDb.
#'
#' @param simplify_terms Si TRUE, colapsa terminos GO redundantes con
#'   `clusterProfiler::simplify()`. GO es jerarquico y un termino padre solapa
#'   con sus hijos, asi que las listas salen dominadas por variaciones del mismo
#'   concepto; simplify usa similitud semantica (GOSemSim, medida de Wang) para
#'   quedarse con un representante por grupo.
#' @param simplify_cutoff Umbral de similitud por encima del cual dos terminos se
#'   consideran redundantes.
run_enrichment_go <- function(genes, universe = NULL, OrgDb = NULL,
                              ont = "BP", keyType = "SYMBOL",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              simplify_terms = FALSE, simplify_cutoff = 0.7) {
  fail <- function(msg, mapping = NULL) list(table = NULL, error = msg, mapping = mapping)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  if (is.null(OrgDb) || (is.character(OrgDb) && !nzchar(OrgDb))) {
    return(fail("Especifica un OrgDb (p.ej. 'org.EcK12.eg.db') para correr GO."))
  }
  if (is.character(OrgDb)) {
    if (!requireNamespace(OrgDb, quietly = TRUE)) {
      return(fail(paste0("El paquete '", OrgDb, "' no esta instalado.")))
    }
    OrgDb <- getFromNamespace(OrgDb, OrgDb)
  }
  if (!length(genes)) return(fail("Lista de genes vacia."))

  # Se calcula antes del enriquecimiento para poder explicar un resultado vacio.
  mapping <- gene_mapping_rate(genes, OrgDb, keyType)

  out <- tryCatch({
    ego <- clusterProfiler::enrichGO(
      gene          = unique(as.character(genes)),
      universe      = if (length(universe)) unique(as.character(universe)) else NULL,
      OrgDb         = OrgDb,
      keyType       = keyType,
      ont           = ont,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      readable      = FALSE
    )
    # simplify() necesita el objeto S4, asi que se aplica antes de convertir.
    # Solo tiene sentido en las tres ontologias concretas (no en "ALL"), porque
    # la similitud semantica se calcula dentro de un mismo grafo.
    if (isTRUE(simplify_terms) && !is.null(ego) && ont %in% c("BP", "MF", "CC") &&
        requireNamespace("GOSemSim", quietly = TRUE)) {
      ego <- tryCatch(
        clusterProfiler::simplify(ego, cutoff = simplify_cutoff,
                                  by = "p.adjust", measure = "Wang"),
        error = function(e) ego
      )
    }
    # Sin `return()`: dentro de tryCatch({...}) un return() sale de la funcion
    # entera y se lleva por delante el resultado (incluido `mapping`).
    # Se cuenta sobre as.data.frame() y no sobre @result porque @result guarda
    # todos los terminos testeados y la conversion aplica los cutoffs.
    df <- if (is.null(ego)) NULL else as.data.frame(ego)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e)
  if (inherits(out, "error")) return(fail(conditionMessage(out), mapping))
  if (is.null(out)) return(fail("Sin terminos GO enriquecidos.", mapping))
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], error = NULL, mapping = mapping)
}

#' Enriquecimiento KEGG. Para E. coli K12 substr MG1655 usa organism = "eco".
#'
#' `universe` es el cambio importante respecto a la version anterior: sin el,
#' enrichKEGG usa todo el genoma como fondo en lugar de los genes testeados.
run_enrichment_kegg <- function(genes, universe = NULL, organism = "eco",
                                keyType = "kegg",
                                pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  fail <- function(msg, mapping = NULL) list(table = NULL, error = msg, mapping = mapping)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  if (!length(genes)) return(fail("Lista de genes vacia."))
  genes_u <- unique(as.character(genes))
  out <- tryCatch({
    ek <- clusterProfiler::enrichKEGG(
      gene          = genes_u,
      universe      = if (length(universe)) unique(as.character(universe)) else NULL,
      organism      = organism,
      keyType       = keyType,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff
    )
    df <- if (is.null(ek)) NULL else as.data.frame(ek)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e)
  if (inherits(out, "error")) {
    return(fail(conditionMessage(out),
                mapping_from_generatio(NULL, length(genes_u), keyType)))
  }
  if (is.null(out)) {
    return(fail("Sin terminos KEGG enriquecidos.",
                mapping_from_generatio(NULL, length(genes_u), keyType)))
  }
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], error = NULL,
       mapping = mapping_from_generatio(out, length(genes_u), keyType))
}

# ── GSEA ────────────────────────────────────────────────────────────────────

#' Construye el vector ordenado de genes para GSEA.
#'
#' La eleccion de metrica es la decision de diseno real:
#'   - `stat` (por defecto): el estadistico moderado del motor. Incorpora la
#'     variabilidad y es la opcion habitual.
#'   - `log2FC`: interpretable, pero ignora la varianza.
#'   - `signed_p`: `sign(log2FC) * -log10(pvalue)`. Es popular pero genera
#'     EMPATES cuando los p-valores saturan, y GSEA no los resuelve: el orden
#'     dentro de un grupo empatado queda arbitrario. Por eso se devuelve
#'     `ties_frac` y la interfaz avisa.
deg_ranking_metric <- function(deg_df, metric = c("stat", "log2FC", "signed_p")) {
  metric <- match.arg(metric)
  if (is.null(deg_df) || !nrow(deg_df)) {
    return(list(ranked = NULL, error = "Sin tabla DEG.", metric = metric))
  }
  v <- switch(metric,
    "stat"     = deg_df$stat,
    "log2FC"   = deg_df$log2FC,
    "signed_p" = sign(deg_df$log2FC) *
                   -log10(pmax(deg_df$pvalue, .Machine$double.xmin))
  )
  if (is.null(v) || all(is.na(v))) {
    return(list(ranked = NULL, metric = metric,
                error = paste0("La metrica '", metric,
                               "' no esta disponible para este motor.")))
  }
  ok <- !is.na(v) & is.finite(v) & !is.na(deg_df$gene)
  v <- v[ok]
  names(v) <- deg_df$gene[ok]
  v <- v[!duplicated(names(v))]
  v <- sort(v, decreasing = TRUE)
  # Fraccion de genes que comparten su valor con al menos otro gen.
  ties_frac <- if (length(v)) sum(duplicated(v) | duplicated(v, fromLast = TRUE)) / length(v) else NA_real_
  list(ranked = v, metric = metric, n = length(v), ties_frac = ties_frac, error = NULL)
}

#' GSEA sobre el ranking completo, via clusterProfiler (backend fgsea).
#'
#' `exponent = 0` es la permutacion clasica NO ponderada. El benchmark con datos
#' curados de 12 tipos de cancer encontro que es la que mejor equilibra
#' sensibilidad y especificidad, y que los estadisticos ponderados (p = 1; 1,5; 2)
#' DEGRADAN el rendimiento: su justificacion original venia de microarrays y no
#' se traslada a RNA-seq. En fgsea el parametro equivalente es `gseaParam`.
#'
#' A diferencia del ORA, GSEA parte del ranking completo sin umbralizar, asi que
#' ve senales debiles pero coordinadas: si veinte genes de una ruta suben un 30 %
#' cada uno, ninguno pasa el corte de la lista y la ruta no aparece en el ORA.
run_gsea <- function(ranked, ont = "BP", OrgDb = NULL, organism = "eco",
                     keyType = "SYMBOL", exponent = 0,
                     pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500) {
  fail <- function(msg) list(table = NULL, error = msg, mapping = NULL)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    return(fail("fgsea no esta instalado."))
  }
  if (is.null(ranked) || !length(ranked)) return(fail("Ranking de genes vacio."))

  is_kegg <- identical(ont, "KEGG")
  if (!is_kegg) {
    if (is.null(OrgDb) || (is.character(OrgDb) && !nzchar(OrgDb))) {
      return(fail("Especifica un OrgDb para correr GSEA sobre GO."))
    }
    if (is.character(OrgDb)) {
      if (!requireNamespace(OrgDb, quietly = TRUE)) {
        return(fail(paste0("El paquete '", OrgDb, "' no esta instalado.")))
      }
      OrgDb <- getFromNamespace(OrgDb, OrgDb)
    }
  }
  mapping <- if (!is_kegg) gene_mapping_rate(names(ranked), OrgDb, keyType) else NULL

  out <- tryCatch({
    gs <- if (is_kegg) {
      clusterProfiler::gseKEGG(
        geneList = ranked, organism = organism, keyType = keyType,
        exponent = exponent, minGSSize = minGSSize, maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff, verbose = FALSE, seed = TRUE
      )
    } else {
      clusterProfiler::gseGO(
        geneList = ranked, ont = ont, OrgDb = OrgDb, keyType = keyType,
        exponent = exponent, minGSSize = minGSSize, maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff, verbose = FALSE, seed = TRUE
      )
    }
    df <- if (is.null(gs)) NULL else as.data.frame(gs)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e)

  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out),
                                          mapping = mapping))
  if (is.null(out)) return(list(table = NULL, error = "Sin conjuntos de genes enriquecidos.",
                                mapping = mapping))
  keep <- intersect(c("ID", "Description", "setSize", "enrichmentScore", "NES",
                      "pvalue", "p.adjust", "qvalue", "core_enrichment"), names(out))
  list(table = out[, keep, drop = FALSE], error = NULL, mapping = mapping)
}

#' Ordena el enriquecimiento por p.adjust y devuelve top_n filas.
enrichment_dotplot_data <- function(enrich_df, top_n = 15) {
  if (is.null(enrich_df) || !nrow(enrich_df)) return(NULL)
  ord <- order(enrich_df$p.adjust, na.last = TRUE)
  out <- enrich_df[ord, , drop = FALSE]
  if (nrow(out) > top_n) out <- out[seq_len(top_n), , drop = FALSE]
  # Convertir GeneRatio "x/y" a numerico
  if ("GeneRatio" %in% names(out)) {
    out$GeneRatioNum <- sapply(strsplit(as.character(out$GeneRatio), "/", fixed = TRUE),
                               function(p) if (length(p) == 2L) suppressWarnings(as.numeric(p[1]) / as.numeric(p[2])) else NA_real_)
  }
  out
}
