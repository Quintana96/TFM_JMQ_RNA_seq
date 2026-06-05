#' utils_enrich.R
#' Funciones puras para enriquecimiento funcional (GO / KEGG) via clusterProfiler.
#' No dependen de Shiny. Tanto si los paquetes no estan instalados como si la
#' query falla, devuelven list(table = NULL, error = mensaje).

#' Enriquecimiento GO. Si OrgDb es NULL, devuelve mensaje pidiendo OrgDb.
run_enrichment_go <- function(genes, universe = NULL, OrgDb = NULL,
                              ont = "BP", keyType = "SYMBOL",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(list(table = NULL, error = "clusterProfiler no esta instalado."))
  }
  if (is.null(OrgDb) || (is.character(OrgDb) && !nzchar(OrgDb))) {
    return(list(table = NULL,
                error = "Especifica un OrgDb (p.ej. 'org.EcK12.eg.db') para correr GO."))
  }
  if (is.character(OrgDb)) {
    if (!requireNamespace(OrgDb, quietly = TRUE)) {
      return(list(table = NULL,
                  error = paste0("El paquete '", OrgDb, "' no esta instalado.")))
    }
    OrgDb <- getFromNamespace(OrgDb, OrgDb)
  }
  if (!length(genes)) {
    return(list(table = NULL, error = "Lista de genes vacia."))
  }
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
    if (is.null(ego) || is.null(ego@result) || !nrow(ego@result)) return(NULL)
    as.data.frame(ego)
  }, error = function(e) e)
  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out)))
  if (is.null(out)) return(list(table = NULL, error = "Sin terminos GO enriquecidos."))
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], error = NULL)
}

#' Enriquecimiento KEGG. Para E. coli K12 substr MG1655 usa organism = "eco".
run_enrichment_kegg <- function(genes, organism = "eco", keyType = "kegg",
                                pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(list(table = NULL, error = "clusterProfiler no esta instalado."))
  }
  if (!length(genes)) return(list(table = NULL, error = "Lista de genes vacia."))
  out <- tryCatch({
    ek <- clusterProfiler::enrichKEGG(
      gene          = unique(as.character(genes)),
      organism      = organism,
      keyType       = keyType,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff
    )
    if (is.null(ek) || is.null(ek@result) || !nrow(ek@result)) return(NULL)
    as.data.frame(ek)
  }, error = function(e) e)
  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out)))
  if (is.null(out)) return(list(table = NULL, error = "Sin terminos KEGG enriquecidos."))
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], error = NULL)
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
