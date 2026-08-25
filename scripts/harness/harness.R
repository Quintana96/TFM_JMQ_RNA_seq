#' scripts/harness/harness.R
#' Etapas del harness de validación.
#'
#' Reproduce lo que hace la aplicación sin abrir la interfaz. Llama a las MISMAS
#' funciones que el servidor de Shiny —`run_deg()`, `run_enrichment_*()`,
#' `run_gsea()`— y al mismo `workflow.sh`, de modo que lo que se valida es SARA
#' y no un script paralelo que podría diferir en un detalle.
#'
#' Que esto sea posible sin tocar la interfaz no es casualidad: la aplicación
#' separa las funciones puras del código que depende de `input`/`session`. El
#' harness consume esa capa pura directamente.
#'
#' Cada etapa escribe su resultado en `validacion/<dataset>/` y se salta si ya
#' está hecho, salvo que se pase `forzar = TRUE`. Así una ejecución interrumpida
#' se reanuda sin repetir lo caro.

RAIZ_VALIDACION <- "validacion"

#' Carga la aplicación completa en el entorno global.
cargar_app <- function() {
  suppressMessages({
    library(shiny); library(bslib); library(shinyjs); library(DT); library(shinyFiles)
  })
  source("global.R")
  for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
    try(source(f), silent = TRUE)
  }
  invisible(TRUE)
}

dir_dataset <- function(d) file.path(RAIZ_VALIDACION, d$id)

#' Directorio de salida de una ruta del pipeline.
#'
#' `pipeline_existente` permite reutilizar ejecuciones ya hechas en otro sitio,
#' que es lo habitual: procesar de nuevo cuatro rutas para probar el harness
#' costaría media hora y no aportaría nada.
dir_pipeline <- function(d, pipeline) {
  ya <- d$pipeline_existente[[pipeline]]
  if (!is.null(ya) && dir.exists(path.expand(ya))) return(path.expand(ya))
  file.path(dir_dataset(d), "pipeline", pipeline)
}
dir_etapa <- function(d, etapa) {
  p <- file.path(dir_dataset(d), etapa)
  dir.create(p, recursive = TRUE, showWarnings = FALSE); p
}
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

#' Samplesheet en el formato que espera `run_deg()`.
samplesheet <- function(d, muestras = NULL) {
  m <- muestras %||% d$muestras
  df <- data.frame(sample_id = names(m), condition = unname(m), stringsAsFactors = FALSE)
  if (!is.null(d$covariables)) df <- merge(df, d$covariables, by = "sample_id", all.x = TRUE)
  df[match(names(m), df$sample_id), , drop = FALSE]
}

#' Matriz de conteos de una ruta del pipeline, o la preexistente si el dataset
#' no tiene FASTQ que procesar.
matriz_de <- function(d, pipeline = NULL) {
  ruta <- if (is.null(pipeline)) d$conteos
          else file.path(dir_pipeline(d, pipeline), "04_counts", "count_matrix.tsv")
  if (!file.exists(ruta)) return(NULL)
  m <- utils::read.delim(ruta, check.names = FALSE)
  rn <- m[[1]]; m <- as.matrix(m[, -1, drop = FALSE]); rownames(m) <- rn
  storage.mode(m) <- "numeric"
  faltan <- setdiff(names(d$muestras), colnames(m))
  if (length(faltan)) {
    warning("Faltan muestras en ", ruta, ": ", paste(faltan, collapse = ", "))
    return(NULL)
  }
  round(m[, names(d$muestras), drop = FALSE])
}

# ── Etapa 1: pipeline ───────────────────────────────────────────────────────

#' Lanza `workflow.sh` para una ruta. La referencia depende de la ruta: las de
#' alineamiento van contra el genoma y las de pseudoalineamiento contra el
#' transcriptoma. Es el error más fácil de cometer a mano, porque el parámetro
#' del workflow se llama igual en los dos casos.
etapa_pipeline <- function(d, pipeline, hilos = 6, forzar = FALSE) {
  destino <- dir_pipeline(d, pipeline)
  estado <- file.path(destino, "exit_status.tsv")
  # Se exigen DOS cosas para saltar una etapa: que el estado diga success y que
  # la matriz de conteos exista. El estado por sí solo no basta —una ejecución
  # interrumpida podía quedar registrada como correcta— y la matriz es el
  # producto real de la etapa: si no está, la etapa no está hecha, diga lo que
  # diga el fichero de estado.
  matriz <- file.path(destino, "04_counts", "count_matrix.tsv")
  if (!forzar && file.exists(estado) &&
      identical(read_exit_status(destino)$status, "success") &&
      file.exists(matriz)) {
    msg("  ", pipeline, ": ya hecho, se salta"); return(invisible("saltado"))
  }
  referencia <- if (pipeline %in% c("bowtie2", "subjunc")) d$rutas$genoma
                else d$rutas$transcriptoma
  if (is.null(referencia) || !file.exists(referencia))
    stop("Falta la referencia para ", pipeline, ": ", referencia %||% "(sin definir)")
  dir.create(destino, recursive = TRUE, showWarnings = FALSE)
  msg("  ", pipeline, ": lanzando workflow.sh")
  args <- c("workflow.sh",
            "--INPUT", d$rutas$fastq, "--OUTPUT", destino,
            "--GENOME_FILE", referencia, "--ANNOTATION_FILE", d$rutas$anotacion,
            "--ALIGNMENT_TYPE", pipeline, "--READ_TYPE", d$read_type,
            "--THREADS", as.character(hilos))
  log_f <- file.path(destino, "harness.log")
  code <- system2("bash", args, stdout = log_f, stderr = log_f)
  st <- read_exit_status(destino)
  msg("  ", pipeline, ": codigo ", code, " | estado ", st$status %||% "?",
      " | ", fmt_duracion(st$duration_seconds), " | ", fmt_memoria(st$peak_rss_mb))
  invisible(st$status %||% "error")
}

# ── Etapa 2: expresión diferencial con los tres motores ─────────────────────

etapa_deg <- function(d, forzar = FALSE) {
  salida <- dir_etapa(d, "deg")
  entradas <- if (length(d$pipelines)) stats::setNames(as.list(d$pipelines), d$pipelines)
              else list(conteos = NULL)
  resultados <- list()
  for (nm in names(entradas)) {
    cm <- matriz_de(d, entradas[[nm]])
    if (is.null(cm)) { msg("  ", nm, ": sin matriz, se salta"); next }
    for (motor in d$motores) {
      etiqueta <- paste0(nm, "__", motor)
      f <- file.path(salida, paste0(etiqueta, ".rds"))
      if (!forzar && file.exists(f)) { resultados[[etiqueta]] <- readRDS(f); next }
      msg("  ", etiqueta, ": ajustando (", nrow(cm), " genes)")
      r <- run_deg(cm, samplesheet(d), method = motor,
                   ref_level = d$contraste$den, contrast_num = d$contraste$num,
                   batch = d$batch, fdr = 0.05, lfc_threshold = 0,
                   # Sin encogido: las tablas publicadas conservan la columna
                   # `stat`, que lfcShrink elimina, así que compararlas con un
                   # log2FC encogido mediría una diferencia de método.
                   shrink = FALSE)
      if (is.null(r$table)) { msg("    fallo: ", r$error %||% "?"); next }
      saveRDS(r, f); resultados[[etiqueta]] <- r
    }
  }
  resultados
}

# ── Etapa 3: enriquecimiento funcional ─────────────────────────────────────

#' Traduce los identificadores con la anotación si el dataset lo pide. Es lo que
#' evita que el enriquecimiento mapee el 0 % cuando la matriz trae locus tags.
traducir_si_toca <- function(d, genes) {
  e <- d$enriquecimiento
  if (is.null(e$traducir_con) || !file.exists(e$traducir_con)) return(list(ids = genes, tasa = NA))
  tr <- translate_ids_with_annotation(genes, e$traducir_con, to = e$traducir_a %||% "gene")
  if (!length(tr$ids)) return(list(ids = genes, tasa = 0))
  list(ids = tr$ids, tasa = tr$mapping$rate)
}

etapa_enriquecimiento <- function(d, degs, forzar = FALSE) {
  salida <- dir_etapa(d, "enriquecimiento")
  e <- d$enriquecimiento
  if (is.null(e)) return(list())
  # Se usa la ruta de referencia: la primera configurada, o la matriz suelta.
  ref <- if (length(d$pipelines)) paste0(d$pipelines[1], "__DESeq2") else "conteos__DESeq2"
  r <- degs[[ref]]
  if (is.null(r)) { msg("  sin resultado DEG de referencia (", ref, ")"); return(list()) }
  tab <- r$table
  sig_g <- tab$gene[!is.na(tab$padj) & tab$padj < 0.05 & abs(tab$log2FC) > 1]
  uni_g <- tab$gene[!is.na(tab$padj)]
  ts <- traducir_si_toca(d, sig_g); tu <- traducir_si_toca(d, uni_g)
  msg("  lista ", length(ts$ids), " | universo ", length(tu$ids),
      if (!is.na(ts$tasa)) sprintf(" | traduccion %.1f %%", 100*ts$tasa) else "")
  out <- list()
  for (col in e$colecciones) {
    f <- file.path(salida, paste0(col, ".rds"))
    if (!forzar && file.exists(f)) { out[[col]] <- readRDS(f); next }
    msg("  ", col, ": calculando ORA")
    res <- switch(col,
      "BP" = , "MF" = , "CC" = run_enrichment_go(ts$ids, universe = tu$ids,
                OrgDb = e$orgdb, ont = col, keyType = e$keytype),
      "KEGG" = run_enrichment_kegg(ts$ids, universe = tu$ids, organism = e$kegg_organismo,
                keyType = e$kegg_keytype %||% "kegg", OrgDb = e$orgdb,
                from_keytype = e$keytype),
      "REACTOME" = run_enrichment_reactome(ts$ids, universe = tu$ids, OrgDb = e$orgdb,
                keyType = e$keytype, organism = e$reactome_organismo %||% "human"),
      NULL)
    if (is.null(res)) next
    saveRDS(res, f); out[[col]] <- res
    msg("    ", col, ": ", if (is.null(res$table)) paste0("sin terminos (", res$error %||% "?", ")")
                           else paste0(nrow(res$table), " terminos"))
  }
  # GSEA sobre la colección principal, que no umbraliza y por tanto recupera la
  # señal direccional que el ORA sobre la lista completa cancela.
  f <- file.path(salida, "GSEA_BP.rds")
  if (forzar || !file.exists(f)) {
    rk <- deg_ranking_metric(tab, "stat")
    if (!is.null(rk$ranked)) {
      nombres <- traducir_si_toca(d, names(rk$ranked))
      if (length(nombres$ids) == length(rk$ranked)) names(rk$ranked) <- nombres$ids
      else {
        trk <- translate_ranking_with_annotation(rk$ranked, e$traducir_con,
                                                to = e$traducir_a %||% "gene")
        if (!is.null(trk$ranked)) rk$ranked <- trk$ranked
      }
      msg("  GSEA: ", length(rk$ranked), " genes en el ranking")
      g <- run_gsea(rk$ranked, ont = "BP", OrgDb = e$orgdb, keyType = e$keytype,
                    pvalueCutoff = 0.05)
      saveRDS(g, f); out$GSEA_BP <- g
      msg("    GSEA: ", if (is.null(g$table)) "sin conjuntos" else paste0(nrow(g$table), " conjuntos"))
    }
  } else out$GSEA_BP <- readRDS(f)
  out
}

# ── Etapa 4: métricas de validación ─────────────────────────────────────────

sig_de <- function(tab, fdr = 0.05, lfc = 1)
  tab$gene[!is.na(tab$padj) & tab$padj < fdr & abs(tab$log2FC) > lfc]

#' Tabla publicada, con los identificadores ya alineados a los nuestros.
tabla_publicada <- function(d) {
  p <- d$publicado
  if (is.null(p) || !file.exists(p$tabla)) return(NULL)
  t <- utils::read.delim(p$tabla, sep = p$sep %||% "\t", stringsAsFactors = FALSE)
  t <- t[!is.na(t[[p$col_id]]) & nzchar(t[[p$col_id]]), ]
  t <- t[!duplicated(t[[p$col_id]]), ]
  data.frame(id = t[[p$col_id]], lfc = t[[p$col_lfc]], padj = t[[p$col_padj]],
             row.names = t[[p$col_id]], stringsAsFactors = FALSE)
}

#' Concordancia de un resultado propio con la tabla publicada.
metricas_vs_publicado <- function(d, r) {
  pub <- tabla_publicada(d); if (is.null(pub)) return(NULL)
  p <- d$publicado; fdr <- p$fdr %||% 0.05; lfc <- p$lfc %||% 1
  t <- r$table; rownames(t) <- t$gene
  com <- intersect(rownames(t), rownames(pub))
  if (length(com) < 50) return(NULL)
  a <- t[com, "log2FC"]; b <- pub[com, "lfc"]; ok <- is.finite(a) & is.finite(b)
  sp <- rownames(pub)[!is.na(pub$padj) & pub$padj < fdr & abs(pub$lfc) > lfc]
  sp <- intersect(sp, com); sn <- intersect(sig_de(t, fdr, lfc), com)
  inter <- intersect(sp, sn)
  # La pendiente por el origen detecta un sesgo sistemático en la magnitud que
  # la correlación no ve: dos series perfectamente correlacionadas pueden
  # diferir en escala.
  pend <- tryCatch(unname(coef(stats::lm(a[ok] ~ 0 + b[ok]))[1]), error = function(e) NA_real_)
  list(genes_comunes = length(com),
       deg_propios = length(sn), deg_publicados = length(sp), coincidentes = length(inter),
       recuperados = if (length(sp)) length(inter)/length(sp) else NA_real_,
       jaccard = if (length(union(sp, sn))) length(inter)/length(union(sp, sn)) else NA_real_,
       pearson = if (sum(ok) > 2) stats::cor(a[ok], b[ok]) else NA_real_,
       spearman = if (sum(ok) > 2) stats::cor(a[ok], b[ok], method = "spearman") else NA_real_,
       pendiente = pend,
       dif_mediana = stats::median(abs(a[ok] - b[ok])),
       mismo_signo = if (length(inter))
         mean(sign(t[inter,"log2FC"]) == sign(pub[inter,"lfc"])) else NA_real_)
}

#' Acuerdo entre dos resultados propios (dos motores, o dos rutas).
metricas_entre <- function(rA, rB) {
  A <- rA$table; rownames(A) <- A$gene; B <- rB$table; rownames(B) <- B$gene
  com <- intersect(rownames(A), rownames(B))
  if (length(com) < 50) return(NULL)
  a <- A[com,"log2FC"]; b <- B[com,"log2FC"]; ok <- is.finite(a)&is.finite(b)
  sA <- intersect(sig_de(A), com); sB <- intersect(sig_de(B), com)
  list(genes_comunes = length(com), deg_A = length(sA), deg_B = length(sB),
       coincidentes = length(intersect(sA,sB)),
       jaccard = if (length(union(sA,sB))) length(intersect(sA,sB))/length(union(sA,sB)) else NA_real_,
       pearson = if (sum(ok)>2) stats::cor(a[ok],b[ok]) else NA_real_,
       mismo_signo = mean(sign(a[ok])==sign(b[ok])))
}

#' Control de permutación: baraja las etiquetas y mide cuántos genes salen.
#'
#' Es la única métrica que mide ESPECIFICIDAD. Las demás comparan con lo
#' publicado, y por tanto no pueden distinguir un pipeline correcto de uno que
#' inventa señal de forma reproducible.
etapa_permutacion <- function(d, forzar = FALSE) {
  salida <- dir_etapa(d, "metricas")
  f <- file.path(salida, "permutacion.rds")
  if (!forzar && file.exists(f)) return(readRDS(f))
  cm <- matriz_de(d, if (length(d$pipelines)) d$pipelines[1] else NULL)
  if (is.null(cm)) return(NULL)
  grupos <- d$muestras; niveles <- unique(unname(grupos))
  if (length(niveles) != 2) return(NULL)
  n1 <- sum(grupos == niveles[1])
  todas <- utils::combn(names(grupos), n1, simplify = FALSE)
  reales <- names(grupos)[grupos == niveles[1]]
  perms <- Filter(function(p) !setequal(p, reales) &&
                    !setequal(p, setdiff(names(grupos), reales)), todas)
  msg("  permutacion: ", length(perms), " reetiquetados posibles")
  n_perm <- integer(0)
  for (p in perms) {
    et <- stats::setNames(ifelse(names(grupos) %in% p, niveles[1], niveles[2]), names(grupos))
    r <- tryCatch(run_deg(cm, samplesheet(d, et), method = "DESeq2",
                          ref_level = niveles[2], contrast_num = niveles[1],
                          fdr = 0.05, lfc_threshold = 0, shrink = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || is.null(r$table)) next
    n_perm <- c(n_perm, length(sig_de(r$table)))
  }
  real <- run_deg(cm, samplesheet(d), method = "DESeq2", ref_level = d$contraste$den,
                  contrast_num = d$contraste$num, fdr = 0.05, lfc_threshold = 0, shrink = FALSE)
  res <- list(deg_real = length(sig_de(real$table)), permutaciones = n_perm,
              mediana = stats::median(n_perm), maximo = max(n_perm),
              ceros = sum(n_perm == 0), n = length(n_perm))
  saveRDS(res, f)
  msg("  permutacion: real ", res$deg_real, " | barajado mediana ", res$mediana,
      " (max ", res$maximo, "; ", res$ceros, "/", res$n, " a cero)")
  res
}

#' Descomposición del error: el mismo análisis sobre los conteos publicados.
#' Separa lo que aporta la etapa estadística de lo que aporta el pipeline.
etapa_descomposicion <- function(d, forzar = FALSE) {
  salida <- dir_etapa(d, "metricas")
  f <- file.path(salida, "descomposicion.rds")
  if (!forzar && file.exists(f)) return(readRDS(f))
  p <- d$publicado
  if (is.null(p) || is.null(p$conteos) || !file.exists(path.expand(p$conteos))) return(NULL)
  m <- utils::read.delim(path.expand(p$conteos), check.names = FALSE)
  rn <- m[[1]]; m <- round(as.matrix(m[, -1, drop = FALSE])); rownames(m) <- rn
  if (!is.null(p$muestras_publicadas))
    colnames(m) <- names(p$muestras_publicadas)[match(colnames(m), p$muestras_publicadas)]
  faltan <- setdiff(names(d$muestras), colnames(m))
  if (length(faltan)) return(NULL)
  m <- m[, names(d$muestras), drop = FALSE]
  r <- run_deg(m, samplesheet(d), method = "DESeq2", ref_level = d$contraste$den,
               contrast_num = d$contraste$num, fdr = 0.05, lfc_threshold = 0, shrink = FALSE)
  if (is.null(r$table)) return(NULL)
  # Si la matriz publicada viene con símbolos y la tabla publicada trae ambos,
  # se pasa a los identificadores de la tabla para poder cruzar.
  if (!is.null(p$mapa_desde_tabla)) {
    bruta <- utils::read.delim(p$tabla, sep = p$sep %||% "\t", stringsAsFactors = FALSE)
    mapa <- stats::setNames(bruta[[p$mapa_desde_tabla[["a"]]]],
                            bruta[[p$mapa_desde_tabla[["de"]]]])
    r$table$gene <- unname(mapa[r$table$gene])
    r$table <- r$table[!is.na(r$table$gene) & !duplicated(r$table$gene), ]
  }
  res <- metricas_vs_publicado(d, r)
  saveRDS(res, f); res
}
