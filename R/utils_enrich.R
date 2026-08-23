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

#' Version y fechas de las fuentes de un OrgDb.
#'
#' Los resultados de enriquecimiento cambian con la version de la anotacion:
#' Wadi et al. (Nature Methods 2016) mostraron que usar anotaciones
#' desactualizadas altera sustancialmente las rutas que salen enriquecidas. Por
#' eso el informe declara con que version se ejecuto, igual que declara la
#' version de los paquetes.
#'
#' Nota practica: el campo KEGG de los OrgDb esta congelado desde 2011, asi que
#' su fecha NO describe la anotacion KEGG que usa la app (que consulta la API en
#' linea); se muestra igualmente para que quede claro que no se usa esa via.
#'
#' @param OrgDb nombre del paquete OrgDb, o NULL
#' @return lista clave-valor con la version y las fechas disponibles
orgdb_source_info <- function(OrgDb) {
  if (is.null(OrgDb) || !nzchar(OrgDb %||% "") ||
      !requireNamespace(OrgDb, quietly = TRUE)) {
    return(list("OrgDb" = "no disponible"))
  }
  out <- list("OrgDb" = paste0(OrgDb, " ", as.character(utils::packageVersion(OrgDb))))
  md <- tryCatch({
    obj <- getFromNamespace(OrgDb, OrgDb)
    AnnotationDbi::metadata(obj)
  }, error = function(e) NULL)
  if (is.data.frame(md) && all(c("name", "value") %in% names(md))) {
    interes <- c("GOSOURCEDATE", "GOEGSOURCEDATE", "EGSOURCEDATE", "ENSOURCEDATE",
                 "KEGGSOURCEDATE", "GOSOURCEVERSION", "DBSCHEMAVERSION")
    hit <- md[md$name %in% interes, , drop = FALSE]
    for (i in seq_len(nrow(hit))) out[[hit$name[i]]] <- hit$value[i]
  }
  out
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
  # Admite tanto el objeto como el nombre del paquete: llamada con el nombre,
  # AnnotationDbi::keys() falla y la tasa se perdia como NA sin explicacion.
  OrgDb <- as_orgdb_object(OrgDb)
  if (is.null(OrgDb)) return(out)
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

#' Convierte los IDs crudos de un resultado de clusterProfiler a simbolos.
#'
#' `setReadable()` reescribe las columnas geneID (ORA) y core_enrichment (GSEA),
#' que son las unicas que el usuario lee gen a gen. Un listado de ENSG00000...
#' obliga a salir de la app para saber de que gen se habla; el mismo listado en
#' simbolos se interpreta directamente. Es cosmetico y puede fallar (keyType que
#' setReadable no reconoce), asi que ante cualquier error se devuelve el objeto
#' original sin tocar en lugar de tumbar el enriquecimiento entero.
enrich_set_readable <- function(obj, OrgDb, keyType) {
  if (is.null(obj) || is.null(OrgDb)) return(obj)
  if (identical(keyType, "SYMBOL")) return(obj)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(obj)
  tryCatch(clusterProfiler::setReadable(obj, OrgDb = OrgDb, keyType = keyType),
           error = function(e) obj)
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
#' @param readable Si TRUE, la columna geneID sale con simbolos en vez de con
#'   los IDs de entrada.
run_enrichment_go <- function(genes, universe = NULL, OrgDb = NULL,
                              ont = "BP", keyType = "SYMBOL",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              simplify_terms = FALSE, simplify_cutoff = 0.7,
                              readable = FALSE) {
  obj_s4 <- NULL
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
    if (isTRUE(readable)) ego <- enrich_set_readable(ego, OrgDb, keyType)
    # Sin `return()`: dentro de tryCatch({...}) un return() sale de la funcion
    # entera y se lleva por delante el resultado (incluido `mapping`).
    # Se cuenta sobre as.data.frame() y no sobre @result porque @result guarda
    # todos los terminos testeados y la conversion aplica los cutoffs.
    obj_s4 <- ego   # el objeto S4 lo necesitan los graficos de enrichplot
    df <- if (is.null(ego)) NULL else as.data.frame(ego)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e)
  if (inherits(out, "error")) return(fail(conditionMessage(out), mapping))
  if (is.null(out)) return(fail("Sin terminos GO enriquecidos.", mapping))
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], obj = obj_s4, error = NULL, mapping = mapping)
}

#' Enriquecimiento KEGG. Para E. coli K12 substr MG1655 usa organism = "eco".
#'
#' `universe` es el cambio importante respecto a la version anterior: sin el,
#' enrichKEGG usa todo el genoma como fondo en lugar de los genes testeados.
#' `from_keytype` es el espacio de identificadores de la matriz de conteos, y
#' `OrgDb` la anotacion con la que traducirlo. Sin ellos, pedir "ncbi-geneid" con
#' una matriz de simbolos mapeaba el 0 %% y devolvia "sin terminos enriquecidos",
#' que no distingue "no hay senal" de "los IDs no eran los que KEGG esperaba".
#' Reactome ya traducia; KEGG no, y era la unica coleccion que no lo hacia.
run_enrichment_kegg <- function(genes, universe = NULL, organism = "eco",
                                keyType = "kegg", OrgDb = NULL,
                                from_keytype = NULL,
                                pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  obj_s4 <- NULL
  fail <- function(msg, mapping = NULL) list(table = NULL, error = msg, mapping = mapping)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  if (!length(genes)) return(fail("Lista de genes vacia."))
  genes_u <- unique(as.character(genes))

  # Traduccion a ENTREZID cuando se pide "ncbi-geneid" y los genes vienen en otro
  # espacio. El universo se traduce con las mismas reglas: un fondo en un espacio
  # y una lista en otro darian un enriquecimiento sin sentido.
  traduccion <- NULL
  if (identical(keyType, "ncbi-geneid") && !is.null(from_keytype) &&
      !identical(from_keytype, "ENTREZID") && !is.null(OrgDb)) {
    tr <- translate_to_entrez(genes_u, OrgDb, keyType = from_keytype)
    if (!length(tr$ids)) {
      return(fail(paste0("No se pudo traducir ningun gen de '", from_keytype,
                         "' a ENTREZID: ", tr$error %||% "sin coincidencias."),
                  tr$mapping))
    }
    traduccion <- tr
    genes_u <- tr$ids
    if (length(universe)) {
      tu <- translate_to_entrez(unique(as.character(universe)), OrgDb,
                                keyType = from_keytype)
      universe <- if (length(tu$ids)) tu$ids else NULL
    }
  }

  out <- tryCatch({
    ek <- clusterProfiler::enrichKEGG(
      gene          = genes_u,
      universe      = if (length(universe)) unique(as.character(universe)) else NULL,
      organism      = organism,
      keyType       = keyType,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff
    )
    obj_s4 <- ek   # el objeto S4 lo necesitan los graficos de enrichplot
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
  tab <- out[, keep, drop = FALSE]
  # Con traduccion, la columna geneID trae ENTREZIDs, que no dicen nada al leer
  # la tabla. Se devuelven a los identificadores de entrada.
  if (!is.null(traduccion) && "geneID" %in% names(tab)) {
    tab$geneID <- vapply(strsplit(as.character(tab$geneID), "/"), function(ids) {
      v <- traduccion$back[ids]
      paste(ifelse(is.na(v), ids, v), collapse = "/")
    }, character(1))
  }
  list(table = tab, obj = obj_s4, error = NULL,
       # Con traduccion la tasa real es la de la traduccion: la deducida del
       # GeneRatio mide otra cosa (cuantos de los traducidos estan en KEGG).
       mapping = if (!is.null(traduccion)) traduccion$mapping
                 else mapping_from_generatio(out, length(genes_u), keyType))
}

# ── Gene sets propios (GMT) ─────────────────────────────────────────────────
#
# Un GMT desacopla el enriquecimiento del OrgDb, y eso resuelve dos casos que
# GO y KEGG dejan fuera:
#   - Organismos sin anotacion rica en Bioconductor. Para la mayoria de
#     procariotas no hay OrgDb, o lo hay pero cubre una fraccion del genoma.
#   - Colecciones curadas que no son ontologias: MSigDB (Hallmark, C2, C7),
#     regulones de RegulonDB, firmas propias del laboratorio.
# El precio es que el fichero manda: si sus identificadores no son los de la
# matriz de conteos el solapamiento es cero y ningun OrgDb lo arregla. Por eso
# se reporta cuantos conjuntos se han leido y de que tamano son.

#' Lee un .gmt y devuelve el TERM2GENE que consumen enricher() y GSEA().
#'
#' Formato GMT: una linea por conjunto, con nombre, descripcion y genes
#' separados por tabuladores.
#'
#' @return list(term2gene, n_sets, sizes, n_genes, error)
read_gene_sets_gmt <- function(path) {
  fail <- function(msg) list(term2gene = NULL, n_sets = 0L, sizes = integer(0),
                             n_genes = 0L, error = msg)
  if (is.null(path) || !length(path) || !nzchar(path[1]) || !file.exists(path[1])) {
    return(fail("No hay fichero GMT cargado."))
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  t2g <- tryCatch(clusterProfiler::read.gmt(path[1]), error = function(e) e)
  if (inherits(t2g, "error")) {
    return(fail(paste0("No se pudo leer el GMT: ", conditionMessage(t2g))))
  }
  if (!is.data.frame(t2g) || !all(c("term", "gene") %in% names(t2g)) || !nrow(t2g)) {
    return(fail("El GMT no tiene el formato esperado (nombre, descripcion, genes)."))
  }
  # read.gmt puede devolver factores segun la version; enricher() los admite,
  # pero el resto del codigo compara con character y un factor rompe intersect().
  t2g$term <- as.character(t2g$term)
  t2g$gene <- as.character(t2g$gene)
  t2g <- t2g[!is.na(t2g$gene) & nzchar(t2g$gene), c("term", "gene"), drop = FALSE]
  t2g <- unique(t2g)
  if (!nrow(t2g)) return(fail("El GMT no contiene genes."))
  sizes <- as.integer(table(t2g$term))
  list(term2gene = t2g, n_sets = length(unique(t2g$term)), sizes = sizes,
       n_genes = length(unique(t2g$gene)), error = NULL)
}

#' Texto de una linea con lo leido del GMT, para mostrarlo junto al selector.
gmt_summary_text <- function(gs) {
  if (is.null(gs)) return(NULL)
  if (!is.null(gs$error)) return(gs$error)
  paste0(fmt_int(gs$n_sets), " conjuntos leidos; tamanos entre ", min(gs$sizes),
         " y ", max(gs$sizes), " genes (mediana ",
         round(stats::median(gs$sizes)), "); ", fmt_int(gs$n_genes),
         " genes distintos.")
}

#' Tasa de mapeo contra los genes de un TERM2GENE.
#'
#' El equivalente de gene_mapping_rate() cuando no hay OrgDb: el "diccionario"
#' es el propio GMT.
mapping_against_sets <- function(genes, term2gene, etiqueta = "GMT") {
  genes <- unique(as.character(genes[!is.na(genes)]))
  out <- list(n_input = length(genes), n_mapped = NA_integer_, rate = NA_real_,
              keytype = etiqueta, source = "GMT")
  if (!length(genes) || is.null(term2gene) || !nrow(term2gene)) return(out)
  out$n_mapped <- length(intersect(genes, unique(term2gene$gene)))
  out$rate <- out$n_mapped / length(genes)
  out
}

#' ORA sobre gene sets propios, via clusterProfiler::enricher().
#'
#' Mismo test hipergeometrico que enrichGO/enrichKEGG y mismo `universe`: el
#' fondo sigue siendo la lista de genes testeados, no todos los genes del GMT.
run_enrichment_gmt <- function(genes, universe = NULL, term2gene = NULL,
                               pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                               minGSSize = 10, maxGSSize = 500) {
  obj_s4 <- NULL
  fail <- function(msg, mapping = NULL) list(table = NULL, error = msg, mapping = mapping)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  if (is.null(term2gene) || !is.data.frame(term2gene) || !nrow(term2gene)) {
    return(fail("Carga un fichero GMT con los conjuntos de genes."))
  }
  if (!length(genes)) return(fail("Lista de genes vacia."))
  mapping <- mapping_against_sets(genes, term2gene)

  out <- tryCatch({
    eg <- clusterProfiler::enricher(
      gene          = unique(as.character(genes)),
      universe      = if (length(universe)) unique(as.character(universe)) else NULL,
      TERM2GENE     = term2gene,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize
    )
    obj_s4 <- eg   # el objeto S4 lo necesitan los graficos de enrichplot
    df <- if (is.null(eg)) NULL else as.data.frame(eg)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e)

  if (inherits(out, "error")) return(fail(conditionMessage(out), mapping))
  if (is.null(out)) return(fail("Sin conjuntos GMT enriquecidos.", mapping))
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], obj = obj_s4, error = NULL, mapping = mapping)
}

# ── ORA direccional ─────────────────────────────────────────────────────────

#' Ejecuta el mismo ORA sobre el conjunto, los genes al alza y los genes a la baja.
#'
#' Mezclar las dos direcciones diluye la senal: una ruta con la mitad de sus
#' genes reprimidos y la otra mitad inducidos aporta el mismo numero de genes al
#' test hipergeometrico que una ruta coherentemente inducida, y sale con el
#' mismo p-valor pese a significar cosas distintas. Separar por signo devuelve
#' esa direccionalidad, que es justo lo que el ORA pierde y GSEA conserva en el
#' NES.
#'
#' `runner` es una funcion de un solo argumento (el vector de genes) que
#' devuelve el mismo list(table, error, mapping) que run_enrichment_*. Asi esta
#' funcion no necesita saber si detras hay GO, KEGG o un GMT.
#'
#' @return list(table, error, mapping, errores, n_genes)
run_ora_directional <- function(deg_df, runner, include_all = TRUE) {
  fail <- function(msg) list(table = NULL, error = msg, mapping = NULL,
                             errores = character(0), n_genes = list())
  if (is.null(deg_df) || !nrow(deg_df)) return(fail("Lista de genes vacia."))
  if (!all(c("gene", "log2FC") %in% names(deg_df))) {
    return(fail("La tabla necesita las columnas 'gene' y 'log2FC' para separar por direccion."))
  }
  lfc <- deg_df$log2FC
  grupos <- list()
  if (isTRUE(include_all)) grupos[["Conjunto"]] <- deg_df$gene
  grupos[["Al alza"]]   <- deg_df$gene[!is.na(lfc) & lfc > 0]
  grupos[["A la baja"]] <- deg_df$gene[!is.na(lfc) & lfc < 0]

  tablas <- list(); errores <- character(0); n_genes <- list(); mapping <- NULL
  for (nm in names(grupos)) {
    g <- unique(as.character(grupos[[nm]]))
    n_genes[[nm]] <- length(g)
    if (!length(g)) {
      errores[nm] <- "Sin genes en esta direccion."
      next
    }
    r <- runner(g)
    # El mapeo del conjunto describe la lista entera; los de las mitades son el
    # mismo diccionario sobre menos genes y no anaden informacion.
    if (is.null(mapping) || identical(nm, "Conjunto")) mapping <- r$mapping %||% mapping
    if (is.null(r$table) || !nrow(r$table)) {
      errores[nm] <- r$error %||% "Sin terminos enriquecidos."
      next
    }
    tb <- r$table
    tb$Direccion <- nm
    tablas[[nm]] <- tb[, c("Direccion", setdiff(names(tb), "Direccion")), drop = FALSE]
  }

  if (!length(tablas)) {
    msg <- if (length(errores)) paste0(names(errores), ": ", errores, collapse = " | ")
           else "Sin terminos enriquecidos."
    return(list(table = NULL, error = msg, mapping = mapping,
                errores = errores, n_genes = n_genes))
  }
  tab <- do.call(rbind, unname(tablas))
  rownames(tab) <- NULL
  list(table = tab, error = NULL, mapping = mapping,
       errores = errores, n_genes = n_genes)
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
#' @param solo_testables si TRUE (por defecto) el ranking se limita a los genes
#'   EVALUABLES, es decir los que tienen p-valor ajustado.
#'
#'   Que GSEA use el ranking completo significa que no se umbraliza por
#'   SIGNIFICACION, y eso es correcto: esa es su ventaja sobre el ORA. Pero el
#'   filtrado independiente es otra cosa: los genes sin `padj` no eran
#'   evaluables, y meterlos en el ranking rompe la coherencia con el universo
#'   del ORA y sesga el resultado.
#'
#'   El sesgo es medible. Con conjuntos de 100 genes tomados AL AZAR del
#'   universo evaluable, 8 de cada 10 salian "significativos" solo por estar
#'   definidos sobre un subconjunto distinto del ranking; tomados del ranking
#'   completo, 0 de 10. Los genes no evaluables (baja expresion) se concentran
#'   en la cola, asi que cualquier conjunto realista queda desplazado hacia la
#'   cabeza por construccion. Ademas los empates pasaban del 18 % al 0,05 %:
#'   los p-valores de los genes sin evaluar saturan y GSEA no resuelve empates.
deg_ranking_metric <- function(deg_df, metric = c("stat", "log2FC", "signed_p"),
                               solo_testables = TRUE) {
  metric <- match.arg(metric)
  if (is.null(deg_df) || !nrow(deg_df)) {
    return(list(ranked = NULL, error = "Sin tabla DEG.", metric = metric))
  }
  n_total <- nrow(deg_df)
  if (isTRUE(solo_testables) && "padj" %in% names(deg_df)) {
    testables <- !is.na(deg_df$padj)
    # Si ningun gen tiene padj (motores que no lo rellenan), se usa la tabla
    # entera: es preferible un ranking imperfecto a no poder correr GSEA.
    if (any(testables)) deg_df <- deg_df[testables, , drop = FALSE]
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
  list(ranked = v, metric = metric, n = length(v), n_total = n_total,
       n_descartados = n_total - nrow(deg_df), ties_frac = ties_frac, error = NULL)
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
#'
#' @param minGSSize,maxGSSize Tamano minimo y maximo de conjunto que se testea.
#'   Los conjuntos muy pequenos dan NES inestables y los muy grandes son
#'   demasiado inespecificos para interpretarse; ademas recortan el numero de
#'   tests y con ello el castigo del ajuste multiple.
#' @param pvalueCutoff Filtro sobre el p-valor AJUSTADO con el que
#'   clusterProfiler recorta la tabla devuelta. Con 1 no filtra nada, que es la
#'   unica forma de distinguir "he testeado y nada llega a 0,05" de un fallo:
#'   con el corte puesto, los dos casos devuelven una tabla vacia.
#' @param eps Cota inferior de los p-valores del multilevel de fgsea. Por
#'   defecto trunca en 1e-10, asi que todos los conjuntos muy significativos
#'   empatan en ese valor y su orden relativo se pierde; con `eps = 0` fgsea
#'   estima el p-valor exacto a costa de mas tiempo de calculo.
#' @param term2gene Data frame (term, gene) para `ont = "GMT"`: enriquecimiento
#'   contra conjuntos propios, sin pasar por el OrgDb.
#' @param seed semilla del remuestreo de fgsea.
#'
#'   `seed = TRUE` de clusterProfiler NO basta: DOSE lo implementa como
#'   `set.seed(.Random.seed)`, y `.Random.seed[1]` es el CODIGO DEL TIPO de
#'   generador (10403 por defecto), no una semilla. El resultado es determinista
#'   solo por accidente —mientras nadie cambie `RNGkind()`— y en una sesion sin
#'   aleatoriedad previa `.Random.seed` ni siquiera existe, lo que convierte un
#'   GSEA perfectamente valido en un "fallo" incomprensible. Se fija aqui de
#'   forma explicita con `withr::with_seed`, que ademas restaura el estado al
#'   salir.
run_gsea <- function(ranked, ont = "BP", OrgDb = NULL, organism = "eco",
                     keyType = "SYMBOL", exponent = 0,
                     pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500,
                     eps = 1e-10, term2gene = NULL, readable = FALSE,
                     seed = 1L) {
  obj_s4 <- NULL
  fail <- function(msg) list(table = NULL, error = msg, mapping = NULL)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(fail("clusterProfiler no esta instalado."))
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    return(fail("fgsea no esta instalado."))
  }
  if (is.null(ranked) || !length(ranked)) return(fail("Ranking de genes vacio."))

  is_kegg <- identical(ont, "KEGG")
  is_gmt  <- identical(ont, "GMT")
  is_reactome <- identical(ont, "REACTOME")
  if (is_gmt) {
    if (is.null(term2gene) || !is.data.frame(term2gene) || !nrow(term2gene)) {
      return(fail("Carga un fichero GMT con los conjuntos de genes."))
    }
    mapping <- mapping_against_sets(names(ranked), term2gene)
    out <- withr::with_seed(seed, tryCatch({
      gs <- clusterProfiler::GSEA(
        geneList = ranked, TERM2GENE = term2gene, exponent = exponent,
        minGSSize = minGSSize, maxGSSize = maxGSSize, eps = eps,
        pvalueCutoff = pvalueCutoff, verbose = FALSE, seed = TRUE  # ver nota de `seed`
      )
      obj_s4 <- gs   # el objeto S4 lo necesitan los graficos de enrichplot
      df <- if (is.null(gs)) NULL else as.data.frame(gs)
      if (is.null(df) || !nrow(df)) NULL else df
    }, error = function(e) e))
    if (inherits(out, "error")) {
      return(list(table = NULL, error = conditionMessage(out), mapping = mapping))
    }
    if (is.null(out)) {
      return(list(table = NULL, error = "Sin conjuntos de genes enriquecidos.",
                  mapping = mapping))
    }
    keep <- intersect(c("ID", "Description", "setSize", "enrichmentScore", "NES",
                        "pvalue", "p.adjust", "qvalue", "core_enrichment"), names(out))
    return(list(table = out[, keep, drop = FALSE], obj = obj_s4, error = NULL, mapping = mapping))
  }
  if (is_reactome) {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      return(fail("ReactomePA no esta instalado."))
    }
    if (!organism %in% REACTOME_ORGANISMOS) {
      return(fail(paste0("Reactome no cubre el organismo '", organism, "'.")))
    }
    # El ranking viaja traducido a ENTREZID, y es el traducido el que hay que
    # devolver: el running score se dibuja sobre el mismo espacio de IDs en el
    # que se calculo el NES, no sobre el original.
    tr <- translate_ranking_to_entrez(ranked, OrgDb, keyType)
    if (is.null(tr$ranked)) {
      return(list(table = NULL, error = tr$error, mapping = tr$mapping))
    }
    args <- list(geneList = tr$ranked, organism = organism, exponent = exponent,
                 minGSSize = minGSSize, maxGSSize = maxGSSize,
                 pvalueCutoff = pvalueCutoff, verbose = FALSE,
                 seed = TRUE)  # ver nota de `seed`
    out <- withr::with_seed(seed, tryCatch({
      # `eps` llego a gsePathway despues que a gseGO: si esta version no lo
      # acepta se repite la llamada sin el, en lugar de perder la coleccion
      # entera por un argumento opcional.
      gs <- tryCatch(do.call(ReactomePA::gsePathway, c(args, list(eps = eps))),
                     error = function(e) do.call(ReactomePA::gsePathway, args))
      if (isTRUE(readable)) gs <- enrich_set_readable(gs, as_orgdb_object(OrgDb), "ENTREZID")
      obj_s4 <- gs   # el objeto S4 lo necesitan los graficos de enrichplot
      df <- if (is.null(gs)) NULL else as.data.frame(gs)
      if (is.null(df) || !nrow(df)) NULL else df
    }, error = function(e) e))
    if (inherits(out, "error")) {
      return(list(table = NULL, error = conditionMessage(out), mapping = tr$mapping,
                  ranked_used = tr$ranked))
    }
    if (is.null(out)) {
      return(list(table = NULL, error = "Sin rutas de Reactome enriquecidas.",
                  mapping = tr$mapping, ranked_used = tr$ranked))
    }
    keep <- intersect(c("ID", "Description", "setSize", "enrichmentScore", "NES",
                        "pvalue", "p.adjust", "qvalue", "core_enrichment"), names(out))
    return(list(table = out[, keep, drop = FALSE], obj = obj_s4, error = NULL,
                mapping = tr$mapping, ranked_used = tr$ranked))
  }
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

  out <- withr::with_seed(seed, tryCatch({
    gs <- if (is_kegg) {
      clusterProfiler::gseKEGG(
        geneList = ranked, organism = organism, keyType = keyType,
        exponent = exponent, minGSSize = minGSSize, maxGSSize = maxGSSize,
        eps = eps, pvalueCutoff = pvalueCutoff, verbose = FALSE, seed = TRUE  # ver nota de `seed`
      )
    } else {
      clusterProfiler::gseGO(
        geneList = ranked, ont = ont, OrgDb = OrgDb, keyType = keyType,
        exponent = exponent, minGSSize = minGSSize, maxGSSize = maxGSSize,
        eps = eps, pvalueCutoff = pvalueCutoff, verbose = FALSE, seed = TRUE  # ver nota de `seed`
      )
    }
    # Aplica a core_enrichment, que es la columna que se lee gen a gen.
    if (isTRUE(readable) && !is_kegg) gs <- enrich_set_readable(gs, OrgDb, keyType)
    obj_s4 <- gs   # el objeto S4 lo necesitan los graficos de enrichplot
    df <- if (is.null(gs)) NULL else as.data.frame(gs)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e))

  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out),
                                          mapping = mapping))
  if (is.null(out)) return(list(table = NULL, error = "Sin conjuntos de genes enriquecidos.",
                                mapping = mapping))
  keep <- intersect(c("ID", "Description", "setSize", "enrichmentScore", "NES",
                      "pvalue", "p.adjust", "qvalue", "core_enrichment"), names(out))
  list(table = out[, keep, drop = FALSE], obj = obj_s4, error = NULL, mapping = mapping)
}

# ── Running score y leading edge ────────────────────────────────────────────

#' Resuelve un OrgDb dado por nombre al objeto, o NULL si no se puede.
as_orgdb_object <- function(OrgDb) {
  if (is.null(OrgDb)) return(NULL)
  if (is.character(OrgDb)) {
    if (!nzchar(OrgDb) || !requireNamespace(OrgDb, quietly = TRUE)) return(NULL)
    return(getFromNamespace(OrgDb, OrgDb))
  }
  OrgDb
}

#' Genes que componen un conjunto, en el mismo espacio de IDs que el ranking.
#'
#' Hace falta para dibujar el running score: el core_enrichment de la tabla solo
#' trae el leading edge (los genes hasta el pico), no el conjunto completo, y sin
#' el conjunto completo no hay curva que dibujar.
#'
#' Para GO se consulta GOALL y no GO: gseGO testea cada termino con todos los
#' genes anotados en el y en sus descendientes (GO2ALLEGS), asi que usar GO
#' devolveria un conjunto mas pequeno que el testeado y una curva que no
#' corresponde al NES de la tabla.
#'
#' @return list(genes, error)
gsea_term_genes <- function(term_id, ont = "BP", OrgDb = NULL, keyType = "SYMBOL",
                            organism = "eco", term2gene = NULL) {
  out <- function(g, err = NULL) list(genes = unique(as.character(g)), error = err)
  if (is.null(term_id) || !length(term_id) || !nzchar(term_id[1])) {
    return(out(character(0), "Sin termino seleccionado."))
  }
  term_id <- term_id[1]

  if (identical(ont, "GMT")) {
    if (is.null(term2gene) || !nrow(term2gene)) return(out(character(0), "Sin GMT cargado."))
    return(out(term2gene$gene[term2gene$term == term_id]))
  }

  if (identical(ont, "REACTOME")) {
    # Los genes de la ruta salen de reactome.db, que es local: a diferencia de
    # KEGG no hay consulta en linea y la curva se puede dibujar sin red. Vienen
    # en ENTREZID, que es el espacio en el que se corrio el GSEA.
    if (!requireNamespace("reactome.db", quietly = TRUE)) {
      return(out(character(0), "reactome.db no esta instalado."))
    }
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      return(out(character(0), "AnnotationDbi no esta instalado."))
    }
    # `ifnotfound` solo admite NA en el mget de AnnotationDbi: pasarle list(NULL)
    # —que es la forma habitual en el mget de base— lanza un error, y envuelto en
    # un tryCatch se convierte en "la ruta no tiene genes". El sintoma no es una
    # excepcion sino una curva vacia sobre un conjunto que si existe.
    g <- tryCatch({
      v <- AnnotationDbi::mget(term_id, reactome.db::reactomePATHID2EXTID,
                               ifnotfound = NA)[[1]]
      if (length(v) == 1L && is.na(v[1])) NULL else as.character(v)
    }, error = function(e) NULL)
    if (is.null(g) || !length(g)) {
      return(out(character(0),
                 paste0("La ruta '", term_id, "' no esta en reactome.db.")))
    }
    return(out(g))
  }

  if (identical(ont, "KEGG")) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      return(out(character(0), "clusterProfiler no esta instalado."))
    }
    # download_KEGG consulta la API en linea; sin red no hay curva, pero eso no
    # debe tumbar la pestana entera.
    kg <- tryCatch(clusterProfiler::download_KEGG(organism), error = function(e) e)
    if (inherits(kg, "error")) {
      return(out(character(0), paste0("No se pudo consultar KEGG: ", conditionMessage(kg))))
    }
    df <- kg$KEGGPATHID2EXTID
    if (is.null(df) || !nrow(df)) return(out(character(0), "KEGG no devolvio conjuntos."))
    return(out(df[[2]][df[[1]] == term_id]))
  }

  db <- as_orgdb_object(OrgDb)
  if (is.null(db) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(out(character(0), "Sin OrgDb para recuperar los genes del termino."))
  }
  res <- tryCatch(
    suppressMessages(AnnotationDbi::select(db, keys = term_id, keytype = "GOALL",
                                           columns = keyType)),
    error = function(e) e
  )
  if (inherits(res, "error")) return(out(character(0), conditionMessage(res)))
  if (is.null(res) || !nrow(res) || !keyType %in% names(res)) {
    return(out(character(0), "El termino no tiene genes anotados en este OrgDb."))
  }
  out(res[[keyType]][!is.na(res[[keyType]])])
}

#' Datos del running score de un conjunto sobre el ranking, via fgsea.
#'
#' `plotEnrichmentData()` devuelve la curva ya calculada (fgsea la usa para su
#' propio grafico), de modo que aqui no se reimplementa el estadistico: se
#' dibuja con plotly, coherente con el resto de la app y sin depender de
#' enrichplot.
#'
#' `gseaParam` debe ser el mismo `exponent` con el que se corrio el GSEA; con 0
#' (permutacion no ponderada, el default de la app) la curva es la de
#' Kolmogorov-Smirnov clasica y no la ponderada por magnitud.
#'
#' @return list(curve, ticks, es, es_pos, es_neg, n_set, n_hits, n_ranked, error)
gsea_running_score <- function(ranked, genes, gseaParam = 0) {
  fail <- function(msg) list(curve = NULL, ticks = NULL, es = NA_real_,
                             es_pos = NA_real_, es_neg = NA_real_,
                             n_set = length(unique(genes %||% character(0))),
                             n_hits = 0L, n_ranked = length(ranked %||% numeric(0)),
                             error = msg)
  if (!requireNamespace("fgsea", quietly = TRUE)) return(fail("fgsea no esta instalado."))
  if (is.null(ranked) || !length(ranked) || is.null(names(ranked))) {
    return(fail("Ranking de genes vacio."))
  }
  genes <- unique(as.character(genes %||% character(0)))
  hits <- intersect(genes, names(ranked))
  if (!length(hits)) {
    return(fail("Ningun gen del conjunto esta en el ranking (revisa el keyType)."))
  }
  d <- tryCatch(fgsea::plotEnrichmentData(pathway = hits, stats = ranked,
                                          gseaParam = gseaParam),
                error = function(e) e)
  if (inherits(d, "error")) return(fail(conditionMessage(d)))
  curve <- as.data.frame(d$curve)
  ticks <- as.data.frame(d$ticks)
  # El ES es el extremo de mayor valor absoluto, con su signo.
  es <- if (abs(d$posES) >= abs(d$negES)) d$posES else d$negES
  list(curve = curve, ticks = ticks, es = es, es_pos = d$posES, es_neg = d$negES,
       n_set = length(genes), n_hits = length(hits), n_ranked = length(ranked),
       error = NULL)
}

#' Separa la cadena core_enrichment ("gen1/gen2/...") en un vector de genes.
leading_edge_genes <- function(core) {
  if (is.null(core) || !length(core) || is.na(core[1]) || !nzchar(core[1])) return(character(0))
  g <- unlist(strsplit(as.character(core[1]), "/", fixed = TRUE))
  g[nzchar(g)]
}

# ── ORA frente a GSEA ───────────────────────────────────────────────────────

#' Solapamiento entre los terminos significativos de un ORA y de un GSEA.
#'
#' Los dos responden a preguntas distintas sobre los mismos datos: el ORA
#' pregunta si la lista umbralizada esta enriquecida en un termino, GSEA si el
#' termino esta desplazado en el ranking completo. Un solapamiento bajo no es un
#' error de ninguno de los dos, y los terminos que solo ve GSEA son justamente
#' los de senal debil pero coordinada que el corte de la lista deja fuera. Es el
#' argumento por el que la app ofrece ambos y no solo el ORA.
#'
#' @return list(n_ora, n_gsea, n_comun, jaccard, solo_gsea, solo_ora, comunes)
compare_ora_gsea <- function(ora_df, gsea_df, padj_cutoff = 0.05) {
  sig_ids <- function(df) {
    if (is.null(df) || !nrow(df) || !"ID" %in% names(df)) return(character(0))
    keep <- if ("p.adjust" %in% names(df)) !is.na(df$p.adjust) & df$p.adjust <= padj_cutoff
            else rep(TRUE, nrow(df))
    unique(as.character(df$ID[keep]))
  }
  ids_o <- sig_ids(ora_df)
  ids_g <- sig_ids(gsea_df)
  comun <- intersect(ids_o, ids_g)
  union_ids <- union(ids_o, ids_g)
  sub <- function(df, ids) {
    if (is.null(df) || !nrow(df) || !length(ids)) return(NULL)
    d <- df[!is.na(df$ID) & df$ID %in% ids, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    d[order(d$p.adjust, na.last = TRUE), , drop = FALSE]
  }
  comunes <- if (length(comun)) {
    o <- ora_df[match(comun, ora_df$ID), , drop = FALSE]
    g <- gsea_df[match(comun, gsea_df$ID), , drop = FALSE]
    data.frame(
      ID          = comun,
      Description = if ("Description" %in% names(g)) as.character(g$Description)
                    else as.character(o$Description),
      padj_ORA    = o$p.adjust,
      padj_GSEA   = g$p.adjust,
      NES         = if ("NES" %in% names(g)) g$NES else NA_real_,
      stringsAsFactors = FALSE
    )
  } else NULL
  list(
    n_ora     = length(ids_o),
    n_gsea    = length(ids_g),
    n_comun   = length(comun),
    # Jaccard sobre los conjuntos de terminos significativos. Sin ninguno de los
    # dos, no es 0 (que sugeriria discrepancia total) sino indefinido.
    jaccard   = if (length(union_ids)) length(comun) / length(union_ids) else NA_real_,
    solo_gsea = sub(gsea_df, setdiff(ids_g, ids_o)),
    solo_ora  = sub(ora_df, setdiff(ids_o, ids_g)),
    comunes   = if (is.null(comunes)) NULL else
                  comunes[order(comunes$padj_GSEA, na.last = TRUE), , drop = FALSE]
  )
}

#' Tabla unica con los terminos de la comparacion y quien los ve.
#'
#' Se etiqueta cada termino en lugar de mostrar dos tablas separadas: lo
#' interesante es leer seguidos los que solo ve GSEA junto a los comunes, y asi
#' se descarga de una vez.
compare_ora_gsea_table <- function(cmp) {
  if (is.null(cmp)) return(NULL)
  fila <- function(df, etiqueta, col_padj) {
    if (is.null(df) || !nrow(df)) return(NULL)
    data.frame(
      Visto       = etiqueta,
      ID          = as.character(df$ID),
      Description = as.character(df$Description %||% NA_character_),
      padj_ORA    = if (identical(col_padj, "ora")) df$p.adjust else NA_real_,
      padj_GSEA   = if (identical(col_padj, "gsea")) df$p.adjust else NA_real_,
      NES         = if ("NES" %in% names(df)) df$NES else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  comunes <- if (is.null(cmp$comunes) || !nrow(cmp$comunes)) NULL else
    data.frame(Visto = "Ambos", ID = cmp$comunes$ID,
               Description = cmp$comunes$Description,
               padj_ORA = cmp$comunes$padj_ORA, padj_GSEA = cmp$comunes$padj_GSEA,
               NES = cmp$comunes$NES, stringsAsFactors = FALSE)
  out <- do.call(rbind, Filter(Negate(is.null), list(
    fila(cmp$solo_gsea, "Solo GSEA", "gsea"),
    comunes,
    fila(cmp$solo_ora, "Solo ORA", "ora")
  )))
  if (is.null(out) || !nrow(out)) return(NULL)
  rownames(out) <- NULL
  out
}

#' Ordena el enriquecimiento por p.adjust y devuelve top_n filas.
#'
#' Con ORA direccional la tabla trae varias direcciones apiladas, y entonces el
#' top se toma DENTRO de cada direccion: si no, la direccion con la senal mas
#' fuerte se lleva las 15 filas y la otra desaparece del grafico, que es
#' precisamente la comparacion que se queria ver.
enrichment_dotplot_data <- function(enrich_df, top_n = 15) {
  if (is.null(enrich_df) || !nrow(enrich_df)) return(NULL)
  ord <- order(enrich_df$p.adjust, na.last = TRUE)
  out <- enrich_df[ord, , drop = FALSE]
  if ("Direccion" %in% names(out)) {
    out <- do.call(rbind, unname(lapply(split(out, out$Direccion), function(d) {
      if (nrow(d) > top_n) d[seq_len(top_n), , drop = FALSE] else d
    })))
    out <- out[order(out$p.adjust, na.last = TRUE), , drop = FALSE]
    # Etiqueta unica: el mismo termino puede salir en dos direcciones y en el eje
    # se colapsarian en una sola fila.
    out$plot_label <- paste0(out$Description, " [", out$Direccion, "]")
  } else {
    if (nrow(out) > top_n) out <- out[seq_len(top_n), , drop = FALSE]
    out$plot_label <- as.character(out$Description)
  }
  rownames(out) <- NULL
  # Convertir GeneRatio "x/y" a numerico
  if ("GeneRatio" %in% names(out)) {
    out$GeneRatioNum <- sapply(strsplit(as.character(out$GeneRatio), "/", fixed = TRUE),
                               function(p) if (length(p) == 2L) suppressWarnings(as.numeric(p[1]) / as.numeric(p[2])) else NA_real_)
  }
  out
}

# ── Reactome ────────────────────────────────────────────────────────────────
#
# Por que existe: GO describe funciones y KEGG mapas metabolicos, pero las rutas
# de senalizacion, el ciclo celular o la respuesta inmune estan mucho mejor
# representadas en Reactome, que es una base curada manualmente y revisada por
# pares (Milacic et al., NAR 2024). Es la tercera coleccion de referencia y la
# que cierra el hueco entre las otras dos.
#
# Dos limites que la interfaz tiene que decir, no esconder:
#
#   1. Reactome NO cubre procariotas. Sus organismos son un catalogo cerrado de
#      eucariotas modelo. Con datos de E. coli la respuesta correcta es "esta
#      coleccion no aplica", no una tabla vacia que parece un fallo.
#   2. Reactome trabaja en ENTREZID. Los identificadores de la matriz casi nunca
#      lo son (simbolos, ENSEMBL, locus tags), asi que hay una TRADUCCION por
#      medio, y toda traduccion pierde genes. Esa perdida se mide y se muestra
#      con el mismo criterio que el resto del modulo: la tasa de mapeo viaja
#      siempre en el resultado.

#' Etiqueta legible de la coleccion de enriquecimiento, para informe y auditoria.
#'
#' El codigo interno ("BP", "REACTOME", "GMT") no dice lo mismo a quien lee el
#' informe seis meses despues que a quien escribio la interfaz.
enrich_collection_label <- function(ont) {
  etiquetas <- c(BP = "GO: procesos biologicos", MF = "GO: funcion molecular",
                 CC = "GO: componente celular", KEGG = "KEGG",
                 REACTOME = "Reactome", GMT = "Gene sets propios (GMT)")
  v <- unname(etiquetas[ont %||% ""])
  if (is.na(v)) ont %||% "—" else v
}

#' Organismos que cubre Reactome, con el nombre que espera ReactomePA.
REACTOME_ORGANISMOS <- c(
  "Homo sapiens"             = "human",
  "Mus musculus"             = "mouse",
  "Rattus norvegicus"        = "rat",
  "Danio rerio"              = "zebrafish",
  "Drosophila melanogaster"  = "fly",
  "Caenorhabditis elegans"   = "celegans",
  "Saccharomyces cerevisiae" = "yeast",
  "Sus scrofa"               = "pig",
  "Bos taurus"               = "cow",
  "Canis familiaris"         = "dog",
  "Gallus gallus"            = "chicken",
  "Xenopus tropicalis"       = "xenopus"
)

#' Organismo de Reactome que corresponde a un OrgDb, o NULL si no hay ninguno.
#'
#' Sirve para preseleccionar el organismo a partir del OrgDb ya elegido, en vez
#' de obligar a declararlo dos veces y arriesgar que no coincidan.
reactome_organism_for_orgdb <- function(pkg) {
  if (is.null(pkg) || !length(pkg) || !nzchar(pkg[1] %||% "")) return(NULL)
  mapa <- c(org.Hs.eg.db = "human", org.Mm.eg.db = "mouse", org.Rn.eg.db = "rat",
            org.Dr.eg.db = "zebrafish", org.Dm.eg.db = "fly",
            org.Ce.eg.db = "celegans", org.Sc.sgd.db = "yeast",
            org.Ss.eg.db = "pig", org.Bt.eg.db = "cow", org.Cf.eg.db = "dog",
            org.Gg.eg.db = "chicken", org.Xl.eg.db = "xenopus")
  v <- unname(mapa[pkg[1]])
  if (is.na(v)) NULL else v
}

#' Traduce identificadores de gen a ENTREZID a traves de un OrgDb.
#'
#' `multiVals = "first"` cuando un identificador de entrada mapea a varios
#' ENTREZID: es el comportamiento de clusterProfiler::bitr y el unico
#' determinista sin criterio biologico adicional.
#'
#' @return list(ids, back, mapping, error) donde `ids` son los ENTREZID unicos y
#'   `back` es un vector con nombres ENTREZID y valores el identificador
#'   original, para poder deshacer la traduccion al mostrar resultados.
translate_to_entrez <- function(genes, OrgDb, keyType = "SYMBOL") {
  genes <- unique(as.character(genes %||% character(0)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  vacio <- function(msg) list(
    ids = character(0), back = character(0), error = msg,
    mapping = list(n_input = length(genes), n_mapped = 0L, rate = 0,
                   keytype = keyType, source = "ENTREZID"))

  if (!length(genes)) return(vacio("Lista de genes vacia."))
  # Ya estan en el espacio de destino: no hay traduccion que hacer ni que medir.
  if (identical(keyType, "ENTREZID")) {
    return(list(ids = genes, back = stats::setNames(genes, genes), error = NULL,
                mapping = list(n_input = length(genes), n_mapped = length(genes),
                               rate = 1, keytype = "ENTREZID", source = "ENTREZID")))
  }
  db <- as_orgdb_object(OrgDb)
  if (is.null(db)) return(vacio("Sin OrgDb con el que traducir a ENTREZID."))
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(vacio("AnnotationDbi no esta instalado."))
  }
  m <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(db, keys = genes, column = "ENTREZID",
                                           keytype = keyType, multiVals = "first")),
    error = function(e) e)
  if (inherits(m, "error")) return(vacio(conditionMessage(m)))

  ok <- !is.na(m) & nzchar(as.character(m))
  entrez <- as.character(m[ok])
  originales <- genes[ok]
  # Varios identificadores de entrada pueden compartir ENTREZID (sinonimos): se
  # queda el primero, y el conteo de mapeados es el de ENTREZID distintos.
  dup <- duplicated(entrez)
  entrez <- entrez[!dup]; originales <- originales[!dup]
  list(ids = entrez, back = stats::setNames(originales, entrez), error = NULL,
       mapping = list(n_input = length(genes), n_mapped = length(entrez),
                      rate = if (length(genes)) length(entrez) / length(genes) else NA_real_,
                      keytype = keyType, source = "ENTREZID"))
}

#' Traduce un ranking con nombres a ENTREZID conservando el orden por valor.
#'
#' Cuando dos identificadores caen en el mismo ENTREZID se conserva el de MAYOR
#' valor absoluto. Quedarse con el primero haria depender el resultado del orden
#' de las filas de la matriz, que es arbitrario; quedarse con el mas extremo es
#' reproducible y conserva la senal que GSEA va a leer.
translate_ranking_to_entrez <- function(ranked, OrgDb, keyType = "SYMBOL") {
  if (is.null(ranked) || !length(ranked) || is.null(names(ranked))) {
    return(list(ranked = NULL, error = "Ranking de genes vacio.", mapping = NULL))
  }
  tr <- translate_to_entrez(names(ranked), OrgDb, keyType)
  if (!length(tr$ids)) {
    return(list(ranked = NULL, mapping = tr$mapping,
                error = tr$error %||% "Ningun gen del ranking se pudo traducir a ENTREZID."))
  }
  # La traduccion se rehace sobre el ranking completo (no sobre `tr$ids`, que ya
  # esta deduplicado) para poder elegir por magnitud.
  db <- as_orgdb_object(OrgDb)
  m <- if (identical(keyType, "ENTREZID")) {
    stats::setNames(names(ranked), names(ranked))
  } else {
    tryCatch(suppressMessages(AnnotationDbi::mapIds(
      db, keys = names(ranked), column = "ENTREZID", keytype = keyType,
      multiVals = "first")), error = function(e) NULL)
  }
  if (is.null(m)) {
    return(list(ranked = NULL, mapping = tr$mapping,
                error = "No se pudo traducir el ranking a ENTREZID."))
  }
  ent <- as.character(m)
  ok <- !is.na(ent) & nzchar(ent)
  v <- ranked[ok]; ent <- ent[ok]
  ord <- order(abs(v), decreasing = TRUE)
  v <- v[ord]; ent <- ent[ord]
  keep <- !duplicated(ent)
  v <- v[keep]; names(v) <- ent[keep]
  v <- sort(v, decreasing = TRUE)
  list(ranked = v, mapping = tr$mapping, error = NULL)
}

#' ORA sobre rutas de Reactome, via ReactomePA::enrichPathway().
#'
#' El universo se traduce con la MISMA funcion que la lista: si se tradujera solo
#' la lista, el fondo dejaria de ser el conjunto de genes testeados y volveriamos
#' al error que Wijesooriya et al. (2022) encontraron en el 95 % de los ORA
#' publicados.
run_enrichment_reactome <- function(genes, universe = NULL, OrgDb = NULL,
                                    keyType = "SYMBOL", organism = "human",
                                    pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                                    minGSSize = 10, maxGSSize = 500,
                                    readable = FALSE) {
  obj_s4 <- NULL
  fail <- function(msg, mapping = NULL) list(table = NULL, error = msg, mapping = mapping)
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    return(fail("ReactomePA no esta instalado."))
  }
  if (!organism %in% REACTOME_ORGANISMOS) {
    return(fail(paste0("Reactome no cubre el organismo '", organism,
                       "'. Organismos disponibles: ",
                       paste(sort(unname(REACTOME_ORGANISMOS)), collapse = ", "), ".")))
  }
  if (!length(genes)) return(fail("Lista de genes vacia."))

  tr <- translate_to_entrez(genes, OrgDb, keyType)
  if (!length(tr$ids)) {
    return(fail(tr$error %||% paste0(
      "Ningun gen se pudo traducir a ENTREZID desde keyType = ", keyType,
      ". Reactome trabaja en ENTREZID; revisa el tipo de identificador."),
      tr$mapping))
  }
  uni <- if (length(universe)) translate_to_entrez(universe, OrgDb, keyType)$ids else NULL

  out <- tryCatch({
    ep <- ReactomePA::enrichPathway(
      gene          = tr$ids,
      universe      = if (length(uni)) uni else NULL,
      organism      = organism,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      readable      = FALSE
    )
    # setReadable con keytype ENTREZID: los IDs del resultado ya estan en ese
    # espacio, no en el de entrada.
    if (isTRUE(readable)) ep <- enrich_set_readable(ep, as_orgdb_object(OrgDb), "ENTREZID")
    obj_s4 <- ep   # el objeto S4 lo necesitan los graficos de enrichplot
    df <- if (is.null(ep)) NULL else as.data.frame(ep)
    if (is.null(df) || !nrow(df)) NULL else df
  }, error = function(e) e)

  if (inherits(out, "error")) return(fail(conditionMessage(out), tr$mapping))
  if (is.null(out)) return(fail("Sin rutas de Reactome enriquecidas.", tr$mapping))
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "Count", "geneID"), names(out))
  list(table = out[, keep, drop = FALSE], obj = obj_s4, error = NULL, mapping = tr$mapping)
}
