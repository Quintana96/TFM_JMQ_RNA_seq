#' utils_enrich_plots.R
#' Representaciones de red y distribución del enriquecimiento funcional.
#'
#' Completan el esquema del pipeline de RNA-seq de la memoria: el dotplot y el
#' running score ya estaban en la aplicación, y aquí se añaden la barra, la red
#' gen-concepto, el mapa de términos, el upset y el ridge.
#'
#' A diferencia del resto de gráficos de la aplicación, estos NO son de plotly.
#' Se apoyan en enrichplot, que devuelve objetos de ggplot2, y necesitan el
#' objeto S4 del enriquecimiento (enrichResult / gseaResult), no la tabla: la
#' tabla no conserva ni la pertenencia gen-término ni el ranking. Por eso los
#' runners de utils_enrich.R devuelven `obj` además de `table`.
#'
#' Convertirlos a plotly con ggplotly() no compensa: cnetplot y emapplot son
#' grafos dibujados con capas que plotly no traduce, y el resultado pierde las
#' aristas. Para la memoria interesa más la descarga en PNG a 300 ppp que la
#' interactividad.

# ── Catalogo ────────────────────────────────────────────────────────────────

#' Los tipos disponibles dependen del enfoque: el ORA no tiene NES que ordenar
#' en un ridge, y el GSEA no tiene una lista umbralizada que contar en barras.
ENRICH_PLOT_TIPOS <- list(
  barra = list(etiqueta = "Barras",            enfoques = "ora"),
  cnet  = list(etiqueta = "Red gen-concepto",  enfoques = c("ora", "gsea")),
  emap  = list(etiqueta = "Mapa de términos",  enfoques = c("ora", "gsea")),
  upset = list(etiqueta = "Solapamiento (upset)", enfoques = c("ora", "gsea")),
  ridge = list(etiqueta = "Distribución (ridge)", enfoques = "gsea")
)

#' Enfoque de un objeto de enriquecimiento: "gsea" si trae NES, "ora" si no.
enrich_obj_enfoque <- function(obj) {
  if (is.null(obj)) return(NA_character_)
  if (inherits(obj, "gseaResult")) "gsea" else "ora"
}

#' Tipos aplicables a un objeto, como vector con nombres para selectInput().
enrich_plot_choices <- function(obj) {
  enf <- enrich_obj_enfoque(obj)
  if (is.na(enf)) return(character(0))
  ok <- vapply(ENRICH_PLOT_TIPOS, function(x) enf %in% x$enfoques, logical(1))
  stats::setNames(names(ENRICH_PLOT_TIPOS)[ok],
                  vapply(ENRICH_PLOT_TIPOS[ok], function(x) x$etiqueta, character(1)))
}

# ── Similitud entre términos ────────────────────────────────────────────────

#' emapplot y treeplot necesitan la matriz de similitud, que pairwise_termsim()
#' calcula y guarda en el propio objeto. No es gratis —compara cada término con
#' todos los demas—, así que se limita a los términos que se van a dibujar.
#'
#' Devuelve el objeto anotado, o una condición si no se pudo calcular.
enrich_con_similitud <- function(obj, top_n = 30) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("enrichplot no está instalado.")
  }
  tryCatch(
    enrichplot::pairwise_termsim(obj, showCategory = top_n),
    error = function(e) {
      # La similitud semantica de GO necesita GOSemSim y el OrgDb; con
      # colecciones sin ontologia (GMT, KEGG) enrichplot cae en Jaccard sobre
      # los genes, que no requiere nada más. Si aun así falla, el mensaje real
      # es más útil que un gráfico vacio.
      stop(paste0("No se pudo calcular la similitud entre términos: ",
                  conditionMessage(e)))
    }
  )
}

# ── Construcción del gráfico ────────────────────────────────────────────────

#' Dibuja uno de los tipos de ENRICH_PLOT_TIPOS sobre el objeto S4.
#'
#' @param obj enrichResult o gseaResult devuelto por los runners.
#' @param tipo clave de ENRICH_PLOT_TIPOS.
#' @param top_n número de categorías a mostrar.
#' @param fold_change vector con nombres (gen -> log2FC) para colorear la red
#'   gen-concepto. Sus nombres tienen que estar en el mismo espacio de
#'   identificadores que el enriquecimiento: si se corrio con readable = TRUE
#'   son simbolos, y si no, el keyType elegido. Un vector que no case no rompe
#'   el gráfico, simplemente deja los genes sin color.
#' @return un objeto de ggplot2, o un error con un mensaje legible.
enrich_make_network_plot <- function(obj, tipo, top_n = 15, fold_change = NULL,
                                     etiquetar_genes = FALSE) {
  if (is.null(obj)) stop("No hay resultado de enriquecimiento sobre el que dibujar.")
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("enrichplot no está instalado.")
  }
  if (!tipo %in% names(ENRICH_PLOT_TIPOS)) stop(paste0("Tipo desconocido: ", tipo))
  enf <- enrich_obj_enfoque(obj)
  if (!enf %in% ENRICH_PLOT_TIPOS[[tipo]]$enfoques) {
    stop(paste0("'", ENRICH_PLOT_TIPOS[[tipo]]$etiqueta, "' no aplica a un resultado de ",
                toupper(enf), "."))
  }
  n_terminos <- nrow(as.data.frame(obj))
  if (!n_terminos) stop("El enriquecimiento no devolvio términos.")
  # as.numeric() y no as.integer(): ridgeplot() de enrichplot 1.30 comprueba el
  # tipo de showCategory de una forma que da "integer" por inválido, avisa de que
  # "should be a number of pathways" y muere después con "objeto 'selected' no
  # encontrado". Con un doble entra por la rama correcta. No es cosmetico.
  top_n <- max(2, min(round(as.numeric(top_n)), n_terminos))

  # enrichplot 1.30 pasa `by = "Count"` a fortify(), que ggplot2 3.5 ya no
  # acepta en `...`: sale un aviso en cada dibujo de barras y de ridge. Es ruido
  # de una incompatibilidad entre las dos, no un problema del dato, y aparece en
  # la consola de Shiny en cada render. Se silencia SOLO ese aviso, por su texto:
  # cualquier otro sigue llegando.
  con_avisos_de_enrichplot <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("must be used|Problematic argument", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
  }

  con_avisos_de_enrichplot(switch(
    tipo,
    # barplot es un método S3 sobre graphics::barplot: se registra al cargar el
    # espacio de nombres de enrichplot, de modo que despacha sin necesidad de
    # adjuntar el paquete.
    barra = graphics::barplot(obj, showCategory = top_n),
    cnet  = {
      # En enrichplot 1.30 cnetplot pasa por ggtangle y recupera el argumento
      # `foldChange` directo. El `color.params = list(foldChange = ...)` de
      # versiones intermedias da "el argumento no fue usado".
      #
      # El tope de 5 categorías no es arbitrario: un término GO de bacteria
      # arrastra facilmente 200 genes, y con ocho términos el gráfico deja de
      # ser legible por muchos pixeles que se le den. Por lo mismo las etiquetas
      # de gen van apagadas: con 500 nodos se solapan hasta tapar la estructura,
      # que es lo único que el gráfico tiene que comunicar.
      args <- list(obj, showCategory = min(top_n, 5),
                   node_label = if (isTRUE(etiquetar_genes)) "all" else "category")
      if (length(fold_change)) args$foldChange <- fold_change
      do.call(enrichplot::cnetplot, args)
    },
    emap  = enrichplot::emapplot(enrich_con_similitud(obj, top_n), showCategory = top_n),
    upset = enrichplot::upsetplot(obj, n = min(top_n, 10L)),
    ridge = enrichplot::ridgeplot(obj, showCategory = top_n)
  ))
}

#' Texto de ayuda de cada tipo, para que el gráfico no haya que interpretarlo
#' de memoria. Se muestra sobre el gráfico en la interfaz.
enrich_plot_ayuda <- function(tipo) {
  # enrichplot 1.30 pasa `by = "Count"` a fortify(), que ggplot2 3.5 ya no
  # acepta en `...`: sale un aviso en cada dibujo de barras y de ridge. Es ruido
  # de una incompatibilidad entre las dos, no un problema del dato, y aparece en
  # la consola de Shiny en cada render. Se silencia SOLO ese aviso, por su texto:
  # cualquier otro sigue llegando.
  con_avisos_de_enrichplot <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("must be used|Problematic argument", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
  }

  switch(
    tipo,
    barra = paste("Número de genes de la lista en cada término, coloreado por",
                  "significacion. Es la lectura más directa del ORA, pero no",
                  "dice nada del tamaño del término en el fondo."),
    cnet  = paste("Cada término se une a los genes que lo sostienen. Sirve para",
                  "ver que un mismo gen aparece en varios términos, que es la",
                  "razón por la que las listas de enriquecimiento parecen más",
                  "largas de lo que son. Si hay log2FC, los genes se colorean",
                  "por su cambio. Se limita a los 5 términos más significativos",
                  "aunque el deslizador pida más: por encima de ahí el gráfico",
                  "deja de leerse."),
    emap  = paste("Cada nodo es un término y cada arista la fraccion de genes",
                  "que comparten. Los grupos que se forman son familias de",
                  "términos redundantes: leerlos como hallazgos independientes",
                  "es el error más comun al interpretar un enriquecimiento."),
    upset = paste("Cuantifica lo que la red sugiere: cuantos genes son",
                  "exclusivos de un término y cuantos caen en intersecciones.",
                  "Las barras de interseccion grandes indican redundancia."),
    ridge = paste("Distribución del estadístico de ordenación de los genes de",
                  "cada conjunto. Un conjunto desplazado a la derecha está",
                  "sobreexpresado en el contraste; uno centrado en cero no",
                  "tiene señal coordinada aunque su p-valor sea bajo."),
    ""
  )
}

# ── Diagrama de ruta KEGG (pathview) ────────────────────────────────────────

#' Pinta los log2FC sobre el diagrama oficial de una ruta KEGG.
#'
#' pathview no devuelve el gráfico: descarga el KGML y el PNG de la ruta desde
#' KEGG y ESCRIBE los ficheros en el directorio de trabajo. Por eso aquí se
#' cambia de directorio a uno temporal y se restaura siempre, y por eso la
#' función devuelve una ruta de fichero y no un objeto.
#'
#' Requiere además que pathview este ADJUNTO, no solo cargado: internamente
#' llama a data(bods) sin declarar el paquete, y con requireNamespace() falla
#' con "objeto 'bods' no encontrado".
#'
#' @param pathway_id identificador con o sin prefijo de organismo (eco00020 o 00020).
#' @param fold_change vector con nombres = identificadores del tipo `gene_idtype`.
#' @param species código KEGG del organismo ("eco", "hsa", ...).
#' @return list(path = ruta al PNG, error = NULL) o list(path = NULL, error = msg).
enrich_pathview_png <- function(pathway_id, fold_change, species = "eco",
                                gene_idtype = "ENTREZ", limite = 3,
                                dir_salida = NULL) {
  if (!requireNamespace("pathview", quietly = TRUE)) {
    return(list(path = NULL, error = "pathview no está instalado."))
  }
  if (is.null(fold_change) || !length(fold_change)) {
    return(list(path = NULL, error = "No hay log2FC que pintar sobre la ruta."))
  }
  pid <- sub(paste0("^", species), "", trimws(pathway_id %||% ""))
  if (!nzchar(pid)) return(list(path = NULL, error = "Selecciona una ruta KEGG."))

  wd <- dir_salida %||% file.path(tempdir(), paste0("pathview_", pid))
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  antiguo <- setwd(wd)
  on.exit(setwd(antiguo), add = TRUE)

  sufijo <- "tfm"
  out <- tryCatch({
    # attachNamespace en lugar de library() para no ensuciar la sesión de Shiny
    # más de lo necesario, pero con el mismo efecto sobre data(bods).
    if (!"package:pathview" %in% search()) {
      suppressMessages(attachNamespace("pathview"))
    }
    suppressMessages(pathview::pathview(
      gene.data = fold_change, pathway.id = pid, species = species,
      gene.idtype = gene_idtype, out.suffix = sufijo, kegg.dir = wd,
      limit = list(gene = limite), low = list(gene = "#2C7BB6"),
      mid = list(gene = "#F7F7F7"), high = list(gene = "#D7191C")
    ))
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(out)) return(list(path = NULL, error = out))

  png <- list.files(wd, pattern = paste0(sufijo, "\\.png$"), full.names = TRUE)
  if (!length(png)) {
    return(list(path = NULL, error = paste0(
      "KEGG devolvio el diagrama pero ningun gen de la lista pudo situarse en el. ",
      "Comprueba que los identificadores son del tipo '", gene_idtype, "'.")))
  }
  list(path = png[1], error = NULL)
}
