#' scripts/harness/datasets.R
#' Configuración declarativa de los datasets de validación.
#'
#' Un dataset se describe una vez aquí y el harness sabe ejecutarlo entero sin
#' abrir la interfaz. La idea es que añadir un dataset sea añadir una entrada a
#' esta lista, no escribir un script.
#'
#' Campos obligatorios:
#'   id            identificador corto; da nombre a la carpeta de resultados
#'   organismo     nombre del organismo, para las tablas
#'   descripcion   qué se compara, en una línea
#'   rutas         directorio de FASTQ y ficheros de referencia
#'   read_type     "pe" o "se"
#'   muestras      vector con nombres = id de muestra, valores = condición
#'   contraste     list(num = , den = ) — el numerador va primero en el log2FC
#'   motores       cuáles de los tres probar
#'   pipelines     rutas a ejecutar: bowtie2, subjunc (splice-aware), salmon, kallisto
#'
#' Campos opcionales:
#'   batch         nombre de la covariable de bloqueo, si el diseño es pareado
#'   covariables   data.frame extra para el samplesheet (donante, lote...)
#'   publicado     tabla de referencia contra la que comparar
#'   enriquecimiento  parámetros del análisis funcional
#'   articulo      cita y PMID, para la tabla de procedencia
#'   conteos       matriz ya calculada, cuando no hay FASTQ que procesar
#'
#' Sobre `publicado`: `col_id` es la columna con el identificador que casa con
#' los nuestros. Si la tabla trae símbolos y nuestra matriz locus tags, hay que
#' indicar `traducir_con` (ruta al GTF) para que el harness resuelva el cruce.

DATASETS <- list()

# ── GSE273773 ───────────────────────────────────────────────────────────────
DATASETS$GSE273773 <- list(
  id = "GSE273773",
  organismo = "Escherichia coli K-12 BW25113",
  descripcion = "Nitrato de galio 1,25 mM frente a control, 10 h",
  read_type = "pe",
  rutas = list(
    fastq       = "~/Desktop/TFM_ejecuciones/GSE273773/fastq",
    genoma      = "~/Desktop/TFM_ejecuciones/GSE273773/referencia/genoma.fna",
    transcriptoma = "~/Desktop/TFM_ejecuciones/GSE273773/referencia/transcriptoma_locustag.fna",
    anotacion   = "~/Desktop/TFM_ejecuciones/GSE273773/referencia/anotacion.gtf"
  ),
  # Correspondencia verificada contra los metadatos de SRA: el orden de los SRR
  # está invertido respecto al de las muestras del GEO.
  muestras = c(SRR30100217 = "Control", SRR30100216 = "Control", SRR30100215 = "Control",
               SRR30100214 = "Galio",   SRR30100213 = "Galio",   SRR30100212 = "Galio"),
  contraste = list(num = "Galio", den = "Control"),
  # El galio es procariota y no tiene intrones, así que las dos rutas son
  # válidas y comparables. En un eucariota con intrones, bowtie2 perdería las
  # lecturas de unión exón-exón y la comparación no sería legítima.
  # Procariota sin intrones: las cuatro rutas son válidas y comparables entre
  # sí. subjunc es además el mismo motor que usó Rsubread en el artículo, de
  # modo que permite separar el efecto del alineador del resto del pipeline.
  pipelines = c("bowtie2", "subjunc", "salmon", "kallisto"),
  # Ejecuciones ya realizadas que se reutilizan en lugar de repetirse. Procesar
  # de nuevo las cuatro rutas costaría media hora y no cambiaría nada.
  pipeline_existente = list(
    bowtie2 = "outputs/GSE273773_bowtie2",
    subjunc = "outputs/GSE273773_subjunc",
    salmon  = "outputs/GSE273773_salmon"
  ),
  motores = c("DESeq2", "edgeR", "limma-voom"),
  publicado = list(
    tabla   = "~/Desktop/TFM_ejecuciones/GSE273773/publicado/deseq_Ga_publicado.csv.gz",
    sep     = ",",
    col_id  = "gene_ID",       # locus_tag: casa directamente con nuestra matriz
    col_lfc = "log2FoldChange",
    col_padj = "padj",
    fdr = 0.05, lfc = 1,
    # Cifras que declara el artículo, para comprobar que interpretamos su
    # criterio antes de comparar nada.
    deg_declarados = c(up = 581, down = 791),
    # Matriz de conteos de los autores. Permite la descomposición del error:
    # corriendo el mismo análisis sobre SUS conteos se separa lo que aporta la
    # etapa estadística de lo que aporta el pipeline.
    conteos = "~/Desktop/TFM_ejecuciones/GSE273773/publicado/counts_publicados.txt.gz",
    # Sus columnas van con nombre de grupo, no con el SRR.
    muestras_publicadas = c(SRR30100217="Control1", SRR30100216="Control2",
                            SRR30100215="Control3", SRR30100214="Gallium1",
                            SRR30100213="Gallium2", SRR30100212="Gallium3"),
    # Y sus filas con símbolo, mientras la tabla de resultados usa locus_tag:
    # el mapa se construye de la propia tabla publicada, que trae los dos.
    mapa_desde_tabla = c(de = "gene", a = "gene_ID")
  ),
  enriquecimiento = list(
    orgdb = "org.EcK12.eg.db", keytype = "SYMBOL",
    traducir_con = "~/Desktop/TFM_ejecuciones/GSE273773/referencia/anotacion.gtf",
    traducir_a = "gene",
    kegg_organismo = "eco", kegg_keytype = "ncbi-geneid",
    reactome = FALSE,   # Reactome no cubre procariotas
    colecciones = c("BP", "KEGG")
  ),
  articulo = list(
    cita = "Salazar-Aleman & Turner (2025) Sci Rep 15:1389",
    pmid = "39789098", doi = "10.1038/s41598-025-85772-y",
    herramientas = "rfastp 1.10, Rsubread 2.14.2, DESeq2 1.40.2"
  )
)

# ── GSE52778 (airway) ───────────────────────────────────────────────────────
# Solo conteos: los FASTQ son ~25 GB de RNA-seq humano y no aportan nada que no
# cubra ya el pipeline. Lo que aporta este dataset es Reactome (no cubre
# procariotas) y un diseño pareado por donante.
DATASETS$GSE52778 <- list(
  id = "GSE52778",
  organismo = "Homo sapiens (músculo liso de vía aérea)",
  descripcion = "Dexametasona 1 uM 18 h frente a sin tratar, 4 donantes",
  read_type = "pe",
  conteos = "~/Desktop/TFM_ejecuciones/GSE52778/conteos/count_matrix.tsv",
  muestras = c(SRR1039508 = "untrt", SRR1039509 = "trt",
               SRR1039512 = "untrt", SRR1039513 = "trt",
               SRR1039516 = "untrt", SRR1039517 = "trt",
               SRR1039520 = "untrt", SRR1039521 = "trt"),
  # Diseño pareado: cada donante aporta una muestra tratada y una sin tratar.
  # Es lo que permite ejercitar el bloqueo por covariable, que con réplicas
  # independientes no se puede enseñar.
  covariables = data.frame(
    sample_id = c("SRR1039508","SRR1039509","SRR1039512","SRR1039513",
                  "SRR1039516","SRR1039517","SRR1039520","SRR1039521"),
    donante = c("N61311","N61311","N052611","N052611",
                "N080611","N080611","N061011","N061011"),
    stringsAsFactors = FALSE),
  batch = "donante",
  contraste = list(num = "trt", den = "untrt"),
  pipelines = character(0),          # sin FASTQ: se entra por la matriz
  motores = c("DESeq2", "edgeR", "limma-voom"),
  publicado = NULL,                  # sin tabla de DEG publicada reutilizable
  enriquecimiento = list(
    orgdb = "org.Hs.eg.db", keytype = "ENSEMBL",
    kegg_organismo = "hsa", kegg_keytype = "ncbi-geneid",
    reactome = TRUE, reactome_organismo = "human",
    colecciones = c("BP", "KEGG", "REACTOME")
  ),
  # Genes cuya inducción por dexametasona está bien establecida. Sirven de
  # control positivo: si no salen, algo está mal configurado.
  marcadores = c("DUSP1", "KLF15", "PER1", "TSC22D3", "CRISPLD2"),
  articulo = list(
    cita = "Himes et al. (2014) PLoS ONE 9(6):e99625",
    pmid = "24926665", doi = "10.1371/journal.pone.0099625",
    herramientas = "TopHat, Cufflinks, DESeq"
  )
)

#' Devuelve la configuración de un dataset, con las rutas ya expandidas.
dataset_config <- function(id) {
  d <- DATASETS[[id]]
  if (is.null(d)) stop("Dataset no configurado: ", id,
                       ". Disponibles: ", paste(names(DATASETS), collapse = ", "))
  expandir <- function(x) if (is.character(x)) path.expand(x) else x
  d$rutas <- lapply(d$rutas %||% list(), expandir)
  if (!is.null(d$conteos)) d$conteos <- path.expand(d$conteos)
  if (!is.null(d$pipeline_existente))
    d$pipeline_existente <- lapply(d$pipeline_existente, path.expand)
  for (k in c("tabla")) if (!is.null(d$publicado[[k]]))
    d$publicado[[k]] <- path.expand(d$publicado[[k]])
  if (!is.null(d$enriquecimiento$traducir_con))
    d$enriquecimiento$traducir_con <- path.expand(d$enriquecimiento$traducir_con)
  d
}

datasets_disponibles <- function() names(DATASETS)
