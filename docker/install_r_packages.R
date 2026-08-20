#' install_r_packages.R
#' Instalacion de las dependencias de R dentro de la imagen.
#'
#' Las versiones se fijan por DOS vias complementarias, y ninguna sobra:
#'
#'   - Bioconductor se fija por RELEASE. Una release de Bioconductor es un
#'     conjunto de paquetes probados entre si contra una version concreta de R;
#'     pedir la release es mas fuerte que pedir versiones sueltas, porque
#'     garantiza tambien la compatibilidad cruzada. La release la fija la imagen
#'     base (bioconductor/bioconductor_docker:RELEASE_x_yy).
#'   - CRAN se fija por FECHA, apuntando a una instantanea del repositorio de
#'     Posit. Sin esto, `install.packages()` traeria la ultima version de shiny o
#'     bslib publicada el dia de la construccion, que es justo la variabilidad
#'     que la imagen debe eliminar.
#'
#' La instalacion falla si falta un paquete REQUERIDO. Es deliberado: una imagen
#' a la que le falta DESeq2 no debe construirse "con avisos", debe no
#' construirse.

cran_snapshot <- Sys.getenv("CRAN_SNAPSHOT", "")
if (nzchar(cran_snapshot)) {
  options(repos = c(CRAN = cran_snapshot))
  message("Repositorio CRAN fijado a: ", cran_snapshot)
}
options(Ncpus = max(1L, as.integer(Sys.getenv("BUILD_CPUS", "4"))))

requeridos_cran <- c(
  "shiny", "bslib", "shinyjs", "shinyFiles", "DT", "ggplot2", "plotly",
  "pheatmap", "ggrepel", "RColorBrewer", "processx", "withr", "htmlwidgets",
  "testthat"
)
requeridos_bioc <- c(
  "DESeq2", "edgeR", "limma", "tximport", "clusterProfiler", "apeglm", "fgsea",
  "AnnotationDbi", "GOSemSim", "enrichplot"
)
# Opcionales de la app: amplian lo que puede hacer, pero su ausencia esta
# contemplada en el codigo (la interfaz oculta la opcion en lugar de fallar).
# Se instalan igualmente en la imagen: el objetivo de contenerizar es que el
# entorno no dependa de lo que cada usuario tenga a mano.
opcionales <- c(
  "ashr", "IHW", "sva", "qvalue", "rtracklayer", "dearseq", "fishpond",
  "RNASeqPower", "ReactomePA", "reactome.db",
  "org.EcK12.eg.db", "org.Hs.eg.db"
)

instalar <- function(pkgs, obligatorio) {
  faltan <- setdiff(pkgs, rownames(installed.packages()))
  if (!length(faltan)) return(invisible(NULL))
  message("Instalando: ", paste(faltan, collapse = ", "))
  BiocManager::install(faltan, ask = FALSE, update = FALSE)
  sin_instalar <- setdiff(faltan, rownames(installed.packages()))
  if (length(sin_instalar)) {
    msg <- paste("No se pudieron instalar:", paste(sin_instalar, collapse = ", "))
    if (obligatorio) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(NULL)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
message("Bioconductor ", as.character(BiocManager::version()),
        " sobre R ", paste0(R.version$major, ".", R.version$minor))

instalar(requeridos_cran, obligatorio = TRUE)
instalar(requeridos_bioc, obligatorio = TRUE)
instalar(opcionales,      obligatorio = FALSE)

# Registro de lo que ha quedado instalado de verdad. Es el artefacto de
# procedencia del lado de R: lo que el informe de un analisis declara a
# posteriori, aqui queda fijado en la propia imagen.
ip <- installed.packages()
todos <- sort(unique(c(requeridos_cran, requeridos_bioc, opcionales)))
presentes <- intersect(todos, rownames(ip))
writeLines(
  c(paste0("R\t", paste0(R.version$major, ".", R.version$minor)),
    paste0("Bioconductor\t", as.character(BiocManager::version())),
    paste0(presentes, "\t", ip[presentes, "Version"])),
  "/opt/env/versions_r.txt"
)
message("Paquetes de R instalados: ", length(presentes), " de ", length(todos))
