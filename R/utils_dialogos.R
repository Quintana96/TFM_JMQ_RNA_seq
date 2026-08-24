#' utils_dialogos.R
#' Dialogo nativo del sistema para elegir una carpeta.
#'
#' El selector de shinyFiles es un widget web: se abre anclado en una raíz, el
#' desplegable para cambiarla es poco visible y sus textos van en ingles dentro
#' de una aplicación en castellano. La consecuencia práctica es que parece que
#' no se puede salir de la carpeta del proyecto, aunque tecnicamente si se
#' pueda.
#'
#' Como SARA se ejecuta en la maquina del usuario —el servidor de Shiny y el
#' navegador son el mismo equipo—, se puede abrir el dialogo del sistema desde
#' el lado servidor y recuperar una RUTA de verdad. Es lo que hace la diferencia
#' con `fileInput`, que si abre el dialogo nativo pero SUBE una copia de los
#' ficheros: con 1,5 GB de FASTQ eso significa duplicarlos en un temporal, y
#' además el pipeline necesita un directorio, no una lista de ficheros subidos.
#'
#' Advertencia: el dialogo BLOQUEA el proceso de R mientras está abierto, de
#' modo que la aplicación queda congelada hasta que se elige o se cancela. En
#' una aplicación local de un solo usuario es aceptable; en un despliegue
#' compartido no lo sería, y por eso solo se ofrece cuando se detecta que la
#' ejecución es local.

#' Orden de preferencia por plataforma. En macOS `osascript` viene de serie; en
#' Linux se prueban los dos dialogos habituales de escritorio.
#' tcltk queda descartado a propósito: en macOS necesita XQuartz, que no viene
#' instalado, y ofrecer una opción que falla es peor que no ofrecerla.
dialogo_nativo_disponible <- function(contenedor = en_contenedor(),
                                      grafico = hay_entorno_grafico()) {
  if (isTRUE(contenedor)) return(FALSE)
  so <- Sys.info()[["sysname"]]
  if (so == "Darwin") return(nzchar(Sys.which("osascript")))
  if (.Platform$OS.type == "windows") return(TRUE)
  if (!isTRUE(grafico)) return(FALSE)
  nzchar(Sys.which("zenity")) || nzchar(Sys.which("kdialog"))
}

#' Abre el dialogo del sistema para elegir una carpeta.
#'
#' @param título texto de la ventana.
#' @param inicio carpeta en la que abrirse; si no existe se ignora.
#' @param ejecutar inyectable para poder probar sin abrir ventanas.
#' @return la ruta elegida, o NULL si se cancelo o no se pudo abrir. Cancelar NO
#'   es un error: es la respuesta más comun después de elegir, y tratarla como
#'   fallo llenaria la interfaz de avisos por usar el boton de cerrar.
elegir_directorio_nativo <- function(titulo = "Selecciona una carpeta",
                                     inicio = path.expand("~"),
                                     ejecutar = system2) {
  so <- Sys.info()[["sysname"]]
  if (!nzchar(inicio %||% "") || !dir.exists(inicio)) inicio <- path.expand("~")

  salida <- tryCatch({
    if (so == "Darwin") {
      # `activate` trae el dialogo al frente. Sin el se abre DETRÁS de la
      # ventana de la aplicación, y como R queda bloqueado el usuario ve una
      # aplicación congelada sin ninguna ventana que explique por qué.
      guion <- sprintf(
        'activate\nPOSIX path of (choose folder with prompt "%s" default location POSIX file "%s")',
        gsub('"', "'", titulo, fixed = TRUE),
        gsub('"', "", inicio, fixed = TRUE))
      ejecutar("osascript", c("-e", shQuote(guion)), stdout = TRUE, stderr = TRUE)
    } else if (.Platform$OS.type == "windows") {
      d <- utils::choose.dir(default = inicio, caption = titulo)
      if (is.na(d)) character(0) else d
    } else if (nzchar(Sys.which("zenity"))) {
      ejecutar("zenity", c("--file-selection", "--directory",
                           paste0("--title=", shQuote(titulo)),
                           paste0("--filename=", shQuote(paste0(inicio, "/")))),
               stdout = TRUE, stderr = FALSE)
    } else {
      ejecutar("kdialog", c("--getexistingdirectory", shQuote(inicio),
                            "--title", shQuote(titulo)),
               stdout = TRUE, stderr = FALSE)
    }
  }, error = function(e) character(0), warning = function(w) character(0))

  estado <- attr(salida, "status") %||% 0L
  # Cancelar devuelve un código distinto de cero (y en macOS el texto "User
  # canceled"): no hay ruta, y tampoco hay nada que reportar.
  if (!identical(as.integer(estado), 0L)) return(NULL)

  ruta <- trimws(paste(salida, collapse = ""))
  if (!nzchar(ruta)) return(NULL)
  # macOS devuelve las carpetas con barra final; normalizePath la quita y de
  # paso resuelve enlaces simbolicos, que es la ruta que el pipeline vera.
  ruta <- sub("/+$", "", ruta)
  if (!dir.exists(ruta)) return(NULL)
  normalizePath(ruta, mustWork = FALSE)
}

#' Válida una ruta escrita o pegada a mano.
#'
#' @return list(ruta, error). `ruta` vacia cuando no vale, con el motivo en
#'   `error`. Se admite `~` porque es lo que la gente teclea.
validar_directorio <- function(texto) {
  t <- trimws(texto %||% "")
  if (!nzchar(t)) return(list(ruta = "", error = NULL))
  # Al arrastrar una carpeta al terminal se pegan las rutas entrecomilladas.
  t <- gsub('^["\']|["\']$', "", t)
  t <- path.expand(t)
  if (!file.exists(t)) {
    return(list(ruta = "", error = "Esa ruta no existe."))
  }
  if (!dir.exists(t)) {
    return(list(ruta = "", error = "Eso es un fichero, no una carpeta."))
  }
  list(ruta = normalizePath(t, mustWork = FALSE), error = NULL)
}
