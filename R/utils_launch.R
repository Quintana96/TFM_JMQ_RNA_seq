#' utils_launch.R
#' Apertura de la interfaz en una ventana sin adornos de navegador.
#'
#' Shiny es una aplicación web: no existe una ventana nativa de verdad sin
#' envolverla en Electron, que obliga a empaquetar R entero y convierte el
#' arranque en un proceso de compilacion. El término medio es el modo
#' aplicación de los navegadores basados en Chromium: `--app=URL` abre una
#' ventana SIN barra de direcciones, pestanas, marcadores ni menu de extensiones.
#' Visualmente es indistinguible de una ventana propia; por dentro sigue siendo
#' el mismo navegador, con sus herramientas de desarrollo si se necesitan.
#'
#' Se invoca el BINARIO directamente en lugar de `open -a` (macOS): `open -a`
#' ignora los argumentos cuando la aplicación ya está abierta, de modo que con
#' el navegador en marcha se limitaria a traerlo al frente sin abrir nada. El
#' binario, en cambio, detecta la instancia existente y le pide la ventana, así
#' que no se duplica el proceso.

#' Candidatos por plataforma, en orden de preferencia.
#' Chrome primero por ser el más extendido; el resto comparten motor y aceptan
#' exactamente los mismos parámetros.
navegadores_de_aplicacion <- function() {
  if (Sys.info()[["sysname"]] == "Darwin") {
    apps <- c("Google Chrome", "Microsoft Edge", "Brave Browser", "Chromium", "Vivaldi")
    rutas <- unlist(lapply(apps, function(a) c(
      file.path("/Applications", paste0(a, ".app"), "Contents/MacOS", a),
      path.expand(file.path("~/Applications", paste0(a, ".app"), "Contents/MacOS", a)))))
    return(rutas)
  }
  if (.Platform$OS.type == "windows") {
    bases <- unique(c(Sys.getenv("PROGRAMFILES"), Sys.getenv("PROGRAMFILES(X86)"),
                      Sys.getenv("LOCALAPPDATA")))
    bases <- bases[nzchar(bases)]
    rel <- c("Google/Chrome/Application/chrome.exe",
             "Microsoft/Edge/Application/msedge.exe",
             "BraveSoftware/Brave-Browser/Application/brave.exe")
    return(as.vector(outer(bases, rel, function(b, r) file.path(b, r))))
  }
  # Linux: por PATH, que es donde los deja cualquier gestor de paquetes.
  bin <- c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
           "microsoft-edge", "brave-browser")
  rutas <- vapply(bin, function(b) unname(Sys.which(b)), character(1))
  unname(rutas[nzchar(rutas)])
}

#' Primer navegador disponible, o NULL.
buscar_navegador_de_aplicacion <- function(rutas = navegadores_de_aplicacion()) {
  for (r in rutas) if (nzchar(r) && file.exists(r)) return(r)
  NULL
}

#' Si el proceso corre dentro de un contenedor.
#'
#' Importa porque ahí no hay escritorio ni navegador: intentar abrir una ventana
#' escribiria un error en el log de arranque sin que nadie pueda hacer nada.
en_contenedor <- function() {
  if (file.exists("/.dockerenv")) return(TRUE)
  # Se comprueba la existencia antes de leer: fuera de Linux no hay /proc, y
  # readLines() sobre un fichero que no existe emite un aviso que tryCatch() no
  # atrapa si solo maneja errores. Acababa en la consola en cada arranque.
  if (!file.exists("/proc/1/cgroup")) return(FALSE)
  cg <- tryCatch(readLines("/proc/1/cgroup", warn = FALSE),
                 error = function(e) character(0))
  any(grepl("docker|kubepods|containerd", cg))
}

#' Si tiene sentido intentar abrir una ventana en esta maquina.
#'
#' Se exige un entorno gráfico. En Linux sin DISPLAY ni WAYLAND_DISPLAY el
#' navegador no arranca; en macOS y Windows siempre lo hay cuando hay sesión.
hay_entorno_grafico <- function() {
  so <- Sys.info()[["sysname"]]
  if (so == "Darwin" || .Platform$OS.type == "windows") return(TRUE)
  nzchar(Sys.getenv("DISPLAY")) || nzchar(Sys.getenv("WAYLAND_DISPLAY"))
}

#' Modo de apertura pedido por el entorno: "app", "navegador" o "ninguno".
#'
#' Por defecto "app" en una maquina de escritorio y "ninguno" dentro de un
#' contenedor o sin entorno gráfico, que es donde abrir algo no puede funcionar.
#' Las dos comprobaciones del entorno entran como argumentos con valor por
#' defecto en lugar de llamarse dentro. Así la decisión se puede probar sin
#' montar un contenedor ni apagar el servidor gráfico, que es la única forma de
#' verificar la rama que más importa: la que evita abrir una ventana donde no
#' puede haberla.
modo_de_apertura <- function(contenedor = en_contenedor(),
                             grafico = hay_entorno_grafico()) {
  pedido <- tolower(trimws(Sys.getenv("SARA_UI", "")))
  if (pedido %in% c("app", "navegador", "browser", "ninguno", "none")) {
    if (pedido == "browser") pedido <- "navegador"
    if (pedido == "none") pedido <- "ninguno"
    return(pedido)
  }
  if (isTRUE(contenedor) || !isTRUE(grafico)) return("ninguno")
  "app"
}

#' Abre `url` en una ventana de aplicación.
#'
#' @return TRUE si se lanzo el navegador, FALSE si no había ninguno compatible.
abrir_ventana_de_aplicacion <- function(url, navegador = buscar_navegador_de_aplicacion(),
                                        tamaño = Sys.getenv("SARA_UI_SIZE", "1440,900"),
                                        ejecutar = system2) {
  if (is.null(navegador)) return(FALSE)
  args <- c(paste0("--app=", url), paste0("--window-size=", tamaño),
            "--no-first-run", "--no-default-browser-check")
  ok <- tryCatch({
    # wait = FALSE es imprescindible: runApp() sigue bloqueado sirviendo la
    # aplicación, así que si se esperase al navegador no habría servidor al que
    # conectarse y la ventana quedaría en blanco.
    ejecutar(navegador, args, wait = FALSE, stdout = FALSE, stderr = FALSE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  isTRUE(ok)
}

#' Función que se le pasa a `runApp(launch.browser = ...)`.
#'
#' Devuelve NULL cuando no hay que abrir nada, que es lo que `runApp` espera
#' para no hacer nada (`launch.browser = FALSE`).
lanzador_de_interfaz <- function(modo = modo_de_apertura(),
                                 abrir = abrir_ventana_de_aplicacion) {
  if (identical(modo, "ninguno")) return(FALSE)
  if (identical(modo, "navegador")) return(TRUE)
  function(url) {
    if (!isTRUE(abrir(url))) {
      # Sin navegador compatible se cae al de siempre en lugar de no abrir
      # nada: es preferible una ventana con barra de direcciones a ninguna.
      message("No se ha encontrado Chrome, Edge, Brave, Chromium ni Vivaldi. ",
              "Se abre en el navegador por defecto.")
      utils::browseURL(url)
    }
  }
}
