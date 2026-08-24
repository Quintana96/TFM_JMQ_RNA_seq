#' ui_helpers.R
#' Helpers de UI compartidos: cabeceras de card con boton de descarga,
#' y fabricas de downloadHandler (CSV / plotly) que se inicializan en server.
#'
#' Nota: csv_download() y plotly_download() devuelven un downloadHandler y por tanto
#' deben usarse dentro de la función server. Se definen aquí para que esten disponibles
#' tanto en el modulo server_tab_results como en cualquier extensión futura.

#' Boton de descarga compacto para la cabecera de una tarjeta.
header_download_btn <- function(output_id, title) {
  downloadButton(
    output_id,
    label = NULL,
    icon  = icon("download"),
    class = "btn-sm btn-outline-secondary header-download",
    title = paste("Descargar", title)
  )
}

#' Card header con título a la izquierda y controles a la derecha.
#'
#' Sustituye a las diez repeticiones de
#' `card_header(tags$div(class = "card-title-download", tags$span(título), div(...)))`
#' que había esparcidas por ui_tab_deg.R y ui_tab_results.R. Cada copia
#' separaba los controles con su propio `gap` y su propio `align-items`, así que
#' la misma cabecera se alineaba distinto en cada pestana. Aquí la alineacion se
#' decide una vez.
#'
#' @param title Título de la tarjeta.
#' @param ... Controles que van a la derecha (selectores, botones).
#' @param download_id Si se indica, añade el boton de descarga al final.
card_header_tools <- function(title, ..., download_id = NULL) {
  tools <- list(...)
  if (!is.null(download_id)) tools <- c(tools, list(header_download_btn(download_id, title)))
  card_header(
    tags$div(
      class = "card-title-download",
      tags$span(title),
      if (length(tools)) tags$div(class = "card-header-tools", tools) else NULL
    )
  )
}

#' Card header con unicamente el boton de descarga a la derecha.
download_header <- function(title, output_id) {
  card_header_tools(title, download_id = output_id)
}

#' Pildora de estado: color + texto, nunca color solo.
#'
#' @param label Texto visible.
#' @param level Uno de "ok", "warn", "bad", "info", "neutral".
#' @param aria Etiqueta accesible completa (por defecto, `label`).
status_pill <- function(label, level = "neutral", aria = NULL) {
  level <- match.arg(level, c("ok", "warn", "bad", "info", "neutral"))
  tags$span(class = paste0("pill pill-", level),
            `aria-label` = aria %||% label,
            label)
}

#' Metrica compacta. `level` colorea el filo izquierdo según el ESTADO, no
#' según la posición de la tarjeta en la fila.
stat_tile <- function(value, label, level = "neutral") {
  level <- match.arg(level, c("ok", "warn", "bad", "info", "neutral"))
  clase <- if (identical(level, "neutral")) "metric-card" else paste0("metric-card is-", level)
  tags$div(class = clase,
    tags$div(class = "metric-card-value", value),
    tags$div(class = "metric-card-label", label)
  )
}

#' Título de sección dentro de una pestana (con subtitulo opcional).
section_title <- function(title, subtitle = NULL) {
  tagList(
    tags$div(class = "section-title", title),
    if (!is.null(subtitle)) tags$div(class = "section-subtitle", subtitle) else NULL
  )
}

#' Fabrica de downloadHandler que escribe un CSV a partir de data_fun()
csv_download <- function(prefix, data_fun) {
  downloadHandler(
    filename = function() paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(data_fun(), f, row.names = FALSE)
  )
}

#' Fabrica de downloadHandler que empaqueta un widget plotly (HTML + PNG si webshot2 disponible)
plotly_download <- function(prefix, plot_fun) {
  downloadHandler(
    filename = function() paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(f) {
      if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
        stop("El paquete htmlwidgets es necesario para descargar figuras Plotly.")
      }
      tmp_dir <- tempfile(prefix)
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
      html_file <- file.path(tmp_dir, "index.html")
      htmlwidgets::saveWidget(
        plot_fun(),
        file = html_file,
        selfcontained = FALSE,
        libdir = "lib"
      )

      image_files <- character(0)
      image_error <- NULL
      if (requireNamespace("webshot2", quietly = TRUE)) {
        png_file <- file.path(tmp_dir, "plot.png")
        image_error <- tryCatch({
          webshot2::webshot(html_file, png_file, vwidth = 1400, vheight = 900, zoom = 1.4)
          NULL
        }, error = function(e) conditionMessage(e))

        if (file.exists(png_file) && file.info(png_file)$size > 0) {
          image_files <- c(image_files, "plot.png")
        } else {
          jpg_file <- file.path(tmp_dir, "plot.jpg")
          image_error <- tryCatch({
            webshot2::webshot(html_file, jpg_file, vwidth = 1400, vheight = 900, zoom = 1.4)
            NULL
          }, error = function(e) conditionMessage(e))
          if (file.exists(jpg_file) && file.info(jpg_file)$size > 0) {
            image_files <- c(image_files, "plot.jpg")
          }
        }
      } else {
        image_error <- "El paquete webshot2 no está instalado."
      }

      if (!length(image_files)) {
        writeLines(
          c(
            "No se pudo generar una imagen PNG/JPEG de la figura.",
            "La figura interactiva está disponible en index.html.",
            "",
            "Para generar imagenes, instala o configura Chrome/Chromium para webshot2/chromote.",
            paste("Detalle:", image_error %||% "sin detalle disponible")
          ),
          file.path(tmp_dir, "README_imagen.txt")
        )
      }

      old_wd <- setwd(tmp_dir)
      on.exit(setwd(old_wd), add = TRUE)
      zip_files <- c("index.html", "lib", image_files)
      if (!length(image_files)) zip_files <- c(zip_files, "README_imagen.txt")
      utils::zip(zipfile = f, files = zip_files, flags = "-r9X")
    }
  )
}
