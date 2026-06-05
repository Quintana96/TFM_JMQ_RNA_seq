#' ui_helpers.R
#' Helpers de UI compartidos: cabeceras de card con boton de descarga,
#' y fabricas de downloadHandler (CSV / plotly) que se inicializan en server.
#'
#' Nota: csv_download() y plotly_download() devuelven un downloadHandler y por tanto
#' deben usarse dentro de la funcion server. Se definen aqui para que esten disponibles
#' tanto en el modulo server_tab_results como en cualquier extension futura.

#' Card header con boton de descarga a la derecha (icono download)
download_header <- function(title, output_id) {
  card_header(
    tags$div(
      class = "card-title-download",
      tags$span(title),
      downloadButton(
        output_id,
        label = NULL,
        icon = icon("download"),
        class = "btn-sm btn-outline-secondary header-download",
        title = paste("Descargar", title)
      )
    )
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
        image_error <- "El paquete webshot2 no esta instalado."
      }

      if (!length(image_files)) {
        writeLines(
          c(
            "No se pudo generar una imagen PNG/JPEG de la figura.",
            "La figura interactiva esta disponible en index.html.",
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
