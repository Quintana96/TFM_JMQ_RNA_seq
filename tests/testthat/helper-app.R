#' helper-app.R
#' Carga la aplicacion para poder testear sus funciones puras.
#'
#' La app no es un paquete de R: Shiny sourcea `global.R` y `R/*.R` por su
#' cuenta al arrancar. testthat no hace eso, asi que se replica aqui el mismo
#' orden de carga. Este helper lo ejecuta testthat automaticamente (los ficheros
#' `helper-*.R` se cargan antes que los `test-*.R`).
#'
#' Se sourcean solo `global.R` y `R/`: `ui.R` y `server.R` construyen objetos
#' Shiny que no hacen falta para testear las funciones puras, que son las que
#' contienen la logica estadistica.

app_root <- local({
  d <- normalizePath(file.path(dirname(dirname(getwd()))), mustWork = FALSE)
  # Segun desde donde se invoque (raiz del proyecto o tests/testthat), la raiz
  # de la app esta a distinta profundidad. Se localiza por un fichero ancla.
  for (cand in c(getwd(), dirname(getwd()), d, file.path(d, ".."))) {
    if (file.exists(file.path(cand, "global.R"))) return(normalizePath(cand))
  }
  stop("No se encuentra la raiz de la aplicacion (global.R).")
})

suppressPackageStartupMessages(
  sys.source(file.path(app_root, "global.R"), envir = globalenv())
)
for (f in sort(list.files(file.path(app_root, "R"), pattern = "[.]R$",
                          full.names = TRUE))) {
  sys.source(f, envir = globalenv())
}

#' Matriz de conteos sintetica con señal diferencial conocida.
#'
#' Se genera con semilla fija para que los tests sean deterministas. `n_de`
#' genes tienen un efecto real entre grupos; el resto es ruido.
make_test_counts <- function(n_genes = 400, n_per_group = 4, n_de = 40,
                             lfc = 2, seed = 42) {
  withr::with_seed(seed, {
    n <- n_per_group * 2
    base <- stats::rnbinom(n_genes, mu = 200, size = 5) + 20
    m <- matrix(0L, nrow = n_genes, ncol = n)
    for (j in seq_len(n)) {
      mu <- base
      if (j > n_per_group) mu[seq_len(n_de)] <- mu[seq_len(n_de)] * 2^lfc
      m[, j] <- stats::rnbinom(n_genes, mu = mu, size = 10)
    }
    rownames(m) <- sprintf("gene%04d", seq_len(n_genes))
    colnames(m) <- c(sprintf("ctrl%d", seq_len(n_per_group)),
                     sprintf("trt%d", seq_len(n_per_group)))
    m
  })
}

#' Samplesheet acompañante de `make_test_counts()`.
#'
#' Los niveles se eligen a proposito de forma que el orden alfabetico NO
#' coincida con el contraste habitual: "trt" vs "ctrl" tiene como denominador
#' "ctrl", que es el primero alfabeticamente, mientras que "ctrl" vs "trt" tiene
#' como denominador "trt", que no lo es. Ese segundo caso es el que destapa los
#' errores de reproducir el contraste por orden alfabetico.
make_test_meta <- function(counts, batch = FALSE) {
  n <- ncol(counts)
  df <- data.frame(
    sample_id = colnames(counts),
    condition = ifelse(grepl("^ctrl", colnames(counts)), "ctrl", "trt"),
    stringsAsFactors = FALSE
  )
  if (batch) df$lote <- rep(c("A", "B"), length.out = n)
  rownames(df) <- df$sample_id
  df
}
