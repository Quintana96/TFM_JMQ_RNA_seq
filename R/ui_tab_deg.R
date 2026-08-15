#' ui_tab_deg.R
#' Tab 4: Expresion diferencial (DEG).
#' Tres motores (DESeq2 / edgeR / limma-voom), metadatos editables y filtros
#' interactivos. Se renderiza via uiOutput("tab_deg_content") desde server_tab_deg.

#' Aviso cuando faltan paquetes Bioconductor necesarios.
#' Se distingue entre los que desactivan un motor (bloqueante) y los que solo
#' desactivan una opcion concreta, cuyo control se oculta en vez de ofrecerse.
ui_tab_deg_missing_pkgs <- function() {
  pkgs <- c()
  if (!isTRUE(HAS_DESEQ2))  pkgs <- c(pkgs, "DESeq2")
  if (!isTRUE(HAS_EDGER))   pkgs <- c(pkgs, "edgeR")
  if (!isTRUE(HAS_LIMMA))   pkgs <- c(pkgs, "limma")

  # clusterProfiler es requerido para el enriquecimiento, no opcional: sin el la
  # pestana entera queda inservible, asi que va en el aviso fuerte.
  if (!isTRUE(HAS_CLUSTERPROFILER)) pkgs <- c(pkgs, "clusterProfiler")

  optional <- c()
  if (!isTRUE(HAS_APEGLM) && !isTRUE(HAS_ASHR))
    optional <- c(optional, "apeglm (encogido de log2FC)")
  if (!isTRUE(HAS_FGSEA))  optional <- c(optional, "fgsea (GSEA)")
  if (!isTRUE(HAS_IHW))    optional <- c(optional, "IHW (ponderacion de hipotesis)")
  if (!isTRUE(HAS_QVALUE)) optional <- c(optional, "qvalue (segunda estimacion de pi0)")
  if (!isTRUE(HAS_ORGECDB)) optional <- c(optional, "org.EcK12.eg.db (GO de E. coli)")
  if (!isTRUE(HAS_PHEATMAP))
    optional <- c(optional, "pheatmap (heatmaps con dendrograma)")
  # rtracklayer sostiene tres cosas a la vez, y sin el fallan las tres en
  # silencio desde el punto de vista del usuario.
  if (!isTRUE(HAS_RTRACKLAYER))
    optional <- c(optional, paste("rtracklayer (mapa transcrito-gen para tximport,",
                                  "longitudes de gen y deteccion de rRNA)"))

  if (!length(pkgs) && !length(optional)) return(NULL)
  tagList(
    if (length(pkgs)) div(class = "alert alert-warning mt-3",
        icon("triangle-exclamation"),
        tags$b(" Paquetes ausentes: "), paste(pkgs, collapse = ", "),
        tags$p(class = "mb-0 small",
               "Instala los paquetes indicados para activar el motor correspondiente.")),
    if (length(optional)) div(class = "alert alert-secondary mt-3 py-2 small",
        icon("circle-info"),
        tags$b(" Funciones desactivadas por falta de paquete: "),
        paste(optional, collapse = "; "),
        tags$span(" Ejecuta requirements.sh para instalarlos."))
  )
}

#' UI del Tab 4 (cards de configuracion + zona de resultados).
#' Panel 1 del acordeon: de donde salen los conteos.
ui_deg_panel_datos <- function() {
  accordion_panel(
    "1 · Datos", value = "datos", icon = icon("table"),
    radioButtons("deg_source", label = NULL,
                 choices = c("Ejecucion actual" = "current",
                             "Ejecucion guardada" = "saved",
                             "Matriz subida" = "upload"),
                 selected = "current"),
    conditionalPanel(
      "input.deg_source === 'saved'",
      uiOutput("deg_saved_run_ui")
    ),
    conditionalPanel(
      "input.deg_source === 'upload'",
      fileInput("deg_counts_upload", "Matriz de conteos (TSV/CSV)",
                accept = c(".tsv", ".csv", ".txt"), width = "100%")
    ),
    tags$small(class = "text-muted",
               "Genes por filas, muestras por columnas. Se autodetectan los IDs de muestra.")
  )
}

#' Panel 2: samplesheet editable.
ui_deg_panel_metadatos <- function() {
  accordion_panel(
    "2 · Metadatos", value = "metadatos", icon = icon("list-check"),
    fileInput("deg_meta_upload", "Samplesheet (CSV/TSV)",
              accept = c(".csv", ".tsv", ".txt"), width = "100%"),
    div(class = "d-flex gap-2 flex-wrap mb-2",
        actionButton("deg_meta_sync_btn",
                     tagList(icon("rotate"), " Sincronizar"),
                     class = "btn-sm btn-outline-secondary",
                     title = "Sincronizar con las muestras detectadas"),
        actionButton("deg_meta_add_row_btn",
                     tagList(icon("plus"), " Anadir fila"),
                     class = "btn-sm btn-outline-secondary")
    ),
    tags$small(class = "text-muted d-block mb-2",
               "Columnas requeridas: sample_id, condition. Opcional: batch."),
    # Seudonimizacion: los identificadores de muestra de un estudio clinico
    # suelen ser identificativos y aqui viajan a graficos, informes y ficheros
    # persistidos.
    checkboxInput("deg_pseudonymize",
                  "Seudonimizar identificadores de muestra", FALSE),
    uiOutput("deg_identifying_cols"),
    DTOutput("deg_meta_table")
  )
}

#' Panel 3: diseno y contraste. El contraste se elige explicitamente (numerador
#' y denominador) en lugar de deducirse del nivel de referencia: con tres o mas
#' niveles, "nivel de referencia" dejaba sin decir cual de las comparaciones
#' posibles se estaba mostrando.
ui_deg_panel_diseno <- function() {
  accordion_panel(
    "3 · Diseno y contraste", value = "diseno", icon = icon("code-branch"),
    selectInput("deg_condition_col", "Columna de condicion",
                choices = c("condition"), selected = "condition"),
    layout_columns(
      col_widths = c(6, 6),
      selectInput("deg_contrast_num", "Numerador", choices = NULL),
      selectInput("deg_contrast_den", "Denominador", choices = NULL)
    ),
    tags$small(class = "text-muted d-block mb-2",
               "log2FC > 0 significa mayor expresion en el numerador."),
    checkboxInput("deg_use_batch", "Incluir variable batch", FALSE),
    conditionalPanel(
      "input.deg_use_batch === true && input.deg_advanced_design !== true",
      selectInput("deg_batch_col", "Columna de batch",
                  choices = c("batch"), selected = "batch")
    ),
    tags$hr(class = "my-2"),
    checkboxInput("deg_advanced_design", "Diseno avanzado (formula libre)", FALSE),
    conditionalPanel(
      "input.deg_advanced_design === true",
      textInput("deg_design_formula", NULL, value = "~ condition",
                placeholder = "~ subject + condition", width = "100%"),
      tags$small(class = "text-muted d-block mb-1",
                 paste("Permite disenos pareados (~ subject + condition),",
                       "covariables continuas (~ edad + condition) e",
                       "interacciones (~ genotipo * condition). Se valida antes",
                       "de ajustar.")),
      uiOutput("deg_design_feedback"),
      selectInput("deg_test_coef", "Coeficiente a testear", choices = NULL),
      tags$small(class = "text-muted d-block",
                 paste("Se rellena con los coeficientes del ultimo ajuste.",
                       "Dejalo en automatico para testear el contraste de la",
                       "condicion."))
    )
  )
}

#' Panel 4: todo lo que forma parte del TEST. Cambiar cualquier cosa de aqui
#' obliga a reajustar el modelo.
ui_deg_panel_analisis <- function() {
  accordion_panel(
    "4 · Test estadistico", value = "analisis", icon = icon("flask"),
    tags$small(class = "text-muted d-block mb-2",
               "Define el test: cambiar algo de aqui exige relanzar el analisis."),
    selectInput(
      "deg_method", "Motor",
      choices = list(
        "Parametricos" = as.list(stats::setNames(DEG_METHODS_PARAMETRIC,
                                                DEG_METHODS_PARAMETRIC)),
        # Los robustos solo tienen sentido con n grande; se ofrecen siempre
        # pero la app los sugiere unicamente cuando el tamano lo justifica.
        "Robustos (n grande)" = as.list(stats::setNames(
          c("Wilcoxon", if (isTRUE(HAS_DEARSEQ)) "dearseq"),
          c("Wilcoxon", if (isTRUE(HAS_DEARSEQ)) "dearseq"))),
        "Con incertidumbre de cuantificacion" = as.list(stats::setNames(
          if (isTRUE(HAS_FISHPOND)) "Swish" else character(0),
          if (isTRUE(HAS_FISHPOND)) "Swish" else character(0)))
      ),
      selected = "DESeq2"),
    uiOutput("deg_method_hint"),
    selectInput("deg_prefilter_mode", "Prefiltrado de genes",
                choices = c("Automatico (filterByExpr)" = "auto",
                            "Manual (umbral explicito)" = "manual"),
                selected = "auto"),
    conditionalPanel(
      "input.deg_prefilter_mode === 'manual'",
      layout_columns(
        col_widths = c(6, 6),
        numericInput("deg_min_count", "Min. cuenta/fila", value = 10, min = 0, step = 1),
        numericInput("deg_min_samples", "Min. muestras", value = NA, min = 1, step = 1)
      ),
      tags$small(class = "text-muted d-block mb-2",
                 "Min. muestras vacio = tamano del grupo mas pequeno.")
    ),
    sliderInput("deg_fdr_target", "FDR objetivo (alpha)",
                min = 0.01, max = 0.5, value = 0.05, step = 0.01),
    numericInput("deg_lfc_threshold", "Umbral |log2FC| del test",
                 value = 0, min = 0, max = 5, step = 0.25),
    tags$small(class = "text-muted d-block mb-2",
               paste("0 = test clasico (H0: log2FC = 0). Un valor > 0 testea",
                     "H0: |log2FC| <= umbral dentro del modelo",
                     "(lfcThreshold / glmTreat / treat), que es la forma de",
                     "exigir un fold-change minimo sin perder el control de la FDR.")),

    # Las opciones que solo aplican a DESeq2 van juntas y anunciadas como tales.
    # Antes se intercalaban con las generales y cada una repetia "Solo DESeq2"
    # en su propia nota al pie.
    tags$details(class = "mb-2",
      tags$summary(class = "small text-muted", "Opciones especificas de DESeq2"),
      tags$div(class = "mt-2",
        checkboxInput("deg_shrink", "Encoger log2FC para visualizacion (apeglm)", TRUE),
        tags$small(class = "text-muted d-block mb-2",
                   "Anade la columna log2FC_shrunk; no altera los p-valores."),
        if (isTRUE(HAS_IHW)) tagList(
          checkboxInput("deg_use_ihw", "Ponderar hipotesis con IHW en lugar de BH", FALSE),
          tags$small(class = "text-muted d-block mb-2",
                     paste("IHW pondera cada gen segun su baseMean y gana potencia",
                           "sin perder el control de la FDR, pero necesita muchos",
                           "tests para poder formar bins: con pocos genes se",
                           "reduce a BH."))
        ) else NULL,
        selectInput("deg_outliers", "Outliers de Cook",
                    choices = c(
                      "Excluirlos del test (por defecto)" = "na",
                      "Sustituir el valor atipico y volver a testear" = "refit",
                      "Ignorar el filtro: tratarlos como biologia real" = "keep"),
                    selected = "na", width = "100%"),
        tags$small(class = "text-muted d-block",
                   paste("Un gen con un valor extremo en una muestra queda sin",
                         "p-valor ajustado. Sustituirlo lo devuelve al test;",
                         "ignorar el filtro es lo apropiado cuando el valor",
                         "extremo ES el hallazgo (un gen que solo se expresa en",
                         "una muestra tratada) y no un artefacto."))
      )
    )
  )
}

#' Panel 5: variacion no deseada. La separacion entre las dos mitades es el
#' punto didactico: corregir-para-testear y corregir-para-visualizar no son la
#' misma operacion.
ui_deg_panel_variacion <- function() {
  accordion_panel(
    "5 · Variacion no deseada", value = "variacion", icon = icon("wave-square"),
    div(class = "alert alert-secondary py-2 px-2 small mb-2",
        icon("circle-info"),
        tags$b(" Dos cosas distintas."), " Para TESTEAR, la variacion se",
        " modela como covariable y los conteos se dejan intactos. Para",
        " VISUALIZAR (PCA, heatmap) hay que corregir la matriz, porque un",
        " grafico no puede incluir covariables. Corregir la matriz y luego",
        " testear sobre ella infla los falsos positivos."),
    tags$b(class = "small d-block", "En el modelo (afecta al test)"),
    if (isTRUE(HAS_SVA)) tagList(
      checkboxInput("deg_use_sva",
                    "Estimar variables sustitutas (sva) y anadirlas al diseno",
                    FALSE),
      conditionalPanel(
        "input.deg_use_sva === true",
        numericInput("deg_n_sv", "Numero de variables sustitutas (0 = automatico)",
                     value = 0, min = 0, max = 10, step = 1),
        tags$small(class = "text-muted d-block mb-2",
                   paste("Cada variable sustituta consume un grado de libertad.",
                         "Se reservan 3 g.l. residuales como minimo."))
      )
    ) else tags$small(class = "text-muted d-block mb-2",
                      "sva no esta instalado."),
    tags$hr(class = "my-2"),
    tags$b(class = "small d-block", "Solo en los graficos (no afecta al test)"),
    selectInput("deg_viz_correction", NULL,
                choices = c("Sin correccion" = "none",
                            "removeBatchEffect (limma)" = "rbe",
                            "ComBat-seq (sva)" = "combat"),
                selected = "none"),
    tags$small(class = "text-muted d-block",
               "Se aplica al PCA, al heatmap y a la matriz de distancias.")
  )
}

#' Panel 6: filtros que solo recortan lo que se ve. Deliberadamente separados
#' del panel 4 para que no se confundan con los umbrales del test.
ui_deg_panel_filtros <- function() {
  accordion_panel(
    "6 · Filtros de visualizacion", value = "filtros", icon = icon("eye"),
    div(class = "alert alert-secondary py-2 px-2 small mb-2",
        icon("eye"),
        tags$b(" Solo afectan a lo que se muestra."),
        " No cambian el modelo ni el control de la FDR. Para exigir un",
        " fold-change minimo con garantia estadistica usa el umbral del test",
        " en el panel 4."),
    sliderInput("deg_log2fc_cutoff", "|log2FC| minimo (visual)",
                min = 0, max = 5, value = 0, step = 0.1),
    sliderInput("deg_basemean_cutoff", "baseMean minimo (visual)",
                min = 0, max = 1000, value = 0, step = 5)
  )
}

#' UI del Tab 4.
#'
#' Antes era un grid rigido de seis tarjetas en tres columnas seguido de la zona
#' de resultados. Los seis bloques tienen tamanos muy distintos —el del test
#' ocupaba mas que los otros cinco juntos— asi que el grid quedaba desigual, y
#' habia que recorrer una pantalla entera de parametros antes de ver el primer
#' grafico, incluso para volver a mirar un resultado ya calculado.
#'
#' Ahora los parametros viven en una barra lateral plegable y los resultados
#' ocupan el area principal: el ciclo real de trabajo es ajustar un parametro y
#' mirar el efecto, no rellenar un formulario de una sola vez.
ui_tab_deg <- function() {
  tagList(
    ui_tab_deg_missing_pkgs(),
    layout_sidebar(
      sidebar = sidebar(
        title = "Parametros",
        width = 400,
        open = "desktop",
        # La accion principal va arriba del todo y fuera del acordeon: estaba
        # enterrada al final del cuarto bloque de seis, de modo que para
        # relanzar el analisis habia que desplegar y recorrer esa tarjeta.
        actionButton("run_deg_btn",
                     tagList(icon("flask"), " Lanzar DEG"),
                     class = "btn-primary btn-lg w-100"),
        tags$div(class = "mt-2",
                 verbatimTextOutput("deg_status_text", placeholder = TRUE)),
        accordion(
          ui_deg_panel_datos(),
          ui_deg_panel_metadatos(),
          ui_deg_panel_diseno(),
          ui_deg_panel_analisis(),
          ui_deg_panel_variacion(),
          ui_deg_panel_filtros(),
          open = c("datos", "diseno"),
          multiple = TRUE
        )
      ),
      # Zona de resultados (cargada al lanzar run_deg)
      uiOutput("deg_results_ui")
    )
  )
}

#' Estado inicial del area de resultados, antes del primer ajuste.
#' Antes esta zona quedaba sencillamente vacia y no habia nada que indicara que
#' faltaba pulsar "Lanzar DEG".
ui_tab_deg_placeholder <- function(msg = NULL) {
  div(class = "alert alert-secondary mt-2",
      icon("circle-info"),
      tags$b(" Sin resultados todavia. "),
      msg %||% paste("Revisa los parametros de la barra lateral y pulsa",
                     "\"Lanzar DEG\" para ajustar el modelo."))
}

#' Panel de diagnosticos post-ajuste.
#'
#' Pall et al. (PLoS Biology 2023) revisaron 4.616 datasets de GEO: solo el 25 %
#' tenia histogramas de p-valores con la forma esperada y el 37 % declaraba
#' implicitamente que mas de la mitad de los genes cambian. Son fallos invisibles
#' en la tabla de resultados. Se organizan en sub-pestanas para que cada
#' diagnostico se lea de uno en uno.
ui_tab_deg_diagnostics <- function() {
  navset_pill(
    id = "deg_diag_tabs",
    nav_panel(
      "p-valores",
      card(
        card_header_tools(
          "Distribucion de p-valores",
          selectInput("deg_diag_pv_subset", NULL,
                      choices = c("Genes que pasan el filtrado independiente" = "tested",
                                  "Todos los genes con p-valor" = "all"),
                      selected = "tested", width = "320px"),
          download_id = "download_deg_pvalue_hist"
        ),
        uiOutput("deg_diag_verdict"),
        plotly::plotlyOutput("deg_pvalue_hist", height = "340px")
      ),
      card(
        card_header("Formas de referencia"),
        tags$p(class = "small text-muted mb-1",
               paste("Comparalas con el histograma de arriba. La forma esperada es un",
                     "pico a la izquierda sobre un suelo aproximadamente plano; las",
                     "otras dos indican que algo no encaja en el modelo.")),
        plotOutput("deg_pvalue_reference", height = "220px")
      )
    ),
    nav_panel(
      "Dispersiones",
      card(
        download_header("Estimacion de la dispersion (equivalente a plotDispEsts)",
                        "download_deg_disp_plot"),
        tags$p(class = "small text-muted mb-1",
               paste("Las estimaciones por gen deben repartirse alrededor de la curva",
                     "ajustada y los valores finales acercarse a ella. Una nube que no",
                     "sigue la curva indica un ajuste de dispersion defectuoso.")),
        plotly::plotlyOutput("deg_disp_plot", height = "420px")
      )
    ),
    nav_panel(
      "Normalizacion (RLE)",
      card(
        download_header("Relative Log Expression por muestra", "download_deg_rle_plot"),
        tags$p(class = "small text-muted mb-1",
               paste("Cada muestra deberia tener la mediana cerca de 0 y un rango",
                     "estrecho. Una mediana desplazada senala fallo de normalizacion o",
                     "una muestra problematica.")),
        plotly::plotlyOutput("deg_rle_plot", height = "420px")
      )
    ),
    nav_panel(
      "Outliers (Cook)",
      card(
        download_header("Reparto de los maximos de distancia de Cook",
                        "download_deg_cooks_plot"),
        tags$p(class = "small text-muted mb-1",
               paste("Si una sola muestra concentra los outliers, el problema es de esa",
                     "muestra y no de los genes. La linea marca el reparto esperado por",
                     "azar (1/n).")),
        uiOutput("deg_cooks_warning"),
        plotly::plotlyOutput("deg_cooks_plot", height = "380px")
      )
    ),
    nav_panel(
      "Sesgo de longitud",
      card(
        download_header("Probabilidad de ser diferencial segun la longitud del gen",
                        "download_deg_lenbias_plot"),
        tags$p(class = "small text-muted mb-1",
               paste("El analisis de sobre-representacion (GO/KEGG) asume que todos",
                     "los genes tienen la misma probabilidad de ser detectados. En",
                     "RNA-seq no es cierto: los transcritos largos acumulan mas",
                     "lecturas y salen diferenciales con mas facilidad, lo que arrastra",
                     "a las categorias que los contienen. Esta curva mide si ocurre en",
                     "TUS datos: si es plana, corregir por longitud no aporta nada.")),
        uiOutput("deg_lenbias_verdict"),
        plotly::plotlyOutput("deg_lenbias_plot", height = "380px")
      )
    )
  )
}

#' tagList con el navset de resultados (tabla, plots, enriquecimiento).
ui_tab_deg_results <- function() {
  tagList(
    # Banner con el contraste realmente testeado. Es informacion imprescindible
    # para interpretar un volcano, y ademas es donde avisamos de que con >2
    # niveles de condition solo se esta mostrando una de las comparaciones.
    uiOutput("deg_contrast_banner"),
    # Once pestanas al mismo nivel no caben en una linea a 1440 px: se partian
    # en dos filas y la barra dejaba de leerse como una barra. Se agrupan por lo
    # que responden —que genes cambian, como se parecen las muestras, se sostiene
    # el ajuste, que significa— y dentro de cada grupo van pildoras.
    navset_tab(
      id = "deg_result_tabs",
      nav_panel(
        "Genes",
        navset_pill(
          id = "deg_genes_tabs",
          nav_panel(
        "Tabla",
        card(
          download_header("Tabla DEG (filtrada)", "download_deg_table"),
          DTOutput("deg_table"),
          # Desglose de por que hay genes sin padj. Colapsarlos a "no
          # significativo" oculta informacion util (un outlier de Cook puede ser
          # el hallazgo mas interesante, o la senal de que una muestra esta mal).
          uiOutput("deg_na_breakdown")
        )
      ),
      nav_panel(
        "Volcano",
        card(
          download_header("Volcano plot", "download_deg_volcano_plot"),
          plotly::plotlyOutput("deg_volcano_plot", height = "480px")
        )
      ),
      nav_panel(
        "MA",
        card(
          download_header("MA plot", "download_deg_ma_plot"),
          plotly::plotlyOutput("deg_ma_plot", height = "480px")
        )
      ),
        )
      ),
      nav_panel(
        "Muestras",
        navset_pill(
          id = "deg_samples_tabs",
          nav_panel(
        "PCA",
        card(
          card_header_tools(
            "PCA (vst/rlog)",
            selectInput("deg_pca_color", NULL, choices = c("condition"),
                        selected = "condition", width = "190px"),
            download_id = "download_deg_pca_plot"
          ),
          # Colorear por covariable es lo que convierte el PCA en un diagnostico:
          # con el color fijado en `condition` es imposible ver que el batch esta
          # confundido con la condicion, que es justo lo que hay que detectar.
          tags$small(class = "text-muted d-block mb-1",
                     paste("Cambia el color a una covariable (batch, lote, sujeto)",
                           "para comprobar si explica la separacion entre muestras",
                           "mejor que la condicion. Si lo hace, tienes un efecto",
                           "confundido con el contraste.")),
          plotly::plotlyOutput("deg_pca_plot", height = "460px")
        )
      ),
      nav_panel(
        "Heatmap top-N",
        card(
          card_header_tools(
            "Heatmap top-N genes",
            numericInput("deg_heatmap_topn", NULL, value = 30, min = 5, max = 200,
                         step = 5, width = "100px"),
            download_id = "download_deg_heatmap"
          ),
          plotOutput("deg_heatmap", height = "560px")
        )
      ),
      nav_panel(
        "Distancia entre muestras",
        card(
          download_header("Matriz de distancias entre muestras", "download_deg_dist_heatmap"),
          plotOutput("deg_dist_heatmap", height = "520px")
        )
      ),
        )
      ),
      nav_panel(
        "Diagnosticos",
        ui_tab_deg_diagnostics()
      ),
      nav_panel(
        "Robustez",
        navset_pill(
          id = "deg_robust_tabs",
          nav_panel(
        "Replicabilidad",
        card(
          card_header_tools(
            "Replicabilidad por bootstrap",
            numericInput("deg_boot_n", NULL, value = 20, min = 5, max = 100,
                         step = 5, width = "90px"),
            actionButton("deg_run_boot_btn", tagList(icon("dice"), " Estimar"),
                         class = "btn-sm btn-outline-secondary")
          ),
          tags$p(class = "small text-muted mb-1",
                 paste("Remuestrea las muestras, repite el analisis y mide si la",
                       "lista aguanta. Es el procedimiento que recomienda el estudio",
                       "de replicabilidad en cohortes pequenas: Spearman > 0,9 indica",
                       "precision alta y < 0,8 avisa de probables falsos positivos.",
                       "Cuesta un reajuste del modelo por remuestreo.")),
          uiOutput("deg_boot_verdict"),
          DTOutput("deg_boot_table"),
          plotly::plotlyOutput("deg_boot_plot", height = "300px")
        )
      ),
      nav_panel(
        "Comparar metodos",
        card(
          card_header_tools(
            "Solapamiento entre metodos",
            selectInput("deg_compare_method", NULL, choices = NULL, width = "180px"),
            actionButton("deg_run_compare_btn",
                         tagList(icon("code-compare"), " Comparar"),
                         class = "btn-sm btn-outline-secondary")
          ),
          tags$p(class = "small text-muted mb-1",
                 paste("Corre un segundo motor sobre los mismos datos y compara las",
                       "listas de significativos. Con n grande, la discrepancia entre",
                       "parametricos y robustos es informativa: es el punto de la",
                       "controversia sobre el control real de la FDR.")),
          uiOutput("deg_compare_summary"),
          plotly::plotlyOutput("deg_compare_plot", height = "300px")
        )
      ),
        )
      ),
      nav_panel(
        "Enriquecimiento",
        card(
          # Los siete selectores de este bloque vivian dentro de la cabecera de
          # la tarjeta, sin etiqueta ninguna, apoyados solo en que su valor por
          # defecto se leyera como titulo ("ORA (lista significativa)",
          # "GO: Procesos biologicos", "Metrica: stat"...). En cuanto se cambiaba
          # uno dejaba de estar claro que era, y los que aparecen y desaparecen
          # segun la ontologia hacian saltar el ancho de la cabecera. Ahora van
          # en el cuerpo, con etiqueta, y en la cabecera queda solo la accion.
          card_header_tools(
            "Enriquecimiento funcional",
            actionButton("deg_run_enrich_btn",
                         tagList(icon("play"), " Calcular"),
                         class = "btn-sm btn-primary")
          ),

          div(class = "border rounded p-2 mb-3",
            layout_columns(
              col_widths = c(4, 4, 4),
              # ORA parte de la lista umbralizada y por tanto es ciego a
              # senales debiles pero coordinadas; GSEA usa el ranking completo.
              selectInput("deg_enrich_approach", "Enfoque",
                          choices = c(c("ORA (lista significativa)" = "ora"),
                                      if (isTRUE(HAS_FGSEA))
                                        c("GSEA (ranking completo)" = "gsea")),
                          selected = "ora"),
              # GMT no es una ontologia mas: es la via para trabajar sin OrgDb
              # (organismos no modelo) y con colecciones curadas como MSigDB.
              selectInput("deg_ontology", "Coleccion",
                          choices = c("GO: Procesos biologicos" = "BP",
                                      "GO: Funcion molecular" = "MF",
                                      "GO: Componente celular" = "CC",
                                      "KEGG" = "KEGG",
                                      "Gene sets propios (GMT)" = "GMT"),
                          selected = "BP"),
              conditionalPanel(
                "input.deg_enrich_approach === 'gsea'",
                selectInput("deg_gsea_metric", "Metrica de ranking",
                            choices = c("stat" = "stat",
                                        "log2FC" = "log2FC",
                                        "signo x -log10(p)" = "signed_p"),
                            selected = "stat")
              )
            ),
            # keyType: los IDs de featureCounts sobre un GFF procariota son
            # locus tags, no simbolos. Fijarlo a SYMBOL hacia fallar el mapeo
            # en silencio, asi que ahora es explicito y seleccionable.
            # Organismo: la app tenia el OrgDb cableado a E. coli, asi que con
            # datos de cualquier otro organismo el enriquecimiento GO no podia
            # funcionar aunque su paquete estuviera instalado.
            conditionalPanel(
              "input.deg_ontology !== 'KEGG' && input.deg_ontology !== 'GMT'",
              layout_columns(
                col_widths = c(8, 4),
                selectInput("deg_orgdb", "Organismo (OrgDb)",
                            choices = if (length(ORGDBS_DISPONIBLES))
                              stats::setNames(ORGDBS_DISPONIBLES,
                                              vapply(ORGDBS_DISPONIBLES, orgdb_label,
                                                     character(1)))
                            else c("(sin OrgDb instalado)" = ""),
                            selected = if (length(ORGDBS_DISPONIBLES))
                              ORGDBS_DISPONIBLES[1] else ""),
                selectInput("deg_go_keytype", "Tipo de identificador",
                            choices = c("SYMBOL"), selected = "SYMBOL")
              )
            ),
            conditionalPanel(
              "input.deg_ontology === 'KEGG'",
              layout_columns(
                col_widths = c(6, 6),
                textInput("deg_kegg_organism", "Codigo de organismo KEGG",
                          value = "eco", placeholder = "eco, hsa, mmu..."),
                selectInput("deg_kegg_keytype", "Tipo de identificador",
                            choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"),
                            selected = "kegg")
              )
            )
          ),

          # Conjuntos propios en formato GMT. Con esto el enriquecimiento deja de
          # depender de que exista un OrgDb del organismo: basta con un fichero
          # de conjuntos cuyos identificadores coincidan con los de la matriz.
          conditionalPanel(
            "input.deg_ontology === 'GMT'",
            div(class = "border rounded p-2 mb-2",
                fileInput("deg_gmt_file", "Fichero .gmt de conjuntos de genes",
                          accept = c(".gmt", ".txt"), width = "100%"),
                tags$small(class = "text-muted d-block mb-1",
                           paste("Un conjunto por linea: nombre, descripcion y genes",
                                 "separados por tabuladores (MSigDB, regulones de",
                                 "RegulonDB, firmas propias). Los identificadores del",
                                 "fichero deben ser los mismos que los de la matriz de",
                                 "conteos; el selector de keyType no se aplica aqui.")),
                uiOutput("deg_gmt_summary"))
          ),

          # simplify() solo aplica al ORA sobre GO: colapsa terminos redundantes
          # por similitud semantica, y necesita el grafo de una sola ontologia.
          conditionalPanel(
            "input.deg_enrich_approach === 'ora' && input.deg_ontology !== 'KEGG' && input.deg_ontology !== 'GMT'",
            checkboxInput("deg_go_simplify",
                          "Colapsar terminos GO redundantes (simplify, similitud de Wang > 0,7)",
                          FALSE)
          ),
          conditionalPanel(
            "input.deg_enrich_approach === 'ora'",
            checkboxInput("deg_ora_directional",
                          "Analizar tambien por separado los genes al alza y a la baja",
                          FALSE),
            tags$small(class = "text-muted d-block mb-2",
                       paste("Mezclar las dos direcciones diluye la senal: una ruta con",
                             "la mitad de sus genes inducidos y la otra mitad reprimidos",
                             "da el mismo p-valor que una ruta coherentemente inducida.",
                             "Multiplica por tres el tiempo de calculo."))
          ),
          conditionalPanel(
            "input.deg_ontology !== 'KEGG' && input.deg_ontology !== 'GMT'",
            checkboxInput("deg_enrich_readable",
                          "Mostrar simbolos de gen en lugar de los IDs de entrada",
                          TRUE)
          ),

          # Parametros de GSEA. Estaban fijados a los valores por defecto de
          # run_gsea() y son justamente los que cambian la lectura del resultado.
          conditionalPanel(
            "input.deg_enrich_approach === 'gsea'",
            div(class = "border rounded p-2 mb-2",
                tags$b(class = "small d-block mb-1", "Parametros de GSEA"),
                layout_columns(
                  col_widths = c(3, 3, 3, 3),
                  numericInput("deg_gsea_min_size", "Tamano minimo del conjunto",
                               value = 10, min = 2, max = 500, step = 5),
                  numericInput("deg_gsea_max_size", "Tamano maximo",
                               value = 500, min = 10, max = 5000, step = 50),
                  numericInput("deg_gsea_pcutoff", "Corte de p ajustado",
                               value = 0.05, min = 0.001, max = 1, step = 0.01),
                  checkboxInput("deg_gsea_eps_exact", "P-valores exactos (eps = 0)",
                                FALSE)
                ),
                tags$small(class = "text-muted d-block",
                           paste("Los conjuntos muy pequenos dan NES inestables y los muy",
                                 "grandes son inespecificos; acotar el tamano tambien",
                                 "reduce el numero de tests y con ello el castigo del",
                                 "ajuste multiple. El corte de p filtra la tabla que",
                                 "devuelve clusterProfiler: ponlo a 1 para ver todos los",
                                 "conjuntos testeados, que es la unica forma de",
                                 "distinguir 'nada llega a 0,05' de un fallo del",
                                 "analisis. Por defecto fgsea trunca los p-valores en",
                                 "1e-10 y los conjuntos mas significativos empatan en ese",
                                 "valor; con eps = 0 se estiman exactos, a costa de",
                                 "bastante mas tiempo de calculo."))
            )
          ),

          # Sub-pestanas: el running score y la comparacion ORA/GSEA son lecturas
          # distintas del mismo calculo y apiladas en una sola vista no se leen.
          navset_pill(
            id = "deg_enrich_tabs",
            nav_panel(
              "Resultados",
              uiOutput("deg_enrich_mapping"),
              plotly::plotlyOutput("deg_enrich_dotplot", height = "440px"),
              tags$hr(),
              tags$div(
                class = "card-title-download mb-2",
                tags$span(class = "fw-semibold", "Tabla de terminos enriquecidos"),
                header_download_btn("download_enrich_table", "la tabla de enriquecimiento")
              ),
              DTOutput("deg_enrich_table")
            ),
            nav_panel(
              "Running score (GSEA)",
              tags$p(class = "small text-muted mb-1",
                     paste("La curva acumula el estadistico a lo largo del ranking:",
                           "sube al encontrar un gen del conjunto y baja con cada gen",
                           "que no lo es. El pico es el enrichment score, y los genes",
                           "que van desde el extremo del ranking hasta el pico son el",
                           "leading edge: el subconjunto que sostiene el resultado, no",
                           "el conjunto entero.")),
              div(style = "display:flex;gap:8px;align-items:flex-end;flex-wrap:wrap;",
                  selectInput("deg_gsea_term", "Conjunto a visualizar",
                              choices = NULL, width = "520px")),
              plotly::plotlyOutput("deg_gsea_runscore", height = "380px"),
              uiOutput("deg_gsea_leading")
            ),
            nav_panel(
              "ORA frente a GSEA",
              tags$p(class = "small text-muted mb-1",
                     paste("Corre los dos enfoques sobre el mismo contraste y la misma",
                           "ontologia y compara los terminos significativos. No compiten:",
                           "el ORA pregunta si la lista umbralizada esta enriquecida y",
                           "GSEA si el conjunto esta desplazado en el ranking completo.",
                           "Los terminos que solo ve GSEA son los de senal debil pero",
                           "coordinada que el corte de la lista deja fuera.")),
              div(style = "display:flex;gap:8px;align-items:flex-end;flex-wrap:wrap;",
                  actionButton("deg_run_enrich_compare_btn",
                               tagList(icon("code-compare"), " Comparar ORA y GSEA"),
                               class = "btn-sm"),
                  downloadButton("download_enrich_compare", "Descargar comparacion",
                                 icon = icon("download"),
                                 class = "btn-sm btn-outline-secondary")),
              uiOutput("deg_enrich_compare_summary"),
              plotly::plotlyOutput("deg_enrich_compare_plot", height = "320px"),
              DTOutput("deg_enrich_compare_table")
            )
          )
        )
      ),
      nav_panel(
        "Reproducibilidad",
        card(
          card_header("Informe y script equivalente"),
          tags$p(class = "small text-muted",
                 paste("Una app grafica no deja rastro de que se hizo. El informe",
                       "recoge todos los parametros, el contraste, la formula del",
                       "diseno, los diagnosticos y las versiones de los paquetes;",
                       "el script reproduce el mismo analisis con Bioconductor,",
                       "fuera de la app.")),
          div(style = "display:flex;gap:10px;flex-wrap:wrap;margin-top:6px;",
              downloadButton("download_deg_report", "Informe HTML",
                             icon = icon("file-code"), class = "btn-sm"),
              downloadButton("download_deg_script", "Script R equivalente",
                             icon = icon("r-project"), class = "btn-sm")),
          tags$hr(),
          # La correspondencia de la seudonimizacion se descarga aparte del
          # informe: exportar los identificadores reales debe ser una decision
          # deliberada, no un efecto colateral.
          uiOutput("deg_pseudonym_ui"),
          tags$b(class = "small", "Vista previa del script"),
          verbatimTextOutput("deg_script_preview", placeholder = TRUE)
        )
      )
    )
  )
}
