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

  optional <- c()
  if (!isTRUE(HAS_APEGLM) && !isTRUE(HAS_ASHR))
    optional <- c(optional, "apeglm (encogido de log2FC)")
  if (!isTRUE(HAS_FGSEA))  optional <- c(optional, "fgsea (GSEA)")
  if (!isTRUE(HAS_IHW))    optional <- c(optional, "IHW (ponderacion de hipotesis)")
  if (!isTRUE(HAS_ORGECDB)) optional <- c(optional, "org.EcK12.eg.db (GO de E. coli)")

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
ui_tab_deg <- function() {
  tagList(
    ui_tab_deg_missing_pkgs(),
    tags$div(
      style = paste(
        "display:grid;",
        "grid-template-columns:repeat(3, minmax(280px, 1fr));",
        "grid-auto-rows:minmax(220px, auto);",
        "gap:12px; align-items:stretch; margin-top:6px;"
      ),

      # Card 1: Datos
      card(
        card_header("1. Datos"),
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
                    accept = c(".tsv", ".csv", ".txt"))
        ),
        tags$small(class = "text-muted",
                   "Genes por filas, muestras por columnas. Se autodetectan los IDs de muestra.")
      ),

      # Card 2: Metadatos
      card(
        card_header("2. Metadatos"),
        fileInput("deg_meta_upload", "Samplesheet (CSV/TSV)",
                  accept = c(".csv", ".tsv", ".txt")),
        div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-bottom:6px;",
            actionButton("deg_meta_sync_btn",
                         tagList(icon("rotate"), " Sincronizar con muestras detectadas"),
                         class = "btn-sm"),
            actionButton("deg_meta_add_row_btn",
                         tagList(icon("plus"), " Anadir fila"),
                         class = "btn-sm")
        ),
        tags$small(class = "text-muted",
                   "Columnas requeridas: sample_id, condition. Opcional: batch."),
        DTOutput("deg_meta_table")
      ),

      # Card 3: Diseno. El contraste se elige explicitamente (numerador y
      # denominador) en lugar de deducirse del nivel de referencia: con tres o
      # mas niveles, "nivel de referencia" dejaba sin decir cual de las
      # comparaciones posibles se estaba mostrando.
      card(
        card_header("3. Diseno y contraste"),
        selectInput("deg_condition_col", "Columna de condicion",
                    choices = c("condition"), selected = "condition"),
        layout_columns(
          col_widths = c(6, 6),
          selectInput("deg_contrast_num", "Numerador", choices = NULL),
          selectInput("deg_contrast_den", "Denominador (referencia)", choices = NULL)
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
                    placeholder = "~ subject + condition"),
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
      ),

      # Card 4: Analisis. Todo lo que hay aqui forma parte del TEST, asi que
      # cambiarlo obliga a reajustar el modelo (boton "Lanzar DEG").
      card(
        card_header("4. Analisis DEG (define el test)"),
        selectInput("deg_method", "Motor",
                    choices = c("DESeq2", "edgeR", "limma-voom"),
                    selected = "DESeq2"),
        selectInput("deg_prefilter_mode", "Prefiltrado de genes",
                    choices = c("Automatico (filterByExpr)" = "auto",
                                "Manual (umbral explicito)" = "manual"),
                    selected = "auto"),
        conditionalPanel(
          "input.deg_prefilter_mode === 'manual'",
          layout_columns(
            col_widths = c(6, 6),
            numericInput("deg_min_count", "Min. cuenta por fila", value = 10, min = 0, step = 1),
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
        checkboxInput("deg_shrink", "Encoger log2FC para visualizacion (apeglm)", TRUE),
        tags$small(class = "text-muted d-block mb-2",
                   "Solo DESeq2. Anade la columna log2FC_shrunk; no altera los p-valores."),
        if (isTRUE(HAS_IHW)) tagList(
          checkboxInput("deg_use_ihw", "Ponderar hipotesis con IHW en lugar de BH", FALSE),
          tags$small(class = "text-muted d-block mb-2",
                     paste("Solo DESeq2. IHW pondera cada gen segun su baseMean y gana",
                           "potencia sin perder el control de la FDR, pero necesita",
                           "muchos tests para poder formar bins: con pocos genes se",
                           "reduce a BH."))
        ) else NULL,
        actionButton("run_deg_btn",
                     tagList(icon("flask"), " Lanzar DEG"),
                     class = "btn-success btn-lg"),
        tags$div(style = "margin-top:8px;",
                 verbatimTextOutput("deg_status_text", placeholder = TRUE))
      ),

      # Card 5: variacion no deseada. La separacion entre las dos mitades de esta
      # tarjeta es el punto didactico: corregir-para-testear y
      # corregir-para-visualizar no son la misma operacion.
      card(
        card_header("5. Variacion no deseada"),
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
      ),

      # Card 6: filtros que solo recortan lo que se ve. Deliberadamente separados
      # de la card 4 para que no se confundan con los umbrales del test.
      card(
        card_header("6. Filtros de visualizacion"),
        div(class = "alert alert-secondary py-2 px-2 small mb-2",
            icon("eye"),
            tags$b(" Solo afectan a lo que se muestra."),
            " No cambian el modelo ni el control de la FDR. Para exigir un",
            " fold-change minimo con garantia estadistica usa el umbral del test",
            " en la tarjeta 4."),
        sliderInput("deg_log2fc_cutoff", "|log2FC| minimo (visual)",
                    min = 0, max = 5, value = 0, step = 0.1),
        sliderInput("deg_basemean_cutoff", "baseMean minimo (visual)",
                    min = 0, max = 1000, value = 0, step = 5)
      )
    ),

    # Zona de resultados (cargada al lanzar run_deg)
    tags$div(class = "mt-3", uiOutput("deg_results_ui"))
  )
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
        card_header(tags$div(
          class = "card-title-download",
          tags$span("Distribucion de p-valores"),
          div(style = "display:flex;gap:6px;align-items:center;",
              selectInput("deg_diag_pv_subset", NULL,
                          choices = c("Genes que pasan el filtrado independiente" = "tested",
                                      "Todos los genes con p-valor" = "all"),
                          selected = "tested", width = "320px"),
              downloadButton("download_deg_pvalue_hist", label = NULL,
                             icon = icon("download"),
                             class = "btn-sm btn-outline-secondary header-download",
                             title = "Descargar histograma"))
        )),
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
    navset_tab(
      id = "deg_result_tabs",
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
        "Diagnosticos",
        ui_tab_deg_diagnostics()
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
      nav_panel(
        "PCA",
        card(
          download_header("PCA (vst/rlog)", "download_deg_pca_plot"),
          plotly::plotlyOutput("deg_pca_plot", height = "480px")
        )
      ),
      nav_panel(
        "Heatmap top-N",
        card(
          card_header(tags$div(
            class = "card-title-download",
            tags$span("Heatmap top-N genes"),
            div(style = "display:flex;gap:6px;align-items:center;",
              numericInput("deg_heatmap_topn", NULL, value = 30, min = 5, max = 200, step = 5, width = "100px"),
              downloadButton("download_deg_heatmap", label = NULL, icon = icon("download"),
                             class = "btn-sm btn-outline-secondary header-download",
                             title = "Descargar heatmap como PNG")
            )
          )),
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
      nav_panel(
        "Enriquecimiento",
        card(
          card_header(tags$div(
            class = "card-title-download",
            tags$span("Enriquecimiento funcional"),
            div(style = "display:flex;gap:6px;align-items:flex-end;flex-wrap:wrap;",
                # ORA parte de la lista umbralizada y por tanto es ciego a
                # senales debiles pero coordinadas; GSEA usa el ranking completo.
                selectInput("deg_enrich_approach", NULL,
                            choices = c(c("ORA (lista significativa)" = "ora"),
                                        if (isTRUE(HAS_FGSEA))
                                          c("GSEA (ranking completo)" = "gsea")),
                            selected = "ora", width = "215px"),
                selectInput("deg_ontology", NULL,
                            choices = c("GO: Procesos biologicos" = "BP",
                                        "GO: Funcion molecular" = "MF",
                                        "GO: Componente celular" = "CC",
                                        "KEGG" = "KEGG"),
                            selected = "BP", width = "200px"),
                conditionalPanel(
                  "input.deg_enrich_approach === 'gsea'",
                  selectInput("deg_gsea_metric", NULL,
                              choices = c("Metrica: stat" = "stat",
                                          "Metrica: log2FC" = "log2FC",
                                          "Metrica: signo x -log10(p)" = "signed_p"),
                              selected = "stat", width = "220px")
                ),
                # keyType: los IDs de featureCounts sobre un GFF procariota son
                # locus tags, no simbolos. Fijarlo a SYMBOL hacia fallar el mapeo
                # en silencio, asi que ahora es explicito y seleccionable.
                conditionalPanel(
                  "input.deg_ontology !== 'KEGG'",
                  selectInput("deg_go_keytype", NULL,
                              choices = c("SYMBOL"), selected = "SYMBOL",
                              width = "150px")
                ),
                conditionalPanel(
                  "input.deg_ontology === 'KEGG'",
                  textInput("deg_kegg_organism", NULL, value = "eco",
                            placeholder = "eco, hsa, mmu...", width = "110px")
                ),
                conditionalPanel(
                  "input.deg_ontology === 'KEGG'",
                  selectInput("deg_kegg_keytype", NULL,
                              choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"),
                              selected = "kegg", width = "150px")
                ),
                actionButton("deg_run_enrich_btn",
                             tagList(icon("play"), " Calcular"),
                             class = "btn-sm"))
          )),
          # simplify() solo aplica al ORA sobre GO: colapsa terminos redundantes
          # por similitud semantica, y necesita el grafo de una sola ontologia.
          conditionalPanel(
            "input.deg_enrich_approach === 'ora' && input.deg_ontology !== 'KEGG'",
            checkboxInput("deg_go_simplify",
                          "Colapsar terminos GO redundantes (simplify, similitud de Wang > 0,7)",
                          FALSE)
          ),
          uiOutput("deg_enrich_mapping"),
          plotly::plotlyOutput("deg_enrich_dotplot", height = "440px"),
          tags$hr(),
          tags$div(
            class = "card-title-download",
            style = "margin-bottom:6px;",
            tags$span("Tabla de terminos enriquecidos"),
            downloadButton("download_enrich_table", label = NULL, icon = icon("download"),
                           class = "btn-sm btn-outline-secondary header-download",
                           title = "Descargar tabla de enriquecimiento")
          ),
          DTOutput("deg_enrich_table")
        )
      )
    )
  )
}
