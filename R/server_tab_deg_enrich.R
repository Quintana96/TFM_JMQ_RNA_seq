#' server_tab_deg_enrich.R
#' Enriquecimiento funcional y GSEA de la pestana 4 (items 13 y 14, B3).
#'
#' Parte del modulo de la pestana 4, separado de server_tab_deg.R por tamano: el
#' fichero original llego a 1.608 lineas en una unica funcion. NO se usa
#' moduleServer(): igual que el resto de la aplicacion, se conservan los IDs
#' Shiny originales y el estado se pasa explicitamente.
#'
#' `ctx` es el contexto compartido del modulo (un environment, como `state`, para
#' que las asignaciones sean visibles entre las partes). Lo crea server_tab_deg()
#' y contiene los reactivos que varias partes necesitan.

server_tab_deg_enrich <- function(input, output, session, state, ctx) {
  # ── Enriquecimiento ────────────────────────────────────────────────────────
  enrich_rv <- reactiveVal(NULL)
  enrich_mapping_rv <- reactiveVal(NULL)

  # El OrgDb disponible hoy es solo el de E. coli K12; se centraliza aqui para
  # que el selector de keyType y el enriquecimiento no se desincronicen.
  deg_orgdb <- reactive({
    if (isTRUE(HAS_ORGECDB)) "org.EcK12.eg.db" else NULL
  })

  # keyType seleccionable, poblado con los keytypes reales del OrgDb. Antes
  # estaba fijado a "SYMBOL" en el codigo, lo que con locus tags de E. coli
  # (b0001) daba un mapeo casi nulo sin ningun aviso.
  observe({
    kt <- orgdb_keytypes(deg_orgdb())
    if (!length(kt)) kt <- c("SYMBOL", "ENTREZID", "ENSEMBL", "ALIAS", "REFSEQ")
    preferred <- if ("SYMBOL" %in% kt) "SYMBOL" else kt[1]
    updateSelectInput(session, "deg_go_keytype", choices = kt,
                      selected = isolate(input$deg_go_keytype) %||% preferred)
  })

  observeEvent(input$deg_run_enrich_btn, {
    req(state$deg_rv$results)
    # El enriquecimiento tarda entre segundos y un minuto largo (mas con
    # simplify(), que calcula similitud semantica entre todos los terminos). Sin
    # feedback la app parece colgada y el usuario vuelve a pulsar, encolando
    # ejecuciones. El resto de acciones largas (ajuste DEG, bootstrap) ya
    # bloquean su boton y muestran progreso; esta era la excepcion.
    shinyjs::disable("deg_run_enrich_btn")
    on.exit(shinyjs::enable("deg_run_enrich_btn"), add = TRUE)

    ont <- input$deg_ontology %||% "BP"
    approach <- input$deg_enrich_approach %||% "ora"
    org_code <- trimws(input$deg_kegg_organism %||% "eco")
    if (!nzchar(org_code)) org_code <- "eco"
    # Fondo = los genes efectivamente EVALUABLES (los que tienen padj), no el
    # genoma completo ni la tabla entera. Un gen descartado por el filtrado
    # independiente nunca podria haber entrado en la lista de significativos, asi
    # que contarlo como fondo infla el enriquecimiento. Aplica a GO y a KEGG.
    universe <- ctx$deg_universe()

    res <- withProgress(
      message = if (identical(approach, "gsea")) "Calculando GSEA..."
                else "Calculando enriquecimiento...",
      value = 0.3, {
    if (identical(approach, "gsea")) {
      # GSEA parte del ranking completo, sin umbralizar: por eso usa la tabla
      # entera y no la lista filtrada.
      metric <- input$deg_gsea_metric %||% "stat"
      rk <- deg_ranking_metric(state$deg_rv$results, metric)
      if (is.null(rk$ranked)) {
        showNotification(paste0("No se pudo construir el ranking: ",
                                rk$error %||% "—"), type = "error", duration = 10)
        enrich_rv(NULL); return()
      }
      if (!is.na(rk$ties_frac) && rk$ties_frac > 0.01) {
        showNotification(
          paste0("El ", round(100 * rk$ties_frac, 1), " % de los genes comparte ",
                 "valor en el ranking. GSEA no resuelve los empates, asi que su ",
                 "orden interno es arbitrario; considera usar la metrica 'stat'."),
          type = "warning", duration = 16
        )
      }
      setProgress(value = 0.6, detail = "permutaciones")
      run_gsea(rk$ranked, ont = ont, OrgDb = deg_orgdb(),
               organism = org_code,
               keyType = if (identical(ont, "KEGG"))
                 input$deg_kegg_keytype %||% "kegg"
               else input$deg_go_keytype %||% "SYMBOL",
               exponent = 0)
    } else {
      # La lista del ORA sale del FDR con el que se AJUSTO el modelo, no de los
      # filtros de la tarjeta 5: esos son de visualizacion y no recortan el
      # universo, asi que dejarlos entrar aqui cambiaba el enriquecimiento al
      # mover un deslizador declarado cosmetico.
      df <- ctx$deg_significant()
      if (is.null(df) || !nrow(df)) {
        showNotification(paste0("No hay genes significativos a FDR <= ",
                                state$deg_rv$fdr %||% 0.05, "."),
                         type = "warning"); return()
      }
      genes <- df$gene
      setProgress(value = 0.6, detail = paste0(length(genes), " genes"))
      if (identical(ont, "KEGG")) {
        run_enrichment_kegg(genes, universe = universe, organism = org_code,
                            keyType = input$deg_kegg_keytype %||% "kegg")
      } else {
        run_enrichment_go(genes, universe = universe,
                          OrgDb = deg_orgdb(), ont = ont,
                          keyType = input$deg_go_keytype %||% "SYMBOL",
                          simplify_terms = isTRUE(input$deg_go_simplify))
      }
    }
    })

    enrich_mapping_rv(res$mapping)
    # Un mapeo bajo hace el resultado no interpretable, asi que se avisa aunque
    # el enriquecimiento haya devuelto terminos.
    mp <- res$mapping
    if (!is.null(mp) && !is.na(mp$rate %||% NA) && mp$rate < 0.5) {
      showNotification(
        paste0("Atencion: solo se ha mapeado el ", round(100 * mp$rate, 1),
               " % de los genes (", mp$n_mapped, "/", mp$n_input,
               "). Revisa el keyType: los IDs de featureCounts suelen ser locus ",
               "tags, no simbolos."),
        type = "warning", duration = 16
      )
    }

    # Se registra el enriquecimiento AUNQUE no haya dado terminos: un resultado
    # negativo tambien forma parte del analisis, y sin sus parametros (ontologia,
    # keyType, universo, tasa de mapeo) no es interpretable ni reproducible.
    state$deg_rv$enrich <- list(
      enfoque    = if (identical(approach, "gsea")) "GSEA" else "ORA (sobre-representacion)",
      ontologia  = ont,
      organismo_kegg = if (identical(ont, "KEGG")) org_code else NA_character_,
      orgdb      = deg_orgdb() %||% "—",
      keytype    = if (identical(ont, "KEGG")) input$deg_kegg_keytype %||% "kegg"
                   else input$deg_go_keytype %||% "SYMBOL",
      metrica    = if (identical(approach, "gsea")) input$deg_gsea_metric %||% "stat" else NA_character_,
      exponent   = if (identical(approach, "gsea")) 0 else NA_real_,
      simplify   = if (identical(approach, "gsea")) NA else isTRUE(input$deg_go_simplify),
      n_lista    = if (identical(approach, "gsea")) NA_integer_
                   else nrow(ctx$deg_significant() %||% data.frame()),
      n_universo = length(universe %||% character(0)),
      mapeo      = res$mapping,
      n_terminos = if (is.null(res$table)) 0L else nrow(res$table),
      error      = res$error %||% NA_character_,
      # Fechas de las fuentes del OrgDb: los resultados de enriquecimiento
      # cambian con la version de la anotacion (Wadi et al., Nat Methods 2016),
      # asi que sin esto el resultado no es reproducible.
      anotacion  = orgdb_source_info(deg_orgdb())
    )

    if (is.null(res$table)) {
      showNotification(paste0("Enriquecimiento sin resultados: ", res$error %||% "—"),
                       type = "warning", duration = 8)
      enrich_rv(NULL); return()
    }
    enrich_rv(res$table)
  })

  output$deg_enrich_mapping <- renderUI({
    mp <- enrich_mapping_rv()
    txt <- mapping_rate_text(mp)
    if (is.null(txt)) return(NULL)
    low <- !is.null(mp) && !is.na(mp$rate %||% NA) && mp$rate < 0.5
    div(class = paste("small py-1 px-2 mb-1 rounded",
                      if (low) "alert alert-warning mb-2" else "text-muted"),
        if (low) tagList(icon("triangle-exclamation"), " ") else NULL,
        tags$b("Tasa de mapeo: "), txt)
  })

  #' El dotplot sirve a dos tipos de tabla: el ORA trae GeneRatio/Count y se
  #' ordena por significacion; GSEA trae NES/setSize, y ahi lo informativo es el
  #' signo del NES (si el conjunto sube o baja), no solo el p-valor.
  make_enrich_dotplot <- function(df) {
    if (is.null(df) || !nrow(df))
      return(plotly_message("Pulsa 'Calcular' para correr el enriquecimiento."))
    top <- enrichment_dotplot_data(df, top_n = 15)
    if (is.null(top) || !nrow(top))
      return(plotly_message("Sin terminos enriquecidos."))

    if ("NES" %in% names(top)) {
      top$size_val <- top$setSize %||% rep(10, nrow(top))
      return(plotly::plot_ly(
        top,
        x = ~NES,
        y = ~stats::reorder(Description, NES),
        type = "scatter", mode = "markers",
        size = ~size_val,
        marker = list(color = ~ifelse(NES > 0, "#7BBF9A", "#F4A6A6"),
                      sizemode = "area"),
        text = ~paste0("Conjunto: ", Description,
                       "<br>NES: ", round(NES, 3),
                       "<br>genes en el conjunto: ", size_val,
                       "<br>p.adjust: ", signif(p.adjust, 3)),
        hoverinfo = "text"
      ) |>
        plotly::layout(
          xaxis = list(title = "NES (score de enriquecimiento normalizado)"),
          yaxis = list(title = "", automargin = TRUE),
          margin = list(l = 220),
          shapes = list(list(type = "line", x0 = 0, x1 = 0, yref = "paper",
                             y0 = 0, y1 = 1,
                             line = list(dash = "dot", color = "#60756A")))
        ))
    }

    plotly::plot_ly(
      top,
      x = ~ -log10(p.adjust),
      y = ~stats::reorder(Description, -log10(p.adjust)),
      type = "scatter", mode = "markers",
      size = ~Count,
      marker = list(color = ~ -log10(p.adjust),
                    colorscale = list(c(0, "#A8DADC"), c(1, "#7BBF9A")),
                    sizemode = "area"),
      text = ~paste0("Termino: ", Description,
                     "<br>Count: ", Count,
                     "<br>p.adjust: ", signif(p.adjust, 3)),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        xaxis = list(title = "-log10(p.adjust)"),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 220)
      )
  }

  output$deg_enrich_dotplot <- plotly::renderPlotly({
    make_enrich_dotplot(enrich_rv())
  })

  output$deg_enrich_table <- renderDT({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df)) return(dt_table(message_df("Sin enriquecimiento calculado.")))
    df_r <- df
    for (nm in intersect(c("pvalue", "p.adjust", "qvalue"), names(df_r)))
      df_r[[nm]] <- signif(df_r[[nm]], 3)
    dt_table(df_r, page_length = 15, filter = "top")
  })

  output$download_enrich_table <- csv_download(
    "enriquecimiento",
    function() {
      df <- enrich_rv()
      if (is.null(df)) return(message_df("Sin enriquecimiento calculado."))
      df
    }
  )

  invisible(NULL)
}
