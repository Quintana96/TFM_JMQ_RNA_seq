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
  # Contexto del ultimo GSEA (ranking y parametros). Lo necesita el running
  # score: la tabla de resultados no contiene el ranking, y recalcularlo al
  # vuelo daria una curva que no corresponde al NES mostrado si el usuario ha
  # tocado la metrica entretanto.
  gsea_ctx_rv <- reactiveVal(NULL)
  compare_rv <- reactiveVal(NULL)

  # El OrgDb disponible hoy es solo el de E. coli K12; se centraliza aqui para
  # que el selector de keyType y el enriquecimiento no se desincronicen.
  # OrgDb elegido en la interfaz. Estaba cableado a org.EcK12.eg.db, de modo que
  # con datos de otro organismo el enriquecimiento GO no podia funcionar aunque
  # su paquete estuviera instalado: el mapeo salia del 0 % y la app avisaba,
  # correctamente, de algo que el usuario no tenia forma de arreglar.
  deg_orgdb <- reactive({
    sel <- input$deg_orgdb %||% ""
    if (nzchar(sel) && sel %in% ORGDBS_DISPONIBLES) return(sel)
    if (length(ORGDBS_DISPONIBLES)) return(ORGDBS_DISPONIBLES[1])
    NULL
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

  # ── Gene sets propios (GMT) ────────────────────────────────────────────────
  gmt_rv <- reactive({
    f <- input$deg_gmt_file
    if (is.null(f) || !nrow(f)) return(NULL)
    read_gene_sets_gmt(f$datapath[1])
  })

  output$deg_gmt_summary <- renderUI({
    gs <- gmt_rv()
    if (is.null(gs)) {
      return(tags$small(class = "text-muted", "Ningun fichero cargado."))
    }
    if (!is.null(gs$error)) {
      return(div(class = "alert alert-warning py-1 px-2 small mb-0",
                 icon("triangle-exclamation"), " ", gs$error))
    }
    div(class = "small text-muted", tags$b("Conjuntos cargados: "),
        gmt_summary_text(gs))
  })

  # ── Parametros de la interfaz ──────────────────────────────────────────────
  # Un unico punto de lectura para que el boton de calcular y el de comparar no
  # puedan divergir en los parametros que usan.
  num_or <- function(x, default) {
    v <- suppressWarnings(as.numeric(x %||% NA))
    if (length(v) != 1L || is.na(v)) default else v
  }

  enrich_inputs <- function() {
    ont <- input$deg_ontology %||% "BP"
    org_code <- trimws(input$deg_kegg_organism %||% "eco")
    if (!nzchar(org_code)) org_code <- "eco"
    gs <- gmt_rv()
    list(
      ont         = ont,
      approach    = input$deg_enrich_approach %||% "ora",
      org_code    = org_code,
      keytype     = if (identical(ont, "KEGG")) input$deg_kegg_keytype %||% "kegg"
                    else input$deg_go_keytype %||% "SYMBOL",
      # setReadable necesita OrgDb: no aplica a KEGG ni a conjuntos propios.
      readable    = isTRUE(input$deg_enrich_readable) && !ont %in% c("KEGG", "GMT"),
      simplify    = isTRUE(input$deg_go_simplify),
      directional = isTRUE(input$deg_ora_directional),
      metric      = input$deg_gsea_metric %||% "stat",
      min_size    = max(2, round(num_or(input$deg_gsea_min_size, 10))),
      max_size    = max(3, round(num_or(input$deg_gsea_max_size, 500))),
      pcut        = min(1, max(1e-4, num_or(input$deg_gsea_pcutoff, 0.05))),
      # eps = 0 pide a fgsea el p-valor exacto del multilevel en lugar de
      # truncarlo en 1e-10, donde todos los conjuntos muy significativos empatan.
      eps         = if (isTRUE(input$deg_gsea_eps_exact)) 0 else 1e-10,
      universe    = ctx$deg_universe(),
      term2gene   = if (is.null(gs)) NULL else gs$term2gene,
      gmt         = gs
    )
  }

  # ── Ejecucion ──────────────────────────────────────────────────────────────
  #' Corre el ORA de la ontologia seleccionada sobre una lista de genes.
  ora_runner <- function(p) {
    function(genes) {
      if (identical(p$ont, "KEGG")) {
        run_enrichment_kegg(genes, universe = p$universe, organism = p$org_code,
                            keyType = p$keytype)
      } else if (identical(p$ont, "GMT")) {
        run_enrichment_gmt(genes, universe = p$universe, term2gene = p$term2gene,
                           minGSSize = p$min_size, maxGSSize = p$max_size)
      } else {
        run_enrichment_go(genes, universe = p$universe, OrgDb = deg_orgdb(),
                          ont = p$ont, keyType = p$keytype,
                          simplify_terms = p$simplify, readable = p$readable)
      }
    }
  }

  #' Ejecuta un enfoque (ORA o GSEA) con los parametros dados.
  #' Devuelve el list(table, error, mapping) de utils_enrich mas, en GSEA, el
  #' ranking usado.
  enrich_compute <- function(p, approach) {
    if (identical(p$ont, "GMT") && is.null(p$term2gene)) {
      return(list(table = NULL, mapping = NULL,
                  error = "Carga un fichero GMT con los conjuntos de genes."))
    }
    if (identical(approach, "gsea")) {
      # GSEA parte del ranking completo, sin umbralizar: por eso usa la tabla
      # entera y no la lista filtrada.
      rk <- deg_ranking_metric(state$deg_rv$results, p$metric)
      if (is.null(rk$ranked)) {
        return(list(table = NULL, mapping = NULL, ranking = rk,
                    error = paste0("No se pudo construir el ranking: ",
                                   rk$error %||% "—")))
      }
      res <- run_gsea(rk$ranked, ont = p$ont, OrgDb = deg_orgdb(),
                      organism = p$org_code, keyType = p$keytype, exponent = 0,
                      pvalueCutoff = p$pcut, minGSSize = p$min_size,
                      maxGSSize = p$max_size, eps = p$eps,
                      term2gene = p$term2gene, readable = p$readable)
      res$ranking <- rk
      return(res)
    }
    # La lista del ORA sale del FDR con el que se AJUSTO el modelo, no de los
    # filtros de la tarjeta 6: esos son de visualizacion y no recortan el
    # universo, asi que dejarlos entrar aqui cambiaba el enriquecimiento al
    # mover un deslizador declarado cosmetico.
    df <- ctx$deg_significant()
    if (is.null(df) || !nrow(df)) {
      return(list(table = NULL, mapping = NULL,
                  error = paste0("No hay genes significativos a FDR <= ",
                                 state$deg_rv$fdr %||% 0.05, ".")))
    }
    runner <- ora_runner(p)
    res <- if (isTRUE(p$directional)) run_ora_directional(df, runner)
           else runner(unique(as.character(df$gene)))
    res$n_lista <- nrow(df)
    res
  }

  observeEvent(input$deg_run_enrich_btn, {
    req(state$deg_rv$results)
    # El enriquecimiento tarda entre segundos y un minuto largo (mas con
    # simplify(), que calcula similitud semantica entre todos los terminos, con
    # el ORA direccional, que lo corre tres veces, y con eps = 0, que sustituye
    # el p-valor truncado por el exacto). Sin feedback la app parece colgada y el
    # usuario vuelve a pulsar, encolando ejecuciones.
    shinyjs::disable("deg_run_enrich_btn")
    on.exit(shinyjs::enable("deg_run_enrich_btn"), add = TRUE)

    p <- enrich_inputs()
    approach <- p$approach

    if (identical(approach, "gsea")) {
      rk_check <- deg_ranking_metric(state$deg_rv$results, p$metric)
      if (!is.null(rk_check$ranked) && !is.na(rk_check$ties_frac) &&
          rk_check$ties_frac > 0.01) {
        showNotification(
          paste0("El ", round(100 * rk_check$ties_frac, 1), " % de los genes comparte ",
                 "valor en el ranking. GSEA no resuelve los empates, asi que su ",
                 "orden interno es arbitrario; considera usar la metrica 'stat'."),
          type = "warning", duration = 16
        )
      }
    }

    res <- withProgress(
      message = if (identical(approach, "gsea")) "Calculando GSEA..."
                else "Calculando enriquecimiento...",
      value = 0.3, {
        setProgress(value = 0.6,
                    detail = if (identical(approach, "gsea")) "permutaciones"
                             else "test hipergeometrico")
        enrich_compute(p, approach)
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
      ontologia  = if (identical(p$ont, "GMT")) "Gene sets propios (GMT)" else p$ont,
      organismo_kegg = if (identical(p$ont, "KEGG")) p$org_code else NA_character_,
      orgdb      = if (identical(p$ont, "GMT")) "no aplica (GMT)" else deg_orgdb() %||% "—",
      keytype    = if (identical(p$ont, "GMT")) "IDs del fichero GMT" else p$keytype,
      metrica    = if (identical(approach, "gsea")) p$metric else NA_character_,
      exponent   = if (identical(approach, "gsea")) 0 else NA_real_,
      simplify   = if (identical(approach, "gsea")) NA else p$simplify,
      # Los parametros de GSEA cambian que conjuntos se testean y cuales se
      # devuelven, asi que sin ellos el resultado no es reproducible.
      gsea_min_size = if (identical(approach, "gsea")) p$min_size else NA_integer_,
      gsea_max_size = if (identical(approach, "gsea")) p$max_size else NA_integer_,
      gsea_pcutoff  = if (identical(approach, "gsea")) p$pcut else NA_real_,
      gsea_eps      = if (identical(approach, "gsea")) p$eps else NA_real_,
      direccional = if (identical(approach, "gsea")) NA else p$directional,
      simbolos   = p$readable,
      gmt        = if (identical(p$ont, "GMT") && !is.null(p$gmt))
                     paste0(p$gmt$n_sets, " conjuntos (", min(p$gmt$sizes), "-",
                            max(p$gmt$sizes), " genes)")
                   else NA_character_,
      n_lista    = if (identical(approach, "gsea")) NA_integer_
                   else res$n_lista %||% nrow(ctx$deg_significant() %||% data.frame()),
      n_universo = length(p$universe %||% character(0)),
      mapeo      = res$mapping,
      n_terminos = if (is.null(res$table)) 0L else nrow(res$table),
      error      = res$error %||% NA_character_,
      # Fechas de las fuentes del OrgDb: los resultados de enriquecimiento
      # cambian con la version de la anotacion (Wadi et al., Nat Methods 2016),
      # asi que sin esto el resultado no es reproducible.
      anotacion  = orgdb_source_info(deg_orgdb())
    )

    # El contexto de GSEA se guarda aunque no haya terminos: es lo que permite
    # que el running score sepa sobre que ranking se calculo todo.
    if (identical(approach, "gsea")) {
      gsea_ctx_rv(list(ranked = res$ranking$ranked, ont = p$ont,
                       keytype = p$keytype, orgdb = deg_orgdb(),
                       organism = p$org_code, term2gene = p$term2gene,
                       exponent = 0))
    } else {
      gsea_ctx_rv(NULL)
    }

    if (is.null(res$table)) {
      # Con GSEA, "sin conjuntos enriquecidos" puede ser simplemente que ninguno
      # baja de 0,05: se sugiere el corte a 1 para poder distinguirlo de un fallo.
      extra <- if (identical(approach, "gsea") && p$pcut < 1)
        " Pon el corte de p ajustado a 1 para ver todos los conjuntos testeados." else ""
      showNotification(paste0("Enriquecimiento sin resultados: ", res$error %||% "—", extra),
                       type = "warning", duration = 10)
      enrich_rv(NULL); return()
    }
    if (isTRUE(p$directional) && length(res$errores %||% character(0))) {
      showNotification(
        paste0("Direcciones sin terminos: ",
               paste0(names(res$errores), " (", res$errores, ")", collapse = "; ")),
        type = "message", duration = 10)
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
        y = ~stats::reorder(plot_label, NES),
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

    # Con ORA direccional el color codifica la direccion: es la lectura que se
    # busca al separar las listas, y sin ella las dos mitades se confunden.
    dir_col <- if ("Direccion" %in% names(top)) {
      unname(c("Conjunto" = "#A8DADC", "Al alza" = "#7BBF9A",
               "A la baja" = "#F4A6A6")[top$Direccion])
    } else NULL

    p <- plotly::plot_ly(
      top,
      x = ~ -log10(p.adjust),
      y = ~stats::reorder(plot_label, -log10(p.adjust)),
      type = "scatter", mode = "markers",
      size = ~Count,
      marker = if (is.null(dir_col))
        list(color = ~ -log10(p.adjust),
             colorscale = list(c(0, "#A8DADC"), c(1, "#7BBF9A")),
             sizemode = "area")
      else list(color = dir_col, sizemode = "area"),
      text = ~paste0("Termino: ", Description,
                     if ("Direccion" %in% names(top)) paste0("<br>Direccion: ", Direccion) else "",
                     "<br>Count: ", Count,
                     "<br>p.adjust: ", signif(p.adjust, 3)),
      hoverinfo = "text"
    )
    p |>
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
    for (nm in intersect(c("enrichmentScore", "NES"), names(df_r)))
      df_r[[nm]] <- round(df_r[[nm]], 3)
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

  # ── Running score / leading edge ───────────────────────────────────────────
  # Solo tiene sentido sobre una tabla de GSEA (la que trae NES).
  observe({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df) || !"NES" %in% names(df)) {
      updateSelectInput(session, "deg_gsea_term", choices = character(0))
      return()
    }
    ord <- order(df$p.adjust, na.last = TRUE)
    ids <- as.character(df$ID[ord])
    etiquetas <- paste0(df$Description[ord], "  (NES ", round(df$NES[ord], 2),
                        ", p.adj ", signif(df$p.adjust[ord], 2), ")")
    updateSelectInput(session, "deg_gsea_term",
                      choices = stats::setNames(ids, etiquetas),
                      selected = isolate(input$deg_gsea_term) %||% ids[1])
  })

  gsea_runscore_rv <- reactive({
    gctx <- gsea_ctx_rv()
    term <- input$deg_gsea_term
    if (is.null(gctx) || is.null(term) || !nzchar(term %||% "")) return(NULL)
    tg <- gsea_term_genes(term, ont = gctx$ont, OrgDb = gctx$orgdb,
                          keyType = gctx$keytype, organism = gctx$organism,
                          term2gene = gctx$term2gene)
    if (!length(tg$genes)) {
      return(list(error = tg$error %||% "No se han podido recuperar los genes del conjunto."))
    }
    # gseaParam = exponent: la curva debe ser la del estadistico con el que se
    # calculo el NES, no la ponderada por defecto de fgsea.
    gsea_running_score(gctx$ranked, tg$genes, gseaParam = gctx$exponent %||% 0)
  })

  output$deg_gsea_runscore <- plotly::renderPlotly({
    if (is.null(gsea_ctx_rv()))
      return(plotly_message("Corre un GSEA para ver el running score."))
    rs <- gsea_runscore_rv()
    if (is.null(rs)) return(plotly_message("Selecciona un conjunto."))
    if (!is.null(rs$error)) return(plotly_message(rs$error))
    cur <- rs$curve
    base <- min(cur$ES, 0) - 0.08 * max(abs(range(cur$ES)), 0.1)
    plotly::plot_ly() |>
      plotly::add_lines(x = cur$rank, y = cur$ES, line = list(color = "#3F7D5C"),
                        name = "Running score", hoverinfo = "x+y") |>
      # Cada marca es un gen del conjunto en su posicion del ranking: es lo que
      # deja ver si el conjunto se concentra en un extremo o esta repartido.
      plotly::add_markers(x = rs$ticks$rank, y = rep(base, nrow(rs$ticks)),
                          marker = list(symbol = "line-ns-open", color = "#60756A",
                                        size = 8),
                          name = "Genes del conjunto", hoverinfo = "x") |>
      plotly::layout(
        xaxis = list(title = "Posicion en el ranking"),
        yaxis = list(title = "Enrichment score acumulado"),
        showlegend = FALSE,
        shapes = list(list(type = "line", x0 = 0, x1 = rs$n_ranked, y0 = 0, y1 = 0,
                           line = list(dash = "dot", color = "#B0BDB5")))
      )
  })

  output$deg_gsea_leading <- renderUI({
    df <- enrich_rv()
    term <- input$deg_gsea_term
    if (is.null(df) || !"core_enrichment" %in% names(df) ||
        is.null(term) || !nzchar(term %||% "")) return(NULL)
    fila <- df[!is.na(df$ID) & df$ID == term, , drop = FALSE]
    if (!nrow(fila)) return(NULL)
    genes <- leading_edge_genes(fila$core_enrichment[1])
    rs <- gsea_runscore_rv()
    tagList(
      div(class = "small text-muted mt-2",
          tags$b("Conjunto: "), fila$Description[1],
          " — NES ", round(fila$NES[1], 3),
          ", p ajustado ", signif(fila$p.adjust[1], 3),
          ", ", fila$setSize[1], " genes",
          if (!is.null(rs) && is.null(rs$error))
            paste0(" (", rs$n_hits, " presentes en el ranking)") else ""),
      div(class = "small mt-1",
          tags$b(paste0("Leading edge (", length(genes), " genes): ")),
          tags$span(class = "text-muted",
                    "el subconjunto que va desde el extremo del ranking hasta el pico ",
                    "de la curva y que sostiene el enriquecimiento.")),
      div(class = "small mt-1 p-2 border rounded",
          style = "max-height:150px;overflow-y:auto;font-family:monospace;",
          paste(genes, collapse = ", "))
    )
  })

  # ── ORA frente a GSEA ──────────────────────────────────────────────────────
  observeEvent(input$deg_run_enrich_compare_btn, {
    req(state$deg_rv$results)
    if (!isTRUE(HAS_FGSEA)) {
      showNotification("fgsea no esta instalado: no se puede correr el GSEA.",
                       type = "error"); return()
    }
    shinyjs::disable("deg_run_enrich_compare_btn")
    on.exit(shinyjs::enable("deg_run_enrich_compare_btn"), add = TRUE)

    p <- enrich_inputs()
    # La comparacion exige que las dos ramas vean lo mismo: misma ontologia,
    # mismo universo y misma lista. El ORA direccional partiria la lista en tres
    # y dejaria de ser comparable con un unico GSEA.
    p$directional <- FALSE

    res <- withProgress(message = "Comparando ORA y GSEA...", value = 0.2, {
      ora <- enrich_compute(p, "ora")
      setProgress(value = 0.6, detail = "GSEA (permutaciones)")
      gsea <- enrich_compute(p, "gsea")
      list(ora = ora, gsea = gsea)
    })

    if (is.null(res$ora$table) && is.null(res$gsea$table)) {
      compare_rv(NULL)
      showNotification(paste0("Ninguno de los dos enfoques devolvio terminos. ORA: ",
                              res$ora$error %||% "—", " | GSEA: ",
                              res$gsea$error %||% "—"),
                       type = "warning", duration = 12)
      return()
    }
    cmp <- compare_ora_gsea(res$ora$table, res$gsea$table, padj_cutoff = 0.05)
    cmp$ont <- p$ont
    cmp$error_ora <- res$ora$error
    cmp$error_gsea <- res$gsea$error
    compare_rv(cmp)
  })

  output$deg_enrich_compare_summary <- renderUI({
    cmp <- compare_rv()
    if (is.null(cmp)) {
      return(div(class = "small text-muted",
                 "Pulsa 'Comparar ORA y GSEA' para correr los dos enfoques."))
    }
    jac <- if (is.na(cmp$jaccard)) "—" else paste0(round(100 * cmp$jaccard, 1), " %")
    tagList(
      div(class = "alert alert-secondary py-2 px-2 small mb-2",
          tags$b("Terminos significativos (p ajustado <= 0,05): "),
          "ORA ", fmt_int(cmp$n_ora), " | GSEA ", fmt_int(cmp$n_gsea),
          " | comunes ", fmt_int(cmp$n_comun),
          " | indice de Jaccard ", jac, ".",
          tags$br(),
          if (cmp$n_gsea > cmp$n_comun)
            paste0("GSEA ve ", cmp$n_gsea - cmp$n_comun, " terminos que el ORA no ve: ",
                   "conjuntos desplazados en el ranking cuyos genes no llegan al corte ",
                   "de la lista.")
          else "GSEA no anade ningun termino sobre los del ORA."),
      if (!is.na(cmp$error_ora %||% NA))
        div(class = "small text-muted", tags$b("ORA: "), cmp$error_ora) else NULL,
      if (!is.na(cmp$error_gsea %||% NA))
        div(class = "small text-muted", tags$b("GSEA: "), cmp$error_gsea) else NULL
    )
  })

  output$deg_enrich_compare_plot <- plotly::renderPlotly({
    cmp <- compare_rv()
    if (is.null(cmp)) return(plotly_message("Sin comparacion calculada."))
    d <- data.frame(
      grupo = factor(c("Solo ORA", "Comunes", "Solo GSEA"),
                     levels = c("Solo ORA", "Comunes", "Solo GSEA")),
      n = c(cmp$n_ora - cmp$n_comun, cmp$n_comun, cmp$n_gsea - cmp$n_comun)
    )
    plotly::plot_ly(d, x = ~grupo, y = ~n, type = "bar",
                    marker = list(color = c("#A8DADC", "#7BBF9A", "#F4C68A")),
                    text = ~n, textposition = "outside", hoverinfo = "x+y") |>
      plotly::layout(xaxis = list(title = ""),
                     yaxis = list(title = "Terminos significativos"),
                     showlegend = FALSE)
  })

  compare_table_rv <- reactive({
    cmp <- compare_rv()
    if (is.null(cmp)) return(NULL)
    tb <- compare_ora_gsea_table(cmp)
    if (is.null(tb)) return(NULL)
    for (nm in c("padj_ORA", "padj_GSEA")) tb[[nm]] <- signif(tb[[nm]], 3)
    tb$NES <- round(tb$NES, 3)
    tb
  })

  output$deg_enrich_compare_table <- renderDT({
    tb <- compare_table_rv()
    if (is.null(tb)) return(dt_table(message_df("Sin comparacion calculada.")))
    dt_table(tb, page_length = 15, filter = "top")
  })

  output$download_enrich_compare <- csv_download(
    "comparacion_ora_gsea",
    function() {
      tb <- compare_table_rv()
      if (is.null(tb)) return(message_df("Sin comparacion calculada."))
      tb
    }
  )

  invisible(NULL)
}
