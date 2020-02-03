
#' @import shiny
#' @import shinyjs
#' @import shinyWidgets
#' @importFrom plotly plot_ly renderPlotly layout add_lines
#' @importFrom DT datatable renderDT
#' @importFrom heatmaply heatmaply
#' @importFrom reshape2 melt
#' @importFrom vegan vegdist
#' @importFrom dendextend seriate_dendrogram
server <- function(input, output, session){

  if (dim(pg$organisms)[2] < 2){
    shinyjs::disable(selector = "[type=radio][value=metadata]")
    shinyjs::runjs("$('[type=radio][value=metadata]').parent().parent().addClass('disabled').css('opacity', 0.4)")
    shinyjs::disable(selector = "[type=check][value=check_color_meta]")
    shinyjs::runjs("$('[type=check][value=check_color_meta]').parent().parent().addClass('disabled').css('opacity', 0.4)")
  }


  updateCoreLevel <- eventReactive(input$core_level, {
    pg$core_level <- input$core_level
    pg
  })

  accs_level <- reactiveVal(value = pg$core_level - 1L)

  observeEvent(input$core_level, {
    accs_level(input$core_level-1L)
  })

  orgs <- reactiveVal(as.character(pg$.__enclos_env__$private$.organisms$org))

  observeEvent(input$organisms, {
    orgs(input$organisms)
  })

  observeEvent(input$meta_dat, {
    updateOrganisms()
    org <- pg$.__enclos_env__$private$.organisms
    meta_column <- input$meta_column
    if (!is.null(meta_column)){
      wh <- which(org[[meta_column]] %in% input$meta_dat)
      sele <- as.character(org$org)[wh]
      orgs(sele)
    }
  })


  updateOrganisms <- reactive({
    ava_orgs <- as.character(pg$organisms[["org"]])
    inp_orgs <- orgs()

    toDrop <- ava_orgs[! ava_orgs %in% inp_orgs]
    toReco <- inp_orgs[! inp_orgs %in% ava_orgs]

    pg$drop(toDrop)
    pg$recover(toReco)
  })


  output$core_evol <- renderPlotly({
    updateOrganisms()
    pm <- pg$pan_matrix
    sq <- seq(1, 0.85, -0.01)
    cs <- colSums(pm)
    ev <- sapply(sq, function(x){
      length(which( cs >= (dim(pm)[1]*x) ))
    })
    df <- data.frame(core_level = sq*100, core_number = ev)
    plot_ly(df, x = ~core_level, y = ~core_number,
            type = "scatter",
            mode = "markers", height = 230) %>%
      layout(xaxis = list(title = "Core Level"),
             yaxis = list(title = "Core Number"))
  })

  output$summary_counts <- DT::renderDT({
    updateOrganisms()
    v1 <- c("Number of organisms",
            "Number of clusters",
            "Number of genes")
    v2 <- c(dim(pg$organisms)[1],
            dim(pg$clusters)[1],
            sum(pg$pan_matrix))
    df <- data.frame(v1, v2)
    datatable(df,
              selection = "none", colnames = "", rownames=F,
              options = list(
                # colnames = "",
                # rownames = FALSE,
                paging = FALSE,
                dom ='t',
                bSort=FALSE
              ))
  })


  output$select_organisms <- renderUI({

    if (input$input_type == "organism"){
      all_orgs <- as.character(pg$.__enclos_env__$private$.organisms$org)
      fluidRow(
        column(
          width = 6,
          pickerInput(inputId = "organisms",
                      label = "Organisms",
                      choices = all_orgs,
                      selected = all_orgs,#selectize = T,
                      multiple = TRUE, options = list(`actions-box` = TRUE),
                      width = "70%")
        )
      )

    }else{

      cols <- colnames(pg$organisms[, -1])
      meta_column <- input$meta_column
      if (is.null(meta_column)){
        col <- colnames(pg$organisms)[2]
      }else{
        col <- meta_column
      }
      chs <- as.character(pg$organisms[[col]])
      fluidRow(
        column(
          width = 3,
          selectInput(inputId = "meta_column",
                      label = "Select metadata:",
                      selected = cols[1],
                      choices = cols,
                      multiple = FALSE)
        ),
        column(
          width = 3,
          pickerInput(inputId = "meta_dat",
                      label = "Select organism group:",
                      selected = chs,
                      choices = chs,
                      multiple = TRUE, options = list(`actions-box` = TRUE))
        )

      )

    }
  })


  output$pie <- renderPlotly({
    # pg$core_level <- input$core_level
    updateCoreLevel()
    updateOrganisms()
    st <- as.data.frame(pg$summary_stats[-1, ])
    st$Category <- factor(st$Category, levels = c("Cloud", "Shell", "Core"))
    plot_ly(st, labels = ~Category, values = ~Number,
            type = "pie",
            showlegend=FALSE,
            textposition = 'inside', height = 350)
  })


  accs_pm <- eventReactive({
    # input$core_level
    # input$accs_level
    input$accs_freq
  }, {
    updateOrganisms()
    pm <- pg$pan_matrix
    norgs <- length(pg$organisms$org)
    pmb <- pm
    pmb[which(pmb>1L, arr.ind = TRUE)] <- 1L
    accs_freq <- input$accs_freq
    accs_num <- round(accs_freq *  norgs / 100)
    clsu <- colSums(pmb)
    wh <- which( clsu >= min(accs_num) & clsu <= max(accs_num))
    colnames(pmb)[wh]
    pm[, wh, drop = FALSE]
  }, ignoreNULL = TRUE)


  output$barplot <- renderPlotly({
    updateOrganisms()
    pm <- pg$pan_matrix
    pm[which(pm > 0, arr.ind = TRUE)] <- 1L
    dd <- as.data.frame(table(colSums(pm)))
    colnames(dd) <- c("Count", "Frequency")
    plot_ly(dd,
            x = ~Count, y = ~Frequency,
            type = "bar", height = 350) %>%
      layout(xaxis = list(title = "Number of genomes"),
             yaxis = list(title = 'Number of genes'))
  })

  output$curves <- renderPlotly({
    updateOrganisms()
    updateCoreLevel()
    what <- c('pangenome', 'coregenome')
    names(what) <- what
    lrar <- lapply(what, function(x) pg$rarefact(what = x))
    lfun <- lapply(what, function(x){
      if (x == 'pangenome'){
        pg$pg_power_law_fit(raref = lrar[[x]])#$formula
      }else{
        pg$cg_exp_decay_fit(raref = lrar[[x]])#$formula
      }
    })
    lrarm <- lapply(lrar, melt)
    ll <- lapply(what, function(x) {
      lrarm[[x]]$category <- x
      lrarm[[x]]
    })
    df <- Reduce(rbind, ll)

    norgs <- dim(pg$organisms)[1]
    interv <- (norgs - 1) / 512
    interp <- seq(1, norgs, interv)

    plot_ly(df,
            x=~Var1,
            y=~value,
            color = ~category,
            type="scatter",
            mode="markers",
            marker = list(opacity = 0.3),
            showlegend = FALSE,
            height = 350) %>%
      add_lines(x = interp, y = lfun$pangenome$formula(interp), color = "pangenome") %>%
      add_lines(x = interp, y = lfun$coregenome$formula(interp), color = "coregenome")
  })


  output$core_clusters <- DT::renderDT({
    updateOrganisms()
    updateCoreLevel()
    df <- as.data.frame(pg$core_clusters)
    datatable(df,
              selection = list(mode = "single", selected = 1),
              options = list(
                rownames = FALSE,
                paging = FALSE,
                # autoWidth = T,
                scrollX = T,
                scrollY = "200px",
                # scrollCollapse = T,
                pageLength = 3,
                dom ='ft'
              ))
  })

  output$core_genes <- DT::renderDT({
    updateOrganisms()
    updateCoreLevel()
    selected_cluster <- input$core_clusters_rows_selected
    df <- as.data.frame(pg$core_genes[[selected_cluster]])
    datatable(df,
              selection = "none",
              options = list(
                rownames = FALSE,
                paging = FALSE,
                scrollX = TRUE,
                scrollY = "200px",
                dom ='ft'
                # pageLength = 5
              ))
  })



  #############
  # ACCESSORY #
  #############

  output$accs_slider <- renderUI({
    # corelev <- input$core_level
    sliderInput("accs_freq",
                "Accessory frequency (%)",
                min = 0,
                max = accs_level(),
                value=c(0, accs_level()))
  })

  output$pca_x <- renderUI({
    updateOrganisms()
    norgs <- dim(pg$organisms)[1]
    selectInput("xpc",
                label = "X-axis PC",
                choices = paste0("PC", seq_len(norgs)),
                selected = "PC1",
                multiple = FALSE)
  })

  output$pca_y <- renderUI({
    updateOrganisms()
    norgs <- dim(pg$organisms)[1]
    selectInput("ypc",
                label = "Y-axis PC",
                choices = paste0("PC", seq_len(norgs)),
                selected = "PC2",
                multiple = FALSE)
  })


  pca <- reactive({
    pm <- accs_pm()
    prcomp(pm)
  })


  output$pca_color_meta <- renderUI({
    updateOrganisms()
    ch <- colnames(pg$organisms)[-1]
    selectInput("color_meta_selected",
                label = "Select metadata",
                choices = ch,
                selected = NULL)
  })



  output$pca <- renderPlotly({
    updateOrganisms()
    df <- as.data.frame(pca()$x)
    xax <- input$xpc
    yax <- input$ypc
    df$xx <- df[[xax]]
    df$yy <- df[[yax]]

    opts <- list()
    if (input$check_color_meta){
      cls <- pg$organisms[[input$color_meta_selected]]
      ln <- length(levels(cls))
      clrs <- colorSet(ln)
      opts <- c(opts, list(color = cls, colors = clrs))
    }

    opts <- c(opts, list(
      data = df,
      x = ~xx,
      y = ~yy,
      text = ~rownames(df),
      mode = "markers",
      marker = list(size = 12),
      type = "scatter"
    ))

    p <- do.call(plot_ly, opts)

    p %>% layout(xaxis = list(title = xax), yaxis = list(title = yax))

  })




  output$pca_summary <- DT::renderDT({
    updateOrganisms()
    PCA <- pca()
    sumpca <- summary(PCA)$importance
    datatable(sumpca,
              selection = "none",
              extensions = "FixedColumns",
              options = list(
                pageLength = 3,
                lengthChange = FALSE,
                searching = FALSE,
                scrollX = TRUE,
                dom = "t",
                # bsort = FALSE,
                fixedColumns = list(leftColumns = 1),
                autoWidth = TRUE,
                columnDefs = list(list(width = '150px', targets = 0),
                                  list(width = '40px', targets = 1:dim(sumpca)[2]))
              )) %>% formatSignif(columns = T, digits = 3)
  }, rownames = TRUE)



  colorSet <- function(n) {
    colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                "#66A61E", "#E6AB02", "#A6761D", "#666666")
    colors[seq_len(n)]
  }


  output$heat_freq <- renderPlotly({
    updateOrganisms()

    opts <- list()
    if (input$check_color_meta){
      cls <- as.data.frame(pg$organisms[input$color_meta_selected])
      rownames(cls) <- as.character(pg$organisms$org)
      opts <- c(opts, list(row_side_colors = cls, row_side_palette=colorSet))
    }

    # wh <- updateClusterList()
    # pm <- pg$pan_matrix[, wh, drop = FALSE]
    pm <- accs_pm()
    pm[which(pm>1L, arr.ind = TRUE)] <- 1L

    # wh <- updateClusterList()
    # pm <- pm[, wh, drop=FALSE]
    csm <- colSums(pm)
    spl <- split(names(csm), csm)
    tpm <- t(pm)
    norder <- lapply(spl, function(x) {
      if (length(x)<2){
        x
      }else{
        d <- vegdist(tpm[x, , drop=F], method = "jaccard", binary = T, na.rm = T)
        hc <- hclust(d, "single")
        x[hc$order]
      }
    })
    norder <- norder[order(as.integer(names(norder)), decreasing = T)]
    forder <- unlist(norder)
    pm <- pm[, forder,  drop = F]

    d <- vegdist(pm, method = "jaccard")
    hc <- as.dendrogram(hclust(d))
    dend <- dendextend::seriate_dendrogram(hc, d)

    ## HEATMAP
    opts <- c(opts,
              list(x=pm,
                   colors = Blues,
                   Colv = FALSE,
                   Rowv = dend,
                   margin = c(0,0,0,0),
                   # distfun = vegan::vegdist,
                   # dist_method = "jaccard",
                   hide_colorbar = TRUE,
                   showticklabels = FALSE))

    p1 <- do.call(heatmaply, opts)

    ## BARPLOT
    rsum <- rowSums(pm)
    odend <- rev(order.dendrogram(dend))
    bar <- data.frame(org = factor(names(rsum),
                                   levels = names(rsum)[odend]),
                      Count = rsum,
                      row.names = NULL)
    # bar <- rsdf[hc$order, ]
    p2 <- plot_ly(#data = bar,
      y = ~bar$org,
      x = ~bar$Count,
      type = "bar",
      orientation="h") %>%
      layout(xaxis = list(autorange='reversed', showticklabels=FALSE),
             yaxis = list(showticklabels = FALSE))

    subplot(p2, p1, widths = c(.3,.7))

  })

  output$accs_clusters <- DT::renderDT({
    updateOrganisms()
    # wh <- updateClusterList()
    wh <- colnames(accs_pm())
    clust <- pg$clusters
    ma <- match(wh, clust$group)
    df <- as.data.frame(clust[ma, ,drop=FALSE])
    datatable(df,
              selection = list(mode = "single", selected = 1),
              options = list(
                rownames = FALSE,
                paging = FALSE,
                # autoWidth = T,
                scrollX = T,
                scrollY = "200px",
                # scrollCollapse = T,
                pageLength = 3,
                dom ='ft'
              ))
  })

  output$accs_genes <- DT::renderDT({
    updateOrganisms()
    wh <- colnames(accs_pm())
    selected_cluster <- wh[input$accs_clusters_rows_selected]
    df <- as.data.frame(pg$genes[[selected_cluster]])
    datatable(df,
              selection = "none",
              options = list(
                rownames = FALSE,
                paging = FALSE,
                scrollX = TRUE,
                scrollY = "200px",
                dom ='ft'
                # pageLength = 5
              ))
  })

}
