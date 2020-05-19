# PgR6M

#' @name PgR6M
#' @title PgR6 class with methods.
#' @description PgR6 with Methods. Final users should use \code{\link{pagoo}}
#' instead of this, since is more easy to understand.
#' Inherits: \code{\link[pagoo]{PgR6}}
#' @importFrom R6 R6Class
#' @importFrom micropan fluidity binomixEstimate
#' @import ggplot2
#' @import ggfortify
#' @import shiny
#' @import shinyWidgets
#' @import shinydashboard
#' @importFrom DT DTOutput renderDT datatable formatSignif JS
#' @importFrom plotly plot_ly plotlyOutput renderPlotly layout add_lines subplot
#' @importFrom heatmaply heatmaply
#' @importFrom reshape2 melt
#' @importFrom vegan vegdist
#' @importFrom dendextend seriate_dendrogram
#' @importFrom magrittr `%>%`
#' @export
PgR6M <- R6Class('PgR6M',

                 inherit = PgR6,

                 public = list(

                   #' @description
                   #' Create a \code{PgR6M} object.
                   #' @param data A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing at least the
                   #' following columns: \code{gene} (gene name), \code{org} (organism name to which the gene belongs to),
                   #' and \code{cluster} (group of orthologous to which the gene belongs to). More columns can be added as metadata
                   #' for each gene.
                   #' @param org_meta (optional) A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
                   #' containging additional metadata for organisms. This \code{data.frame} must have a column named "org" with
                   #' valid organisms names (that is, they should match with those provided in \code{data}, column \code{org}), and
                   #' additional columns will be used as metadata. Each row should correspond to each organism.
                   #' @param cluster_meta (optional) A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
                   #' containging additional metadata for clusters. This \code{data.frame} must have a column named "cluster" with
                   #' valid organisms names (that is, they should match with those provided in \code{data}, column \code{cluster}), and
                   #' additional columns will be used as metadata. Each row should correspond to each cluster.
                   #' @param core_level The initial core_level (thats the percentage of organisms a core cluster must be in to be
                   #' considered as part of the core genome). Must be a number between 100 and 85, (default: 95). You can change it
                   #' later by using the \code{$core_level} field once the object was created.
                   #' @param sep A separator. By default is '__'(two underscores). It will be used to
                   #' create a unique \code{gid} (gene identifier) for each gene. \code{gid}s are created by pasting
                   #' \code{org} to \code{gene}, separated by \code{sep}.
                   #' @param verbose \code{logical}. Whether to display progress messages when loading class.
                   #' @param DF Deprecated. Use \code{data} instead.
                   #' @param group_meta Deprecated. Use \code{cluster_meta} instead.
                   #' @return An R6 object of class PgR6M. It contains basic fields and methods for analyzing a pangenome. It also
                   #' contains additional statistical methods for analize it, and methods to make basic
                   #' exploratory plots.
                   initialize = function(data,
                                         org_meta,
                                         cluster_meta,
                                         core_level = 95,
                                         sep = '__',
                                         verbose = TRUE,
                                         DF,
                                         group_meta){

                     # Deprecated args
                     if (!missing(DF)){
                       data <- DF
                       warning('Argument "DF" is deprecated. You should use "data" instead. ',
                               'This still exists for compatibility with previous versions, ',
                               'but will throw an error in the future.',
                               immediate. = TRUE, noBreaks. = TRUE)
                     }

                     if (!missing(group_meta)){
                       cluster_meta <- group_meta
                       warning('Argument "group_meta" is deprecated. You should use "cluster_meta" instead. ',
                               'This still exists for compatibility with previous versions, ',
                               'but will throw an error in the future.',
                               immediate. = TRUE, noBreaks. = TRUE)
                     }

                     super$initialize(data = data,
                                      org_meta,
                                      cluster_meta,
                                      core_level = core_level,
                                      sep = sep,
                                      verbose = verbose)

                   },

                   #Analyses Methods
                   #' @description
                   #' Rarefact pangenome or corgenome. Compute the number of genes which belong to
                   #' the pangenome or to the coregenome, for a number of random permutations of
                   #' increasingly bigger sample of genomes.
                   #' @param what One of \code{"pangenome"} or \code{"coregenome"}.
                   #' @param n.perm The number of permutations to compute (default: 10).
                   #' @return A \code{matrix}, rows are the number of genomes added, columns are
                   #' permutations, and the cell number is the number of genes in ethier category.
                   rarefact = function(what = 'pangenome', n.perm = 10){
                     what <- match.arg(what, c('pangenome', 'coregenome'))
                     pm <- self$pan_matrix
                     pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                     norgs <- dim(self$organisms)[1]
                     rmat <- matrix(0L, nrow = norgs, ncol = n.perm)
                     if (what=='pangenome'){
                       for (i in seq_len(n.perm)){
                         cm <- apply(pm[sample(norgs), ], 2, cumsum)
                         rmat[, i] <- rowSums(cm > 0)
                       }
                     }else{
                       sq <- seq_len(norgs)
                       for (i in seq_len(n.perm)){
                         cm <- apply(pm[sample(norgs), ], 2, cumsum)
                         cr <- apply(cm, 2, `==`, sq)
                         rmat[, i] <- rowSums(cr)
                       }
                     }
                     rownames(rmat) <- seq_len(norgs)
                     colnames(rmat) <- paste0('permut_', seq_len(n.perm))
                     rmat
                   },

                   #' @description
                   #' Compute distance between all pairs of genomes. The default dist method is
                   #' \code{"bray"} (Bray-Curtis distance). Annother used distance method is \code{"jaccard"},
                   #' but you should set \code{binary = FALSE} (see below) to obtain a meaningful result.
                   #' See \code{\link[vegan]{vegdist}} for details, this is just a wrapper function.
                   #' @param method The distance method to use. See \link[vegan]{vegdist}
                   #' for available methods, and details for each one.
                   #' @param binary Transform abundance matrix into a presence/absence
                   #' matrix before computing distance.
                   #' @param diag Compute diagonals.
                   #' @param upper Return only the upper diagonal.
                   #' @param na.rm Pairwise deletion of missing observations when
                   #' computing dissimilarities.
                   #' @param ... Other parameters. See \link[vegan]{vegdist} for details.
                   #' @return A \code{dist} object containing all pairwise dissimilarities between genomes.
                   dist = function(method = 'bray',
                                   binary = FALSE,
                                   diag = FALSE,
                                   upper = FALSE,
                                   na.rm = FALSE,
                                   ...){

                     METHODS <- c("manhattan","euclidean","canberra","bray",
                                  "kulczynski","jaccard","gower","altGower",
                                  "morisita","horn","mountford","raup",
                                  "binomial","chao","cao","mahalanobis")
                     method <- match.arg(method, METHODS)

                     if (method == 'jaccard' & binary == FALSE){
                       warning('It is recommended to set binary = TRUE when running dist(method = "jaccard")')
                     }

                     #vegan::vegdist()
                     vegan::vegdist(self$pan_matrix,
                                    method = method,
                                    diag = diag,
                                    upper = upper,
                                    na.rm = na.rm,
                                    ...)
                   },

                   #' @description
                   #' Performs a principal components analysis on the panmatrix
                   #' @param center a logical value indicating whether the variables should be shifted
                   #' to be zero centered. Alternately, a vector of length equal the number of columns of x can be
                   #' supplied. The value is passed to scale.
                   #' @param scale.  a logical value indicating whether the variables should be scaled
                   #' to have unit variance before the analysis takes place. The default is TRUE.
                   #' @param ... Other arguments. See \link[stats]{prcomp}
                   #' @return Returns a list with class "prcomp". See \link[stats]{prcomp} for more information.
                   pan_pca = function(center = TRUE, scale. = FALSE, ...){
                     prcomp(self$pan_matrix, center = center, scale. = scale., ...)
                   },

                   #' @description
                   #' Fits a power law curve for the pangenome rarefaction simulation.
                   #' @param raref (Optional) A rarefaction matrix, as returned by \code{rarefact()}.
                   #' @param ... Further arguments to be passed to \code{rarefact()}. If \code{raref}
                   #' is missing, it will be computed with default arguments, or with the ones provided here.
                   #' @return A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
                   #' with fitted parameters. An attribute \code{"alpha"} is also returned (If
                   #' \code{alpha>1}, then the pangenome is closed, otherwise is open.
                   pg_power_law_fit = function(raref, ...){
                     # #micropan::heaps()
                     # heaps(self$pan_matrix ,n.perm = n.perm)
                     if (missing(raref)){
                       raref <- self$rarefact(...)
                     }
                     rm <- melt(raref)
                     # Power law linearization:
                     # y = K * x ^ delta ==> log(y) = log(K) + delta * log(x)
                     fitHeaps <- lm(log(rm$value) ~ log(rm$Var1))
                     logK <- summary(fitHeaps)$coef[1]
                     delta <- summary(fitHeaps)$coef[2]
                     K <- exp(logK)
                     alpha <- 1-delta
                     ret <- list(formula=NULL,
                                 params=NULL)
                     ret$formula <- function(x) K * x ^ delta
                     ret$params <- c(K = K, delta = delta)
                     attr(ret, 'alpha') <- alpha
                     ret
                   },

                   #' @description
                   #' Fits an exponential decay curve for the coregenome rarefaction simulation.
                   #' @param raref (Optional) A rarefaction matrix, as returned by \code{rarefact()}.
                   #' @param pcounts An integer of pseudo-counts. This is used to better fit the function
                   #' at small numbers, as the linearization method requires to substract a constant C, which is the
                   #' coregenome size, from \code{y}. As \code{y} becomes closer to the coregenome size, this operation
                   #' tends to 0, and its logarithm goes crazy. By default \code{pcounts=10}.
                   #' @param ... Further arguments to be passed to \code{rarefact()}. If \code{raref}
                   #' is missing, it will be computed with default arguments, or with the ones provided here.
                   #' @return A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
                   #' with fitted intercept and decay parameters.
                   cg_exp_decay_fit = function(raref, pcounts = 10, ...){
                     # Exponential decay linearization:
                     # y = A * exp(K * t) + C ==> y - C = K*t + log(A)
                     # C == core size
                     if (missing(raref)){
                       raref <- self$rarefact(what = 'core', ...)
                     }
                     rm <- melt(raref)
                     C <- min(rm$value)
                     fitExp <- lm(log(value - C + pcounts) ~ Var1, data= rm)
                     A <- exp(summary(fitExp)$coef[1])
                     B <- summary(fitExp)$coef[2]
                     ret <- list(formula=NULL,
                                 params=NULL)
                     ret$formula <- function(x) A * exp(B * x) + C
                     ret$params <- c(A = A, B = B, C = C)
                     attr(ret, 'pseudo_counts') <- pcounts
                     ret
                   },

                   #' @description
                   #' Computes the genomic fluidity, which is a measure of population
                   #' diversity. See \code{\link[micropan]{fluidity}} for more details.
                   #' @param n.sim An integer specifying the number of random samples
                   #' to use in the computations.
                   #' @return A list with two elements, the mean fluidity and its sample standard
                   #' deviation over the n.sim computed values.
                   fluidity = function(n.sim = 10){
                     #micropan::fluidity()
                     fluidity(self$pan_matrix, n.sim = n.sim)
                   },

                   #' @description
                   #' Fits binomial mixture models to the data given as a pan-matrix. From the
                   #' fitted models both estimates of pan-genome size and core-genome size are
                   #' available. See \code{\link[micropan]{binomixEstimate}} for more details.
                   #' @param K.range The range of model complexities to explore. The
                   #' vector of integers specify the number of binomial densities to combine in the
                   #' mixture models.
                   #' @param core.detect.prob The detection probability of core genes.
                   #' This should almost always be 1.0, since a core gene is by definition always
                   #' present in all genomes, but can be set fractionally smaller.
                   #' @param verbose  Logical indicating if textual output should be
                   #' given to monitor the progress of the computations.
                   #' @return A \code{Binomix} object (\code{micropan} package), which is a small (S3)
                   #' extension of a \code{list} with two components. These two components are named
                   #' \code{BIC.table} and \code{Mix.list}. Refer to the \code{micropan} function
                   #' \code{\link[micropan]{binomixEstimate}} for a detailed explanation of this
                   #' method.
                   binomix_estimate = function(K.range = 3:5,
                                               core.detect.prob = 1,
                                               verbose = TRUE){
                     # micropan::binomxEstimate()
                     binomixEstimate(pan.matrix = self$pan_matrix,
                                     K.range = K.range,
                                     core.detect.prob = core.detect.prob,
                                     verbose = verbose)
                   },


                   #Plot Methods
                   #' @description
                   #' Plot a barplot with the frequency of genes within the total number of
                   #' genomes.
                   #' @return A barplot, and a \code{gg} object (\code{ggplot2} package) invisibly.
                   gg_barplot = function(){
                     pm <- self$pan_matrix
                     pm[which(pm > 0, arr.ind = TRUE)] <- 1L
                     dd <- as.data.frame(table(colSums(pm)))
                     dd$Var1 <- as.integer(as.character(dd$Var1))
                     ggplot(dd, aes(x=Var1, y=Freq)) +
                       geom_bar(stat='identity') +
                       ylab('Number of genes') +
                       xlab('Number of genomes')
                   },

                   #' @description
                   #' Plot a pangenome binary map representing the presence/absence of each
                   #' gene within each organism.
                   #' @return A binary map (\code{ggplot2::geom_raster()}), and a \code{gg} object (\code{ggplot2}
                   #' package) invisibly.
                   gg_binmap = function(){
                     tpm <- t(self$pan_matrix)
                     tpm[which(tpm > 0, arr.ind = TRUE)] <- 1L
                     bm <- as.data.frame(tpm)
                     or <- order(rowSums(bm), decreasing = TRUE)
                     lvls <- rownames(bm)[or]
                     bm$Cluster <- factor(rownames(bm), levels = lvls)
                     bm <- melt(bm, 'Cluster')
                     bm$value <- factor(bm$value, levels = c(1, 0))
                     colnames(bm)[which(colnames(bm) == 'variable')] <- "Organism"
                     ggplot(bm, aes(Cluster, Organism, fill=value)) +
                       geom_raster() +
                       theme(axis.ticks.x = element_blank(),
                             axis.text.x = element_blank()) +
                       scale_fill_grey(start = .2, end = .9)
                   },


                   #' @description
                   #' Plot a heatmap showing the computed distance between all pairs of organisms.
                   #' @param method Distance method. One of "Jaccard" (default), or "Manhattan",
                   #' see above.
                   #' @param ... More arguments to be passed to \code{\link[micropan]{distManhattan}}.
                   #' @return A heatmap (\code{ggplot2::geom_tile()}), and a \code{gg} object (\code{ggplot2}
                   #' package) invisibly.
                   gg_dist = function(method = 'bray', ...){
                     m <- as.matrix(self$dist(method = method, ...))
                     me <- melt(m)
                     ggplot(me, aes(Var1, Var2, fill=value)) +
                       geom_tile()
                   },


                   #' @description
                   #' Plot a scatter plot of a Principal Components Analysis.
                   #' @param colour The name of the column in \code{$organisms} field from which points will take
                   #' colour (if provided). \code{NULL} (default) renders black points.
                   #' @param ... More arguments to be passed to \code{ggplot2::autoplot()}.
                   #' @return A scatter plot (\code{ggplot2::autoplot()}), and a \code{gg} object (\code{ggplot2}
                   #' package) invisibly.
                   gg_pca = function(colour = NULL, ...){
                     ocol <- colnames(self$organisms)
                     ocol <- ocol[ocol!='org']
                     color <- match.arg(colour, c(ocol, NULL))
                     pca <- self$pan_pca(...)
                     autoplot(pca,
                              data = as.data.frame(self$organisms),
                              colour = colour,
                              ...)
                   },


                   #' @description
                   #' Plot a pie chart showing the number of clusters of each pangenome category: core,
                   #' shell, or cloud.
                   #' @return A pie chart (\code{ggplot2::geom_bar() + coord_polar()}), and a \code{gg} object
                   #' (\code{ggplot2} package) invisibly.
                   gg_pie = function(){
                     st <- self$summary_stats[-1, ]
                     # st <- st[-1, ]
                     st$Category <- factor(st$Category,
                                           levels = c("Cloud", "Shell", "Core"))
                     ggplot(as.data.frame(st), aes(x='', y=Number, fill=Category)) +
                       geom_bar(width = 1, stat = "identity") +
                       coord_polar("y", start=0)
                   },

                   #' @description
                   #' Plot pangenome and/or coregenome curves with the fitted functions returned by \code{pg_power_law_fit()}
                   #'  and \code{cg_exp_decay_fit()}. You can add points by adding \code{+ geom_points()}, of ggplot2 package
                   #' @param what One of \code{"pangenome"} or \code{"coregenome"}.
                   #' @param ... ????
                   #' @return A scatter plot, and a \code{gg} object (\code{ggplot2} package) invisibly.
                   gg_curves = function(what = c('pangenome', 'coregenome'),
                                        ...){
                     what <- match.arg(what, c('pangenome', 'coregenome'), several.ok = TRUE)
                     names(what) <- what
                     lrar <- lapply(what, function(x) self$rarefact(what = x))
                     lfun <- lapply(what, function(x){
                       if (x == 'pangenome'){
                         self$pg_power_law_fit(raref = lrar[[x]])#$formula
                       }else{
                         self$cg_exp_decay_fit(raref = lrar[[x]])#$formula
                       }
                     })
                     lrarm <- lapply(lrar, melt)
                     ll <- lapply(what, function(x) {
                       olst <- list(data = NULL, formula = NULL)
                       lrarm[[x]]$Category <- x
                       lrarm[[x]]
                      })
                     data <- Reduce(rbind, ll)

                     #plot
                     g <- ggplot(data, aes(x=Var1, y=value, colour=Category)) +
                       xlab('Number of genomes') +
                       ylab('Number of clusters')
                     for (i in seq_along(what)){
                       g <- g +
                         stat_function(data = data[which(data$Category == what[i]), ],
                                       fun = lfun[[what[i]]]$formula)
                     }
                     g
                   },





                   #################
                   ##             ##
                   ##  SHINY APP  ##
                   ##             ##
                   #################


                   #' @description
                   #' Launch an interactive shiny app. It contains a sidebar
                   #' with controls and switches to interact with the pagoo
                   #' object. You can drop/recover organisms from the dataset,
                   #' modify the core_level, visualize statistics, plots, and
                   #' browse cluster and gene information. In the main body, it
                   #' contains 2 tabs to switch between summary statistics plots
                   #' and core genome information on one side, and accessory
                   #' genome plots and information on the other.
                   #'
                   #' The lower part of each tab contains two tables, side by
                   #' side. On the "Summary" tab, the left one contain
                   #' information about core clusters, with one cluster per row.
                   #' When one of them is selected (click), the one on the right
                   #' is updated to show information about its genes (if
                   #' provided), one gene per row. On the "Accessory" tab, a
                   #' similar configuration is shown, but on this case only
                   #' accessory clusters/genes are displayed. There is a slider
                   #' on the sidebar where one can select the accessory
                   #' frequency range to display.
                   #'
                   #' Give it a try!
                   #'
                   #' Take into account that big pangenomes can slow down the
                   #' performance of the app. More than 50-70 organisms often
                   #' leads to a delay in the update of the plots/tables.
                   runShinyApp = function(){

                     pg <- self

                     app <- list(ui = dashboardPage(

                       header = dashboardHeader(
                         title = "Pagoo Shiny App",
                         titleWidth = 250
                       ),


                       sidebar = dashboardSidebar(width = 250,

                                                  sidebarMenu(

                                                    # uiOutput("select_organisms"),
                                                    pickerInput(inputId = "organisms",
                                                                label = "Organisms",
                                                                choices = as.character(pg$.__enclos_env__$private$.organisms$org),
                                                                selected = as.character(pg$organisms$org),
                                                                multiple = TRUE,
                                                                options = list(`actions-box` = TRUE),
                                                                width = "100%"),

                                                    # uiOutput("select_variable"),
                                                    pickerInput(inputId = "variable",
                                                                label = "Variable",
                                                                choices = c("\a", colnames(pg$.__enclos_env__$private$.organisms)[-1]),
                                                                # selected = all_vars[1],
                                                                multiple = FALSE,
                                                                width = "100%"),

                                                    uiOutput("select_category"),

                                                    sliderInput("core_level",
                                                                "Core Level",
                                                                min=85,
                                                                max=100,
                                                                value=95),

                                                    uiOutput(outputId = "accs_slider"),

                                                    selectInput("color_meta_selected",
                                                                label = "Color by (select metadata):",
                                                                choices = c("\a", colnames(pg$.__enclos_env__$private$.organisms)[-1]))

                                                  )
                       ),



                       body = dashboardBody(

                         tags$head(tags$style(HTML('.info-box {min-height: 45px;}
                                                   .info-box-icon {height: 45px; line-height: 45px;}
                                                   .info-box-content {padding-top: 0px; padding-bottom: 0px;}
                                                   .selectize-dropdown {z-index: 10000 !important;}'))),

                         tags$script(HTML("$('body').addClass('fixed');")),

                         fluidRow(
                           infoBoxOutput("num_orgs"),
                           infoBoxOutput("num_clus"),
                           infoBoxOutput("num_gene")
                         ),
                         fluidRow(
                           tabBox(width = 12,

                                  tabPanel(
                                    title = "Summary",

                                    fluidRow(
                                      box(
                                        title = "Core number",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 3,
                                        height = 500,
                                        # align = "center",
                                        # DT::DTOutput("summary_counts"),
                                        addSpinner(plotlyOutput("core_evol", width = "auto"))
                                      ),
                                      box(
                                        title = "Summary",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 3,
                                        height = 500,
                                        #offset = 0,
                                        addSpinner(plotlyOutput("pie"))
                                      ),
                                      box(
                                        title = "Frequency plot",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 3,
                                        height = 500,
                                        #offset = 0,
                                        addSpinner(plotlyOutput("barplot"))
                                      ),
                                      box(
                                        title = "Pangenome curves",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 3,
                                        height = 500,
                                        #offset = 0,
                                        addSpinner(plotlyOutput("curves"))
                                      )
                                    ),

                                    # hr(),

                                    fluidRow(
                                      box(
                                        title = "Core clusters",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 6,
                                        addSpinner(DT::DTOutput("core_clusters"))
                                      ),
                                      box(
                                        title = "Core genes",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 6,
                                        addSpinner(DT::DTOutput("core_genes"))
                                      )

                                    )


                                  ), # Ends tabPanel "Summary"

                                  tabPanel(
                                    title = "Accessory",

                                    fluidRow(

                                      box(
                                        title = "Gene presence/absence",
                                        status = "primary",
                                        solidHeader = T,
                                        height = 500,
                                        addSpinner(plotlyOutput("heat_freq", height = 420))
                                      ),

                                      box(
                                        title = "PCA",
                                        status = "primary",
                                        solidHeader = T,
                                        height = 500,
                                        splitLayout(
                                          addSpinner(DT::DTOutput("pca_summary")),
                                          column(width = 12,
                                                 uiOutput("pca_x"),
                                                 uiOutput("pca_y"),
                                                 style="z-index:1002;"
                                          ),
                                          cellWidths = c("80%", "20%")
                                        ),
                                        box(addSpinner(plotlyOutput("pca", height = 250)), width = 12)
                                      ),

                                    ),

                                    fluidRow(
                                      box(
                                        title = "Accessory clusters",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 6,
                                        addSpinner(DT::DTOutput("accs_clusters"))
                                      ),
                                      box(
                                        title = "Accessory genes",
                                        status = "primary",
                                        solidHeader = T,
                                        width = 6,
                                        DT::DTOutput("accs_genes")
                                      )

                                    )

                                  ) # Ends tabPanel "Accessory"

                           )
                         )
                       )

                       # Ends tabPanel

                     ),


                     server = function(input, output, session){


                       #############
                       ## OPTIONS ##
                       #############

                       dataorg <- as.data.frame(pg$.__enclos_env__$private$.organisms)
                       all_orgs <- as.character(dataorg$org)
                       orgs <- reactiveVal(as.character(pg$organisms$org))
                       all_vars <- colnames(dataorg)[-1]

                       output$select_category <- renderUI(
                         pickerInput(inputId = "category",
                                     label = "Category",
                                     choices = "\a",
                                     selected = "\a",
                                     multiple = TRUE,
                                     options = list(`actions-box` = TRUE),
                                     width = "100%")
                       )


                       observeEvent(input$organisms, {
                         orgs(input$organisms)
                       })

                       updateOrganisms <- reactive({
                         ava_orgs <- as.character(pg$organisms[["org"]])
                         inp_orgs <- orgs()

                         toDrop <- ava_orgs[! ava_orgs %in% inp_orgs]
                         toReco <- inp_orgs[! inp_orgs %in% ava_orgs]

                         pg$drop(toDrop)
                         pg$recover(toReco)
                       })

                       observeEvent(input$variable, {
                         # updateOrganisms()
                         meta <- unique(as.character(dataorg[[req(input$variable)]]))
                         updatePickerInput(session = session,
                                           inputId = "category",
                                           label = "Category",
                                           choices = meta,
                                           selected = meta)
                       })

                       observeEvent(input$category, {
                         # updateOrganisms()
                         var <- input$variable
                         # orgs <- as.character(pg$organisms$org)
                         nw <- all_orgs[dataorg[[var]] %in% req(input$category)]
                         updatePickerInput(
                           session = session,
                           inputId = "organisms",
                           label = "Organisms",
                           choices = all_orgs,
                           selected = nw
                         )
                       })

                       updateCoreLevel <- eventReactive(input$core_level, {
                         pg$core_level <- input$core_level
                         pg
                       })

                       accs_level <- reactiveVal(value = pg$core_level - 1L)

                       observeEvent(input$core_level, {
                         accs_level(input$core_level-1L)
                       })

                       output$num_orgs <- renderInfoBox({
                         updateOrganisms()
                         x <- dim(pg$organisms)[1]
                         infoBox("Number of organisms", x)
                       })

                       output$num_clus <- renderInfoBox({
                         updateOrganisms()
                         x <- dim(pg$clusters)[1]
                         infoBox("Number of clusters", x, color = "purple")
                       })

                       output$num_gene <- renderInfoBox({
                         updateOrganisms()
                         x <- sum(pg$pan_matrix)
                         infoBox("Number of genes", x, color = "yellow")
                       })



                       #############
                       ## SUMMARY ##
                       #############



                       output$core_evol <- renderPlotly({
                         updateOrganisms()
                         pm <- pg$pan_matrix
                         pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                         sq <- seq(1, 0.85, -0.01)
                         cs <- colSums(pm)
                         norg <- dim(pm)[1]
                         ev <- sapply(sq, function(x){
                           length(which( cs >= floor(norg*x) ))
                         })
                         data <- data.frame(core_level = sq*100, core_number = ev)
                         plot_ly(data, x = ~core_level, y = ~core_number,
                                 type = "scatter",
                                 mode = "markers", height = 420) %>%
                           layout(xaxis = list(title = "Core Level"),
                                  yaxis = list(title = "Core Number"))
                       })


                       output$pie <- renderPlotly({
                         # pg$core_level <- input$core_level
                         updateCoreLevel()
                         updateOrganisms()
                         st <- as.data.frame(pg$summary_stats[-1, ])
                         st$Category <- factor(st$Category, levels = c("Cloud", "Shell", "Core"))
                         plot_ly(st, labels = ~Category, values = ~Number,
                                 type = "pie",
                                 showlegend=FALSE,height = 400,
                                 textposition = 'inside') %>%
                           layout(margin = list(
                             l = 20,
                             r = 20,
                             b = 0,
                             t = 20,
                             pad = 4
                           ))
                       })

                       # Initialize accs_pm in order to not have to wait until user clicks on
                       # accessory tab.
                       pm <- pg$pan_matrix
                       norgs <- length(pg$organisms$org)
                       pmb <- pm
                       pmb[which(pmb>1L, arr.ind = TRUE)] <- 1L
                       accs_freq <- c(0, pg$core_level-1)
                       accs_num <- round(accs_freq *  norgs / 100)
                       clsu <- colSums(pmb)
                       wh <- which( clsu >= min(accs_num) & clsu <= max(accs_num))
                       colnames(pmb)[wh]
                       accs_pm_i <- pm[, wh, drop = FALSE]
                       accs_pm <- reactiveVal(value = accs_pm_i)

                       observeEvent({
                         input$accs_freq
                         input$organisms
                       },{
                         updateOrganisms()
                         pm <- pg$pan_matrix
                         norgs <- length(pg$organisms$org)
                         pmb <- pm
                         pmb[which(pmb>1L, arr.ind = TRUE)] <- 1L
                         accs_freq <- req(input$accs_freq)
                         accs_num <- floor(accs_freq *  norgs / 100)
                         clsu <- colSums(pmb)
                         wh <- which( clsu >= min(accs_num) & clsu <= max(accs_num))
                         colnames(pmb)[wh]
                         accs_pm(pm[, wh, drop = FALSE])
                       }, ignoreNULL = TRUE)


                       output$barplot <- renderPlotly({
                         updateOrganisms()
                         pm <- pg$pan_matrix
                         pm[which(pm > 0, arr.ind = TRUE)] <- 1L
                         dd <- as.data.frame(table(colSums(pm)))
                         colnames(dd) <- c("Count", "Frequency")
                         plot_ly(dd,
                                 x = ~Count, y = ~Frequency,
                                 type = "bar", height = 420) %>%
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
                         data <- Reduce(rbind, ll)

                         norgs <- dim(pg$organisms)[1]
                         interv <- (norgs - 1) / 512
                         interp <- seq(1, norgs, interv)

                         plot_ly(data,
                                 x=~Var1,
                                 y=~value,
                                 color = ~category,
                                 type="scatter",
                                 mode="markers",
                                 marker = list(opacity = 0.3),
                                 showlegend = FALSE,
                                 height = 420) %>%
                           add_lines(x = interp, y = lfun$pangenome$formula(interp), color = "pangenome") %>%
                           add_lines(x = interp, y = lfun$coregenome$formula(interp), color = "coregenome") %>%
                           layout(xaxis = list(title = "Number of genomes"),
                                  yaxis = list(title = "Number of genes"))
                       })


                       output$core_clusters <- DT::renderDT({
                         updateOrganisms()
                         updateCoreLevel()
                         data <- as.data.frame(pg$core_clusters)
                         datatable(data,
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
                         selected_cluster <- req(input$core_clusters_rows_selected)
                         data <- as.data.frame(pg$core_genes[[selected_cluster]])
                         tgt <- which(sapply(data, function(x) max(nchar(as.character(x), allowNA = T, type = "width")) )>=26)
                         datatable(data,
                                   selection = "none",
                                   options = list(
                                     rownames = FALSE,
                                     paging = FALSE,
                                     scrollX = TRUE,
                                     scrollY = "200px",
                                     dom ='ft',
                                     autoWidth = TRUE,
                                     columnDefs = list(list(
                                       targets = unname(tgt),
                                       render = JS(
                                         "function(data, type, row, meta) {",
                                         "return type === 'display' && data.length > 30 ?",
                                         "'<span title=\"' + data + '\">' + data.substr(0, 26) + '...</span>' : data;",
                                         "}"),
                                       width = "200px"
                                     ))
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
                                     multiple = FALSE,
                                     width = "100px")
                       })

                       output$pca_y <- renderUI({
                         updateOrganisms()
                         norgs <- dim(pg$organisms)[1]
                         selectInput("ypc",
                                     label = "Y-axis PC",
                                     choices = paste0("PC", seq_len(norgs)),
                                     selected = "PC2",
                                     multiple = FALSE,
                                     width = "100px")
                       })


                       pca <- reactive({
                         pm <- accs_pm()
                         prcomp(pm)
                       })



                       output$pca <- renderPlotly({
                         updateOrganisms()
                         data <- as.data.frame(pca()$x)
                         xax <- req(input$xpc)
                         yax <- req(input$ypc)
                         data$xx <- data[[xax]]
                         data$yy <- data[[yax]]

                         opts <- list()
                         # color_by <- input$color_meta_selected
                         if (colorBy() != "\a"){
                           cls <- pg$organisms[[colorBy()]]
                           ln <- length(levels(cls))
                           clrs <- colorSet(ln)
                           opts <- c(opts, list(color = cls, colors = clrs))
                         }

                         opts <- c(opts, list(
                           data = data,
                           x = ~xx,
                           y = ~yy,
                           text = ~rownames(data),
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


                       colorBy <- reactiveVal("\a")
                       observeEvent(input$color_meta_selected, {
                         colorBy(input$color_meta_selected)
                       })

                       output$heat_freq <- renderPlotly({
                         updateOrganisms()
                         pm <- accs_pm()
                         pm[which(pm>1L, arr.ind = TRUE)] <- 1L

                         opts <- list()
                         # color_by <- input$color_meta_selected
                         # if (is.null(color_by)) color_by <- "\a"
                         if (colorBy() != "\a"){
                           cls <- as.data.frame(pg$organisms[colorBy()])
                           rownames(cls) <- as.character(pg$organisms$org)
                           opts <- c(opts, list(row_side_colors = cls, row_side_palette=colorSet))
                         }

                         # wh <- updateClusterList()
                         # pm <- pg$pan_matrix[, wh, drop = FALSE]

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
                                   list(x = pm,
                                        colors = c("#F7FBFF", "#08306B"),
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
                         # bar <- rsdata[hc$order, ]
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
                         ma <- match(wh, clust$cluster)
                         data <- as.data.frame(clust[ma, ,drop=FALSE])
                         datatable(data,
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
                         accs_rows <- req(input$accs_clusters_rows_selected)
                         wh <- colnames(accs_pm())
                         selected_cluster <- wh[accs_rows]
                         data <- as.data.frame(pg$genes[[selected_cluster]])
                         # chf <- which(lapply(data, class) %in% c("character", "factor"))
                         tgt <- which(sapply(data, function(x) max(nchar(as.character(x), allowNA = T, type = "width")) )>=26)
                         datatable(data,
                                   selection = "none",
                                   options = list(
                                     rownames = FALSE,
                                     paging = FALSE,
                                     scrollX = TRUE,
                                     scrollY = "200px",
                                     dom ='ft',
                                     autoWidth = TRUE,
                                     columnDefs = list(list(
                                       targets = unname(tgt),
                                       render = JS(
                                         "function(data, type, row, meta) {",
                                         "return type === 'display' && data.length > 30 ?",
                                         "'<span title=\"' + data + '\">' + data.substr(0, 26) + '...</span>' : data;",
                                         "}"),
                                       width = "200px"
                                     ))
                                   ))
                       })

                     })


                     runApp(app)

                   }


                 )

)
