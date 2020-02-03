# PgR6M

#' @name PgR6M
#' @title PgR6 class with methods.
#' @description PgR6 with Methods.
#' Inherits: \code{\link[pagoo]{PgR6}}
#'
#'
#' @section Class Constructor:
#' \describe{
#'     \item{\code{new(DF, org_meta, group_meta, sep = "__")}}{
#'         \itemize{
#'             \item{Create a \code{PgR6M} object.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{DF}}: A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing at least the
#'                     following columns: \code{gene} (gene name), \code{org} (organism name to which the gene belongs to),
#'                     and \code{group} (group of orthologous to which the gene belongs to). More columns can be added as metadata
#'                     for each gene.
#'                  }
#'                     \item{\bold{\code{org_meta}}: (optional) A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'                     containging additional metadata for organisms. This \code{data.frame} must have a column named "org" with
#'                     valid organisms names (that is, they should match with those provided in \code{DF}, column \code{org}), and
#'                     additional columns will be used as metadata. Each row should correspond to each organism.
#'
#'                  }
#'                     \item{\bold{\code{group_meta}}: (optional) A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'                     containging additional metadata for clusters. This \code{data.frame} must have a column named "group" with
#'                     valid organisms names (that is, they should match with those provided in \code{DF}, column \code{group}), and
#'                     additional columns will be used as metadata. Each row should correspond to each cluster.
#'
#'                  }
#'                     \item{\bold{\code{sep}}: A separator. By default is '__'(two underscores). It will be used to
#'                     create a unique \code{gid} (gene identifier) for each gene. \code{gid}s are created by pasting
#'                     \code{org} to \code{gene}, separated by \code{sep}.
#'                  }
#'                     \item{\bold{\code{verbose}}: \code{logical}. Whether to display progress messages when loading class.
#'
#'                  }
#'
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{An R6 object of class PgR6M. It contains basic fields and methods for analyzing a pangenome. It also
#'                     contains additional statistical methods for analize it, and methods to make basic exploratory plots.}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#' @section Public Methods:
#' \describe{
#'     \item{\code{add_metadata(map = 'org', df)}}{
#'         \itemize{
#'             \item{Add metadata to the object. You can add metadata to each organism, to each
#'             group of orthologous, or to each gene. Elements with missing data should be filled
#'             by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'             data).}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{map}}: \code{character} identifiying the metadata to map. Can
#'                     be one of \code{"org"}, \code{"group"}, or \code{"gid"}.}
#'                     \item{\bold{\code{df}}: \code{data.frame} or \code{DataFrame} with the metadata to
#'                     add. For ethier case, a column named as \code{"map"} must exists, which should
#'                     contain identifiers for each element. In the case of adding gene (\code{gid})
#'                     metadata,each gene should be referenced by the name of the organism and the name
#'                     of the gene as provided in the \code{"DF"} data.frame, separated by the
#'                     \code{"sep"} argument.}
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{\code{self} invisibly, but with additional metadata.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{drop(x)}}{
#'         \itemize{
#'             \item{Drop an organism from the dataset. This method allows to hide an organism
#'             from the real dataset, ignoring it in downstream analyses. All the fields and
#'             methods will behave as it doesn't exist. For instance, if you decide to drop
#'             organism 1, the \code{$pan_matrix} field (see below) would not show it when
#'             called.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{x}}: \code{character} or \code{numeric}. The name of the
#'                     organism wanted to be dropped, or its numeric id as returned in
#'                     \code{$organism} field (see below).
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{\code{self} invisibly, but with \code{x} dropped. It isn't necessary
#'                     to assign the function call to a new object, nor to re-write it as R6 objects
#'                     are mutable.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{recover(x)}}{
#'         \itemize{
#'             \item{Recover a previously \code{$drop()}ped organism (see above). All fields
#'             and methods will start to behave considering this organism again.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{x}}: \code{character} or \code{numeric}. The name of the
#'                     organism wanted to be recover, or its numeric id as returned in
#'                     \code{$dropped} field (see below).
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{\code{self} invisibly, but with \code{x} recovered. It isn't necessary
#'                     to assign the function call to a new object, nor to re-write it as R6 objects
#'                     are mutable.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{rarefact()}}{
#'         \itemize{
#'             \item{Rarefact pangenome or corgenome. Compute the number of genes which belong to
#'             the pangenome or to the coregenome, for a number of random permutations of
#'             increasingly bigger sample of genomes.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{what}}: One of \code{"pangenome"} or \code{"coregenome"}.}
#'                     \item{\bold{\code{n.perm}}: The number of permutations to compute}
#'                 }
#'              }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A \code{matrix}, rows are the number of genomes added, columns are
#'                     permutations, and the cell number is the number of genes in ethier category.
#'                     }
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{dist(method = 'bray', binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE, ...)}}{
#'         \itemize{
#'             \item{Compute distance between all pairs of genomes. The default dist method is
#'             \code{"bray"} (Bray-Curtis distance). Annother used distance method is \code{"jaccard"},
#'             but you should set \code{binary = FALSE} (see below) to obtain a meaningful result.
#'             See \code{\link[vegan]{vegdist}} for details, this is just a wrapper function.
#'             }
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{method}}: The distance method to use. See \link[vegan]{vegdist}
#'                     for available methods, and details for each one.
#'                     }
#'                     \item{\bold{\code{binary}}: Transform abundance matrix into a presence/absence
#'                     matrix before computing distance.}
#'                     \item{\bold{\code{diag}}: Compute diagonals.}
#'                     \item{\bold{\code{upper}}: Return only the upper diagonal.}
#'                     \item{\bold{\code{na.rm}}: Pairwise deletion of missing observations when
#'                     computing dissimilarities.}
#'                     \item{\bold{\code{...}}: Other parameters. See \link[vegan]{vegdist} for details.}
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A \code{dist} object containing all pairwise dissimilarities between genomes.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{pan_pca(center = TRUE, scale. = FALSE, ...)}}{
#'         \itemize{
#'             \item{Performs a principal components analysis on the panmatrix}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{center}}: a logical value indicating whether the variables should be shifted
#'                     to be zero centered. Alternately, a vector of length equal the number of columns of x can be
#'                     supplied. The value is passed to scale.
#'                     }
#'                     \item{\bold{\code{scale.}}: a logical value indicating whether the variables should be scaled
#'                     to have unit variance before the analysis takes place. The default is TRUE.
#'                     }
#'                     \item{\bold{\code{...}}: Other arguments. See \link[stats]{prcomp}.
#'                     }
#'                  }
#'
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{Returns a list with class "prcomp". See \link[stats]{prcomp} for more information.
#'                     }
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{pg_power_law_fit(raref, ...)}}{
#'         \itemize{
#'             \item{Fits a power law curve for the pangenome rarefaction simulation.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{raref}}: (Optional) A rarefaction matrix, as returned by \code{rarefact()}.
#'                     }
#'                     \item{\bold{\code{...}}: Further arguments to be passed to \code{rarefact()}. If \code{raref}
#'                     is missing, it will be computed with default arguments, or with the ones provided here.
#'                     }
#'                  }
#'
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
#'                     with fitted intercept and decay parameters. An attribute \code{"alpha"} is also returned (If
#'                     \code{alpha>1}, then the pangenome is closed, otherwise is open.)
#'                     }
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{cg_exp_decay_fit(raref, pcounts = 10, ...)}}{
#'         \itemize{
#'             \item{Fits an exponential decay curve for the coregenome rarefaction simulation.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{raref}}: (Optional) A rarefaction matrix, as returned by \code{rarefact()}.
#'                     }
#'                     \item{\bold{\code{pcounts}}: An integer of pseudo-counts. This is used to better fit the function
#'                     at small numbers, as the linearization method requires to substract a constant C, which is the
#'                     coregenome size, from \code{y}. As \code{y} becomes closer to the coregenome size, this operation
#'                     tends to 0, and its logarithm goes crazy. By default \code{pcounts=10}.
#'                     }
#'                     \item{\bold{\code{...}}: Further arguments to be passed to \code{rarefact()}. If \code{raref}
#'                     is missing, it will be computed with default arguments, or with the ones provided here.
#'                     }
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{fluidity(n.sim = 10)}}{
#'         \itemize{
#'             \item{Computes the genomic fluidity, which is a measure of population
#'             diversity. See \code{\link[micropan]{fluidity}} for more details.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{n.sim}}: An integer specifying the number of random samples
#'                     to use in the computations.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A list with two elements, the mean fluidity and its sample standard
#'                     deviation over the n.sim computed values.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{binomix_estimate(K.range = 3:5, core.detect.prob = 1, verbose = TRUE)}}{
#'         \itemize{
#'             \item{Fits binomial mixture models to the data given as a pan-matrix. From the
#'             fitted models both estimates of pan-genome size and core-genome size are
#'             available. See \code{\link[micropan]{binomixEstimate}} for more details.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{K.range}}: The range of model complexities to explore. The
#'                     vector of integers specify the number of binomial densities to combine in the
#'                     mixture models.
#'                  }
#'                     \item{\bold{\code{core.detect.prob}}: The detection probability of core genes.
#'                     This should almost always be 1.0, since a core gene is by definition always
#'                     present in all genomes, but can be set fractionally smaller.
#'                  }
#'                     \item{\bold{\code{verbose}}: Logical indicating if textual output should be
#'                     given to monitor the progress of the computations.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A \code{Binomix} object (\code{micropan} package), which is a small (S3)
#'                     extension of a \code{list} with two components. These two components are named
#'                     \code{BIC.table} and \code{Mix.list}. Refer to the \code{micropan} function
#'                     \code{\link[micropan]{binomixEstimate}} for a detailed explanation of this
#'                     method.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{gg_barplot()}}{
#'         \itemize{
#'             \item{Plot a barplot with the frequency of genes within the total number of
#'             genomes.}
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A barplot, and a \code{gg} object (\code{ggplot2} package) invisibly.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{gg_binmap()}}{
#'         \itemize{
#'             \item{Plot a pangenome binary map representing the presence/absence of each
#'             gene within each organism.}
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A binary map (\code{ggplot2::geom_raster()}), and a \code{gg} object (\code{ggplot2}
#'                     package) invisibly.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{gg_dist(dist = "Jaccard", ...)}}{
#'         \itemize{
#'             \item{Plot a heatmap showing the computed distance between all pairs of organisms.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{dist}}: Distance method. One of "Jaccard" (default), or "Manhattan",
#'                     see above.
#'                  }
#'                     \item{\bold{\code{...}}: More arguments to be passed to \code{\link[micropan]{distManhattan}}.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A heatmap (\code{ggplot2::geom_tile()}), and a \code{gg} object (\code{ggplot2}
#'                     package) invisibly.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{gg_pca()}}{
#'         \itemize{
#'             \item{Plot a scatter plot of a Principal Components Analysis.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{colour}}: The name of the column in \code{$organisms} field from which points will take
#'                     colour (if provided). \code{NULL} (default) renders black points.
#'                     }
#'                     \item{\bold{...}}: Arguments to be passed to \code{$pan_pca(), or to \link[ggplot2]{autoplot} generic for
#'                     class \code{prcomp}.}
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A scatter plot (\code{ggplot2::autoplot()}), and a \code{gg} object (\code{ggplot2}
#'                     package) invisibly.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{gg_pie()}}{
#'         \itemize{
#'             \item{Plot a pie chart showing the number of clusters of each pangenome category: core,
#'             shell, or cloud.}
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A pie chart (\code{ggplot2::geom_bar() + coord_polar()}), and a \code{gg} object
#'                     (\code{ggplot2} package) invisibly.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{gg_curves()}}{
#'         \itemize{
#'             \item{Plot pangenome and/or coregenome curves with the fitted functions returned by \code{pg_power_law_fit()}
#'             and \code{cg_exp_decay_fit()}. You can add points by adding \code{+ geom_points()}, of ggplot2 package}
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A scatter plot, and a \code{gg} object (\code{ggplot2} package) invisibly.}
#'                 }
#'             }
#'          }
#'     }
#' }
#'
#'
#' @section Public Fields:
#' \describe{
#'     \item{\bold{\code{$pan_matrix}}}{: The panmatrix. Rows are organisms, and
#'     columns are groups of orthologous. Cells indicates the presence (>=1) or
#'     absence (0) of a given gene, in a given organism. Cells can have values
#'     greater than 1 if contain in-paralogs.}
#'     \item{\bold{\code{$organisms}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with available
#'     organism names, and organism number identifier as \code{rownames()}. (Dropped
#'     organisms will not be displayed in this field, see \code{$dropped} below).
#'     Additional metadata will be shown if provided, as additional columns.}
#'     \item{\bold{\code{$clusters}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with the groups
#'     of orthologous (clusters). Additional metadata will be shown as additional columns,
#'     if provided before. Each row corresponds to each cluster.}
#'     \item{\bold{\code{$genes}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}link[IRanges:DataFrameList-class]{SplitDataFrameList}} object with
#'     one entry per cluster. Each element contains a \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'     with gene ids (\code{<gid>}) and additional metadata, if provided. \code{gid} are
#'     created by \code{paste}ing organism and gene names, so duplication in gene names
#'     are avoided.}
#'     \item{\bold{\code{$core_level}}}{: The percentage of organisms a gene must be in
#'     to be considered as part of the coregenome. \code{core_level = 95} by default.
#'     Can't be set above 100, and below 85 raises a warning.}
#'     \item{\bold{\code{$core_genes}}}{: Like \code{genes}, but only showing core genes.}
#'     \item{\bold{\code{$core_clusters}}}{: Like \code{$clusters}, but only showing core
#'     clusters.}
#'     \item{\bold{\code{$cloud_genes}}}{: Like \code{genes}, but only showing cloud genes.
#'     These are defined as those clusters which contain a single gene (singletons), plus
#'     those which have more than one but its organisms are probably clonal due to identical
#'     general gene content. Colloquially defined as strain-specific genes.}
#'     \item{\bold{\code{$cloud_clusters}}}{: Like \code{$clusters}, but only showing cloud
#'     clusters as defined above.}
#'     \item{\bold{\code{$shell_genes}}}{: Like \code{genes}, but only showing shell genes.
#'     These are defined as those clusters than don't belong nethier to the core genome,
#'     nor to cloud genome. Colloquially defined as genes that are present in some but not
#'     all strains, and that aren't strain-specific.}
#'     \item{\bold{\code{$shell_clusters}}}{: Like \code{$clusters}, but only showing shell
#'     clusters, as defined above.}
#'     \item{\bold{\code{$summary_stats}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with
#'     information about the number of core, shell, and cloud clusters, as well as the
#'     total number of clusters.}
#'     \item{\bold{\code{$random_seed}}}{: The last \code{.Random.seed}. Used for
#'     reproducibility purposes only.}
#'     \item{\bold{\code{$dropped}}}{: A \code{character} vector with dropped organism
#'     names, and organism number identifier as \code{names()}}
#' }
#'
#'
#' @section Special Methods:
#' \describe{
#'     \item{\code{clone(deep = FALSE)}}{
#'         \itemize{
#'             \item{Method for copying an object. See \href{https://adv-r.hadley.nz/r6.html#r6-semantics}{emph{Advanced R}} for the intricacies of R6 reference semantics.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{deep}}: logical. Whether to recursively clone nested R6 objects.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{Cloned object of this class.}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#' @importFrom R6 R6Class
#' @importFrom micropan fluidity binomixEstimate
#' @import ggplot2
#' @import ggfortify
#' @import shiny
#' @import shinyjs
#' @import shinyWidgets
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom plotly plot_ly plotlyOutput renderPlotly layout add_lines
#' @importFrom heatmaply heatmaply
#' @importFrom reshape2 melt
#' @importFrom vegan vegdist
#' @importFrom dendextend seriate_dendrogram
#' @export
PgR6M <- R6Class('PgR6M',

                 inherit = PgR6,

                 public = list(

                   initialize = function(DF,
                                         org_meta,
                                         group_meta,
                                         sep = '__',
                                         verbose = TRUE){

                     super$initialize(DF = DF,
                                      org_meta,
                                      group_meta,
                                      sep = sep,
                                      verbose = TRUE)

                   },

                   #Analyses Methods

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

                   pan_pca = function(center = TRUE, scale. = FALSE, ...){
                     prcomp(self$pan_matrix, center = center, scale. = scale., ...)
                   },

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

                   fluidity = function(n.sim = 10){
                     #micropan::fluidity()
                     fluidity(self$pan_matrix, n.sim = n.sim)
                   },

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
                   gg_barplot = function(){
                     pm <- self$pan_matrix
                     pm[which(pm > 0, arr.ind = TRUE)] <- 1L
                     dd <- as.data.frame(table(colSums(pm)))
                     ggplot(dd, aes(x=Var1, y=Freq)) +
                       geom_bar(stat='identity') +
                       ylab('Number of genes') +
                       xlab('Number of genomes')
                   },

                   gg_binmap = function(){
                     tpm <- t(self$pan_matrix)
                     tpm[which(tpm > 0, arr.ind = TRUE)] <- 1L
                     bm <- as.data.frame(tpm)
                     or <- order(rowSums(bm), decreasing = TRUE)
                     lvls <- rownames(bm)[or]
                     bm$gene <- factor(rownames(bm), levels = lvls)
                     bm <- melt(bm, 'gene')
                     bm$value <- factor(bm$value, levels = c(1, 0))
                     ggplot(bm, aes(gene, variable, fill=value)) +
                       geom_raster() +
                       theme(axis.ticks.x = element_blank(),
                             axis.text.x = element_blank()) +
                       scale_fill_grey(start = .2, end = .9)
                   },

                   gg_dist = function(method = 'bray', ...){
                     m <- as.matrix(self$dist(method = method, ...))
                     me <- melt(m)
                     ggplot(me, aes(Var1, Var2, fill=value)) +
                       geom_tile()
                   },

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

                   gg_pie = function(){
                     st <- self$summary_stats[-1, ]
                     # st <- st[-1, ]
                     st$Category <- factor(st$Category,
                                           levels = c("Cloud", "Shell", "Core"))
                     ggplot(as.data.frame(st), aes(x='', y=Number, fill=Category)) +
                       geom_bar(width = 1, stat = "identity") +
                       coord_polar("y", start=0)
                   },

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
                       lrarm[[x]]$category <- x
                       lrarm[[x]]
                      })
                     df <- Reduce(rbind, ll)

                     #plot
                     g <- ggplot(df, aes(x=factor(Var1), y=value, colour=category)) +
                       xlab('Number of genomes') +
                       ylab('Number of clusters')
                     for (i in seq_along(what)){
                       g <- g +
                         stat_function(data = df[which(df$category == what[i]), ],
                                       fun = lfun[[what[i]]]$formula)
                     }
                     g
                   },


                   runShinyApp = function(){

                     pg <- self

                     app <- list(ui = fluidPage(

                       useShinyjs(),
                       # setBackgroundColor(),

                       titlePanel("Pagoo Shiny App"),

                       tabsetPanel(

                         tabPanel(
                           "Summary",

                           fluidRow(
                             column(width = 3, sliderInput("core_level",
                                                           "Core Level",
                                                           min=85,
                                                           max=100,
                                                           value=95)),

                             column(width = 2,
                                    radioButtons("input_type",
                                                 label = "Selection type:",
                                                 choices = c("organism", "metadata"),
                                                 selected = "organism")
                             ),

                             uiOutput("select_organisms")

                           ),

                           fluidRow(
                             column(width = 2,
                                    align = "center",
                                    DT::DTOutput("summary_counts"),
                                    plotlyOutput("core_evol", width = "auto")
                             ),
                             column(width = 3, offset = 0, plotlyOutput("pie")),
                             column(width = 3, offset = 0, plotlyOutput("barplot")),
                             column(width = 4, offset = 0, plotlyOutput("curves"))
                           ),

                           # hr(),

                           fluidRow(
                             column(width = 6,style = "margin-top:-5em;",
                                    align = "center",
                                    DT::DTOutput("core_clusters"), br()
                             ),
                             column(
                               width = 6,style = "margin-top:-5em;",
                               align = "center",
                               DT::DTOutput("core_genes"), br()
                             )
                           )


                         ), # Ends tabPanel "Summary"

                         tabPanel(
                           "Accessory",

                           fluidRow(
                             column(
                               width = 3,
                               uiOutput(outputId = "accs_slider")
                             ),

                             column(
                               width = 1,
                               br(),
                               uiOutput("pca_x"),
                             ),

                             column(
                               width = 1,
                               br(),
                               uiOutput("pca_y")
                             ),

                             column(
                               width = 2,
                               # offset = 1,
                               # br(),
                               checkboxInput("check_color_meta",
                                             label = "Color by",
                                             value = FALSE),
                               uiOutput("pca_color_meta")

                             ),

                             column(
                               width = 4,
                               # offset = 1,
                               # br(),
                               DT::DTOutput("pca_summary")
                             )
                           ),

                           hr(),

                           fluidRow(
                             splitLayout(
                               plotlyOutput("heat_freq", height = "350px"),
                               plotlyOutput("pca", height = "350px"),
                               cellWidths = c("60%", "40%")
                             )
                           ),

                           hr(),

                           fluidRow(
                             column(
                               width = 6,align = "center",
                               DT::DTOutput("accs_clusters"),
                             ),
                             column(
                               width = 6, align = "center",
                               DT::DTOutput("accs_genes")
                             )

                           )

                         ) # Ends tabPanel "Accessory"

                       ) # Ends tabPanel

                     ),


                     server = function(input, output, session){

                       if (dim(pg$organisms)[2] < 2){
                         shinyjs::disable(selector = "[type=radio][value=metadata]")
                         shinyjs::runjs("$('[type=radio][value=metadata]').parent().parent().addClass('disabled').css('opacity', 0.4)")
                         shinyjs::disable(selector = "[type=check][value=check_color_meta]")
                         shinyjs::runjs("$('[type=check][value=check_color_meta]').parent().parent().addClass('disabled').css('opacity', 0.4)")
                       }


                       #############
                       ## SUMMARY ##
                       #############

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

                     })

                     runApp(app)

                   }


                 )

)
