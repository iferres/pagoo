# PgR6M

#' @name PgR6M
#' @title PgR6 class with methods.
#' @description PgR6 with Methods.
#' Inherits: \code{\link[pagoo]{PgR6}}
#'
#'
#' @section Class Constructor:
#' \describe{
#'     \item{\code{new(cluster_df, sep = "__")}}{
#'         \itemize{
#'             \item{Create a \code{PgR6} object.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{cluster_df}}: A \code{data.frame} or \code{data.table} containing at least the
#'                     following columns: \code{gene} (gene name), \code{org} (organism name to which the gene belongs to),
#'                     and \code{group} (group of orthologous to which the gene belongs to). More columns are allowed but
#'                     this basic class do not contain any methods to handle them.
#'                  }
#'                     \item{\bold{\code{sep}}: A separator. By default is '__'(two underscores). It will be used to
#'                     create a unique \code{gid} (gene identifier) for each gene. \code{gid}s are created by pasting
#'                     \code{org} to \code{gene}, separated by \code{sep}.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{An R6 object of class PgR6. It contains basic fields and methods for
#'                     analyzing a pangenome.}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#' @section Public Methods:
#' \describe{
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
#'     \item{\code{dist_jaccard()}}{
#'         \itemize{
#'             \item{Computes Jaccard's distance between all pairs of genomes. See
#'             \code{\link[micropan]{distJaccard}} for details, this is just a wrapper
#'             function.}
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A \code{dist} object (see \link[stats]{dist}) containing all pairwise
#'                     Jaccard distances between genomes.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{dist_manhattan(scale = 0, weights = rep(1, dim(self$pan_matrix)[2]))}}{
#'         \itemize{
#'             \item{Computes the (weighted) Manhattan distances beween all pairs of genomes.
#'             See \code{\link[micropan]{distManhattan}} for details. }
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{scale}}: An optional scale to control how copy numbers should
#'                     affect the distances.
#'                  }
#'                     \item{\bold{\code{weights}}: Vector of optional weights of gene clusters.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A \code{dist} object (see \link[stats]{dist}) containing all pairwise Manhattan
#'                     distances between genomes.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{heaps(n.perm = 100)}}{
#'         \itemize{
#'             \item{Estimating if a pan-genome is open or closed based on a Heaps law model.
#'             See \code{\link[micropan]{heaps}} for more details.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{n.perm}}: The number of random permutations of genome
#'                     ordering.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A vector of two estimated parameters: The Intercept and the decay
#'                     parameter alpha. If alpha<1.0 the pan-genome is open, if alpha>1.0 it is closed.}
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
#'     \item{\code{gg_dendro(dist_method = "Jaccard", hclust_method = "complete", ...)}}{
#'         \itemize{
#'             \item{Plot a dendrogram showing the clustering between organisms.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{dist_method}}: The \code{dist} method to use. Default:
#'                     \code{"Jaccard"}.
#'                  }
#'                     \item{\bold{\code{hclust_method}}: The \code{hclust} method to use. Default:
#'                     \code{"complete"}.
#'                  }
#'                     \item{\bold{\code{...}}: More arguments to be passed to
#'                     \code{\link[micropan]{distManhattan}}, or to \code{\link[ggdendro]{ggdendrogram}}.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{A dendrogram plot (\code{ggdendro::ggdendrogram()}), and a \code{gg} object
#'                     (\code{ggplot2} package) invisibly.}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#' @section Public Fields:
#' \describe{
#'     \item{\bold{\code{pan_matrix}}}{: The panmatrix. Rows are organisms, and
#'     columns are groups of orthologous. Cells indicates the presence (>=1) or
#'     absence (0) of a given gene, in a given organism. Cells can have values
#'     greater than 1 if contain in-paralogs.}
#'     \item{\bold{\code{organisms}}}{: A \code{character} vector with available
#'     organism names, and organism number identifier as \code{names()}. (Dropped
#'     organisms will not be displayed in this field, see \code{$dropped} below).}
#'     \item{\bold{\code{clusters}}}{: A named \code{list} of clusters. Clusters are
#'     shown as \code{data.table} objects containing 3 columns: \code{gid}, \code{org},
#'     and \code{gene}, as explained before.}
#'     \item{\bold{\code{core_level}}}{: The percentage of organisms a gene must be in
#'     to be considered as part of the coregenome. \code{core_level = 95} by default.
#'     Can't be set above 100, and below 85 raises a warning.}
#'     \item{\bold{\code{core_clusters}}}{: A \code{character} vector with core
#'     cluster names, as defined by \code{$core_level}.}
#'     \item{\bold{\code{cloud_clusters}}}{: A \code{character} vector with
#'     cloud clusters. These are defined as those clusters which contain a single
#'     gene (singletons), plus those which have more than one but its organisms are
#'     probably clonal due to identical general gene content. Colloquially defined as
#'     strain-specific genes.}
#'     \item{\bold{\code{shell_clusters}}}{: A \code{character} vector with shell
#'     clusters. These are defined as those clusters than don't belong nethier to the
#'     core genome, nor to cloud genome. Colloquially defined as genes that are
#'     present in some but not all strains, and that aren't strain-specific.}
#'     \item{\bold{\code{summary_stats}}}{: A \code{data.frame} with information about
#'     the number of core, shell, and cloud clusters, as well as the total number of
#'     clusters.}
#'     \item{\bold{\code{random_seed}}}{: The last \code{.Random.seed}. Used for
#'     reproducibility purposes only.}
#'     \item{\bold{\code{dropped}}}{: A \code{character} vector with dropped organism
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
#' @importFrom micropan distJaccard fluidity binomixEstimate
#' @importFrom ggplot2 ggplot aes geom_bar geom_raster geom_tile theme element_blank scale_fill_grey xlab ylab coord_polar
#' @importFrom reshape2 melt
#' @importFrom vegan vegdist
#' @export
PgR6M <- R6Class('PgR6M',

                 inherit = PgR6,

                 public = list(

                   initialize = function(cluster_df,
                                         sep = '__'){

                     super$initialize(cluster_df = cluster_df,
                                      sep = sep)

                   },

                   #Analyses Methods

                   rarefact = function(what = 'pangenome', n.perm = 10){
                     what <- match.arg(what, c('pangenome', 'coregenome'))
                     pm <- self$pan_matrix
                     pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                     norgs <- length(self$organisms)
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

                     method <- match.arg(method, c("manhattan",
                                                   "euclidean",
                                                   "canberra",
                                                   "bray",
                                                   "kulczynski",
                                                   "jaccard",
                                                   "gower",
                                                   "altGower",
                                                   "morisita",
                                                   "horn",
                                                   "mountford",
                                                   "raup" ,
                                                   "binomial",
                                                   "chao",
                                                   "cao",
                                                   "mahalanobis"))

                     if (method == 'jaccard' & binary = FALSE){
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

                   power_law_fit = function(raref, ...){
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

                   exp_decay_fit = function(raref, pcounts = 10, ...){
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

                   gg_dist = function(dist){
                     if (missing(dist)){
                       m <- as.matrix(self$dist_jaccard())
                     }else{
                       cls <- class(dist)
                       if (cls=='character'){
                         dist <- match.arg(dist, c('euclidean',
                                                   'maximum',
                                                   'manhattan',
                                                   'canberra',
                                                   'binary',
                                                   'minkowski'))
                         m <- as.matrix(dist(self$pan_matrix, method = dist))
                       }else if (cls=='function'){
                         m <- as.matrix(dist(self$pan_matrix))
                       }else{
                         stop('Invalid dist argument.')
                       }
                     }
                     me <- melt(m)
                     ggplot(me, aes(Var1, Var2, fill=value)) +
                       geom_tile()
                   },

                   gg_pie = function(){
                     st <- self$summary_stats
                     st <- st[-1, ]
                     st$Category <- factor(st$Category,
                                           levels = c("Cloud", "Shell", "Core"))
                     ggplot(st, aes(x='', y=Number, fill=Category)) +
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
                         self$power_law_fit(raref = lrar[[x]])#$formula
                       }else{
                         self$exp_decay_fit(raref = lrar[[x]])#$formula
                       }
                     })
                     lrarm <- lapply(lrar, melt)
                     ll <- lapply(what, function(x) {
                       olst <- list(data = NULL, formula = NULL)
                       lrarm[[x]]$category <- x
                       lrarm[[x]]
                      })

                     df <- Reduce(rbind, ll)

                     # eq <- lapply(lfun, function(x){
                     #    exp <- paste0('y == ', format(x$formula)[2])
                     #    params <- x$params
                     #    for (i in seq_along(params)){
                     #      exp <- sub(names(params)[i],
                     #                 format(params[i], digits = 3),
                     #                 exp,
                     #                 fixed = TRUE)
                     #    }
                     #    exp
                     # })

                     #plot
                     g <- ggplot(df, aes(x=factor(Var1), y=value, colour=category)) +
                       xlab('Number of genomes') + ylab('Number of clusters')
                     for (i in seq_along(what)){
                       g <- g +
                         stat_function(data = df[which(df$category == what[i]), ],
                                       fun = lfun[[what[i]]]$formula) # +
                         # annotate('text', x = 3, y=4000, label=eq[[what[i]]], parse=T)
                     }
                     g
                   }


                 )

)
