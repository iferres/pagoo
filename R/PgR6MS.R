#PgR6MS

# sequences <- structure(list(org1 = structure(c("jwnxoarfiq", "wajpcygeox", "xvgfkqacpb", "rewbxafzdj"),
#                                              .Names = c("gene1", "gene2", "gene3", "gene4")),
#                             org2 = structure(c("vqphkgajem", "ictkgfpmvh"),
#                                              .Names = c("gene1", "gene2")),
#                             org3 = structure(c("ixwflohqzs", "fuivpjtehs"),
#                                              .Names = c("gene1", "gene2"))),
#                        .Names = c("org1", "org2", "org3"))

##' Expect 1) a named list of named character vector Names of list is names of
##' organisms, names of character vector are gene names; 2) a named list of
##' BStringSet objects (same requirements as (1), but with BStringSet names as
##' gene names); 3) a BStringSetList (same requirements as (3) but
##' BStringSetList names are organisms names).
##'
##'




#' @name PgR6MS
#' @description PgR6 with Methods and Sequences.
#'
#' #'  Inherits: \code{\link[pagoo]{PgR6M}}
#'
#'
#'  @section Class Constructor:
#' \describe{
#'     \item{\code{new(cluster_df, sep = "__", sequences)}}{
#'         \itemize{
#'             \item{~~DESCRIBE THE METHOD~~}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{cluster_df}}: ~~DESCRIBE THIS PARAMETER~~
#'                  }
#'                     \item{\bold{\code{sep}}: ~~DESCRIBE THIS PARAMETER~~
#'                  }
#'                     \item{\bold{\code{sequences}}: ~~DESCRIBE THIS PARAMETER~~
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{~~WHAT DOES THIS RETURN~~}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#'  @section Public Methods:
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
#'     \item{\code{get_core_seqs(max_per_org = 1, fill = TRUE)}}{
#'         \itemize{
#'             \item{~~DESCRIBE THE METHOD~~}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{max_per_org}}: ~~DESCRIBE THIS PARAMETER~~
#'                  }
#'                     \item{\bold{\code{fill}}: ~~DESCRIBE THIS PARAMETER~~
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{~~WHAT DOES THIS RETURN~~}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#'  @section Public Fields:
#' \describe{
#'      \item{\bold{\code{pan_matrix}}}{: The panmatrix. Rows are organisms, and
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
#'     \item{\bold{\code{sequences}}}{: ~~DESCRIBE THIS FIELD~~}
#' }
#'
#'
#'  @section Special Methods:
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
#' @importFrom GenomicRanges mcols mcols<-
#' @importFrom Biostrings BStringSet DNAStringSet
#' @importFrom S4Vectors split
#' @export
PgR6MS <- R6Class('PgR6MS',

                  inherit = PgR6M,

                  private = list(
                    .sequences = NULL
                  ),

                  public = list(

                    initialize = function(cluster_df,
                                          sep = '__',
                                          sequences){

                      super$initialize(cluster_df = cluster_df,
                                       sep = sep)

                      # Check input sequences
                      if (class(sequences)%in%c('list',
                                                'BStringSetList',
                                                'DNAStringSetList')){

                        norgs <- length(sequences)
                        orgNames <- names(sequences)

                        if (length(orgNames)==0){
                          stop('Unnamed list.')
                        }

                        clss <- unique(sapply(sequences, class))
                        genNames <- lapply(sequences, names)

                        if (clss%in%c('character',
                                      'BStringSet',
                                      'DNAStringSet')){

                          if (length(genNames)==0){
                            stop('Unnamed sequences.')
                          }

                          gids <- mapply(paste,
                                         orgNames, genNames,
                                         MoreArgs = list(sep = sep))
                          gids <- unlist(gids, use.names = FALSE)
                          sqs <- unlist(sequences, use.names = FALSE)
                          names(sqs) <- gids
                          private$.sequences <- DNAStringSet(sqs)

                        }else{
                          stop('Unrecognized sequences format.')
                        }

                      }else{
                        stop('Unrecognized sequences format')
                      }

                      dtgid <- private$.dt[, gid]

                      if (!all(dtgid%in%gids)){
                        stop('Missing sequences: some gid do not match any sequence name.')
                      }

                      if (any(!gids%in%dtgid)){
                        warning('Missing gid: some sequence names do not match to any gid. Continuing anyway..\n', immediate. = TRUE)
                      }

                      # Add metadata to sequences
                      spl <- strsplit(gids, sep)
                      mcols(private$.sequences)$org <- vapply(spl, '[', 1,
                                                              FUN.VALUE = NA_character_)
                      mcols(private$.sequences)$group <- vapply(gids,
                                                                function(x){
                                                                  wh <- which(dtgid==x)
                                                                  if (length(wh))
                                                                    as.character(private$.dt$group[wh])
                                                                  else
                                                                    NA_character_
                                                                }, FUN.VALUE = NA_character_)

                      # # Expect a named list (names = organisms names), with a
                      # # named vector of strings. Names of strings are gene
                      # # names. Applies separator prior to check.
                      # nmsq <- names(sequences)
                      # if (!all(self$organisms %in% nmsq))
                      #   stop('organism names in cluster_list do not match names of sequences.')
                      #
                      # # Rename genes applying separator. This is to later make
                      # # them hashable.
                      # sequences <- lapply(seq_along(sequences), function(i){
                      #   sq <- sequences[i]
                      #   norg <- names(sq)
                      #   sqv <- sq[[1]]
                      #   names(sqv) <- paste(norg, sep, names(sqv), sep = '')
                      #   sqv
                      # })
                      # names(sequences) <- nmsq
                      #
                      # # Extract for each organism the gene names in each
                      # # cluster.
                      # namesInClusters <- lapply(self$organisms, function(x){
                      #   unlist(sapply(self$clusters, function(y){
                      #    rr <- y[which(names(y)==x)]
                      #    rr[!is.na(rr)]
                      #   }), use.names = FALSE)
                      # })
                      # names(namesInClusters) <- self$organisms
                      #
                      #
                      # # Check if gene names in clusters match gene names in
                      # # sequences. Stop if not.
                      # genesOK <- all(sapply(self$organisms, function(x){
                      #   all(namesInClusters[[x]] %in% names(sequences[[x]]))
                      # }))
                      # genesOK2 <- all(sapply(self$organisms, function(x){
                      #   all(names(sequences[[x]] %in% namesInClusters[[x]]))
                      # }))
                      # if (!genesOK){
                      #   stop('gene names in sequences are not contained in gene names in cluster_list')
                      # }else if (!genesOK2) {
                      #   warning('Not all gene names in sequences are present in cluster_list. Moving on.')
                      # }
                      #
                      # # Create private env and populate it with sequences.
                      # names(sequences) <- NULL
                      # sequences <- unlist(sequences, use.names = TRUE)
                      # private$.sequences <- BStringSet(sequences, use.names = TRUE)
                      # ptt <- paste0(sep, '\\w+$')
                      # mcols(private$.sequences)$organisms <- sub(ptt, '', names(private$.sequences))
                      # # private$.sequences <- list2env(as.list(sequences),
                      # #                                 hash = TRUE)
                    },

                    # Methods for sequences
                    #' @importFrom S4Vectors mcols
                    #' @importFrom BiocGenerics table
                    get_core_seqs = function(max_per_org = 1, fill = TRUE){

                      if (!class(max_per_org)%in%c('logical', 'NULL', 'numeric', 'integer'))
                        stop('"max_per_org" is not numeric, NULL, or NA')
                      if (!is.logical(fill))
                        stop('"fill" is not logical')

                      ccl <- self$core_clusters
                      orgs <- self$organisms
                      sqs <- private$.sequences
                      mcls <- mcols(sqs)
                      whcore <- which(mcls$group %in% ccl)
                      sqs <- sqs[whcore]
                      mcls <- mcols(sqs)

                      if (!(is.na(max_per_org) | is.null(max_per_org))){
                        wma <- which(self$pan_matrix[, ccl] > max_per_org, arr.ind = T)
                        if (dim(wma)[1]){
                          ogs2fix <- ccl[wma[, 2]]
                          orgs2fix <- rownames(wma)
                          remrw <- apply(cbind(ogs2fix, orgs2fix), 1, function(x){
                            wh <- which(mcls$group==x[1] & mcls$org==x[2])
                            wh[-seq_len(max_per_org)]
                          })
                          torm <- unlist(remrw)
                          if (length(torm)){
                            sqs <- sqs[-torm]
                            mcls <- mcols(sqs)
                          }
                        }
                      }

                      # Fill '0' cells with empthy DNAStringSet
                      # spl <- split(sqs, mcls$group)
                      # nro <- elementNROWS(spl)
                      if (fill){
                        tbl <- table(mcls)
                        tofll <- which(tbl==0, arr.ind = TRUE)
                        if (dim(tofll)[1]){
                          sqs <- append(sqs, DNAStringSet(rep('', dim(tofll)[1])))
                          mcls <- mcols(sqs)
                          mfill <- which(is.na(mcls$org))
                          mcols(sqs)$org[mfill] <- rownames(tofll)
                          mcols(sqs)$group[mfill] <- ccl[tofll[, 2]]
                          mcls <- mcols(sqs)
                        }
                      }

                      ## Order ! ! !
                      sqs <- sqs[order(mcls$org)]
                      mcls <- mcols(sqs)

                      # Return
                      split(sqs, mcls$group)
                    }

                  ),

                  active = list(

                    sequences = function(){
                      dn <- dimnames(self$pan_matrix)
                      ogs <- dn[[2]]
                      orgs <- dn[[1]]
                      sqs <- private$.sequences
                      sset <- which(mcols(sqs)$group %in% ogs &
                                      mcols(sqs)$org %in% orgs)
                      split(sqs[sset], mcols(sqs[sset])$group)
                    }
                  )
)



# # Subset method, matrix-like.
# `[.SequenceList` <- function(x, i, j){
#
#   Narg <- nargs()
#   orgs <- attr(x, 'organisms')
#   sep <- attr(x, 'separator')
#
#   # simple subset
#   if (Narg<3){
#     if (missing(i)) {return(x)} else {rr <- unclass(x)[i]}
#   }
#
#   #double arg
#   #Take all orgs, subset clusters
#   if (missing(i)){
#     rr <- unclass(x)[j]
#   }
#
#   pattern <- paste0('^',orgs, sep)
#   #Take all clusters, subset orgs
#   if (missing(j)){
#     if (is.character(i)){
#       i <- which(orgs%in%i)
#     }
#     sorgs <- orgs[i]
#     rr <- lapply(x, function(y){
#       ss <- y[which(names(y)%in%sorgs)]
#       ss[!is.na(ss)]
#     })
#   }
#
#
#
# }
