#' @name PgR6MS
#' @title PgR6 class with Methods and Sequences. Final users should use \code{\link{pagoo}}
#' instead of this, since is more easy to understand.
#' @description PgR6 with Methods and Sequences.
#'  Inherits: \code{\link[pagoo]{PgR6M}}
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

                    #' @description
                    #' Create a \code{PgR6MS} object.
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
                    #' @param DF Deprecated. Use \code{data} instead.
                    #' @param group_meta Deprecated. Use \code{cluster_meta} instead.
                    #' @param sequences Can accept: 1) a named \code{list} of named \code{character} vector. Name of list
                    #' are names of organisms, names of character vector are gene names; or 2) a named \code{list} of \code{\link[Biostrings:XStringSet-class]{DNAStringSetList}}
                    #' objects (same requirements as (1), but with BStringSet names as gene names); or 3) a \code{\link[Biostrings:XStringSetList-class]{DNAStringSetList}}
                    #' (same requirements as (2) but \code{DNAStringSetList} names are organisms names).
                    #' @param verbose \code{logical}. Whether to display progress messages when loading class.
                    #' @return An R6 object of class PgR6MS. It contains basic fields and methods for analyzing a pangenome. It also
                    #' contains additional statistical methods for analize it, methods to make basic exploratory plots, and methods
                    #' for sequence manipulation.
                    initialize = function(data,
                                          org_meta,
                                          cluster_meta,
                                          core_level = 95,
                                          sep = '__',
                                          DF,
                                          group_meta,
                                          sequences,
                                          verbose = TRUE){

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

                      # Check input sequences
                      if (verbose) message('Checking input sequences.')
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
                          if (clss=='character'){
                            sqs <- unlist(sequences, use.names = FALSE)
                            sqs <- DNAStringSet(sqs)
                            names(sqs) <- gids
                          }else{
                            sqs <- DNAStringSetList(sequences)
                            sqs <- unlist(sqs, use.names = FALSE)
                            names(sqs) <- gids
                          }
                          private$.sequences <- sqs

                        }else{
                          stop('Unrecognized sequences format.')
                        }

                      }else{
                        stop('Unrecognized sequences format')
                      }

                      if (verbose) message('Checking that sequence names matches with DataFrame.')

                      dfgid <- private$.data[, 'gid']

                      if (!all(dfgid%in%gids)){
                        stop('Missing sequences: some gid do not match any sequence name.')
                      }

                      if (any(!gids%in%dfgid)){
                        warning('Missing gid: some sequence names do not match to any gid. Continuing anyway..\n', immediate. = TRUE)
                      }

                      # Add metadata to sequences
                      if (verbose) message('Adding metadata to sequences.')
                      spl <- strsplit(gids, sep)
                      mcols(private$.sequences)$org <- vapply(spl, '[', 1,
                                                              FUN.VALUE = NA_character_)
                      mcols(private$.sequences)$cluster <- as.character(private$.data$cluster[match(gids, dfgid)])
                    },

                    # Methods for sequences
                    #' @description
                    #' A field for obtaining core gene sequences is available (see below), but for creating a phylogeny with this
                    #' sets is useful to: 1) have the possibility of extracting just one sequence of each organism on each cluster, in
                    #' case paralogues are present, and 2) filling gaps with empthy sequences in case the core_level was set below 100\%,
                    #' allowing more genes (some not in 100\% of organisms) to be incorporated to the phylogeny. That is the purpose of
                    #' this special function.
                    #' @param max_per_org Maximum number of sequences of each organism to be taken from each cluster.
                    #' @param fill \code{logical}. If fill \code{DNAStringSet} with emphty \code{DNAString} in cases where
                    #' \code{core_level} is set below 100\%, and some clusters with missing organisms are also considered.
                    #' @return A \code{DNAStringSetList} with core genes. Order of organisms on each cluster is conserved, so it is easier
                    #' to concatenate them into a super-gene suitable for phylogenetic inference.
                    #' @importFrom S4Vectors mcols
                    #' @importFrom BiocGenerics table
                    core_seqs_4_phylo = function(max_per_org = 1, fill = TRUE){

                      if (!class(max_per_org)%in%c('logical', 'NULL', 'numeric', 'integer'))
                        stop('"max_per_org" is not numeric, NULL, or NA')
                      if (!is.logical(fill))
                        stop('"fill" is not logical')

                      sqs <- unlist(self$core_sequences, use.names = FALSE)
                      mcls <- mcols(sqs)
                      ccl <- unique(mcls$cluster)

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
                      split(sqs, mcls$cluster)
                    }

                  ),

                  active = list(

                    #' @field sequences A \code{\link[Biostrings:XStringSetList-class]{DNAStringSetList}} with the
                    #' set of sequences grouped by cluster. Each group is accessible as were a list. All
                    #' \code{Biostrings} methods are available.
                    sequences = function(){
                      dn <- dimnames(self$pan_matrix)
                      ogs <- dn[[2]]
                      orgs <- dn[[1]]
                      sqs <- private$.sequences
                      sset <- which(mcols(sqs)$cluster %in% ogs &
                                      mcols(sqs)$org %in% orgs)
                      split(sqs[sset], mcols(sqs[sset])$cluster)
                    },

                    #' @field core_sequences Like \code{$sequences}, but only showing core
                    #' sequences.
                    core_sequences = function(){
                      orgs <- self$organisms$org
                      ogs <- self$core_clusters$cluster
                      sqs <- private$.sequences
                      sset <- which(mcols(sqs)$cluster %in% ogs &
                                      mcols(sqs)$org %in% orgs)
                      split(sqs[sset], mcols(sqs[sset])$cluster)

                    },

                    #' @field cloud_sequences Like \code{$sequences}, but only showing cloud
                    #' sequences as defined above.
                    cloud_sequences = function(){
                      orgs <- self$organisms$org
                      ogs <- self$cloud_clusters$cluster
                      sqs <- private$.sequences
                      sset <- which(mcols(sqs)$cluster %in% ogs &
                                      mcols(sqs)$org %in% orgs)
                      split(sqs[sset], mcols(sqs[sset])$cluster)
                    },

                    #' @field shell_sequences Like \code{$sequences}, but only showing shell
                    #' sequences, as defined above.
                    shell_sequences = function(){
                      orgs <- self$organisms$org
                      ogs <- self$shell_clusters$cluster
                      sqs <- private$.sequences
                      sset <- which(mcols(sqs)$cluster %in% ogs &
                                      mcols(sqs)$org %in% orgs)
                      split(sqs[sset], mcols(sqs[sset])$cluster)
                    }

                  )
)

