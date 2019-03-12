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
