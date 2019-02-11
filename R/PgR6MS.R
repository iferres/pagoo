#PgR6MS

# sequences <- structure(list(org1 = structure(c("jwnxoarfiq", "wajpcygeox", "xvgfkqacpb", "rewbxafzdj"),
#                                              .Names = c("gene1", "gene2", "gene3", "gene4")),
#                             org2 = structure(c("vqphkgajem", "ictkgfpmvh"),
#                                              .Names = c("gene1", "gene2")),
#                             org3 = structure(c("ixwflohqzs", "fuivpjtehs"),
#                                              .Names = c("gene1", "gene2"))),
#                        .Names = c("org1", "org2", "org3"))

PgR6MS <- R6Class('PgR6MS',

                  inherit = PgR6M,

                  private = list(
                    p_sequences = NULL
                  ),

                  public = list(

                    initialize = function(cluster_list,
                                          prefix = 'group',
                                          sequences){

                      super$initialize(cluster_list = cluster_list,
                                       prefix = prefix)

                      # Check input sequences
                      # Expect a named list (names = organisms names), with a
                      # named vector of strings. Names of strings are gene
                      # names.
                      if (!all(self$organisms %in% names(sequences)))
                        stop('organism names in cluster_list do not match names of sequences.')

                      namesInClusters <- lapply(self$organisms, function(x){
                        unlist(sapply(self$clusters, function(y){
                         rr <- y[which(names(y)==x)]
                         rr[!is.na(rr)]
                        }), use.names = FALSE)
                      })
                      names(namesInClusters) <- self$organisms

                      genesOK <- all(sapply(self$organisms, function(x){
                        all(namesInClusters[[x]]%in%names(sequences[[x]]))
                      }))
                      genesOK2 <- all(sapply(self$organisms, function(x){
                        all(names(sequences[[x]]%in%namesInClusters[[x]]))
                      }))

                      if (!genesOK){
                        stop('gene names in sequences are not contained in gene names in cluster_list')
                      }else if (!genesOK2) {
                        warning('Not all gene names in sequences are present in cluster_list. Moving on.')
                      }

                      private$p_sequences <- lapply(sequences, function(x){
                        list2env(as.list(x), hash = TRUE)
                      })
                    }

                    # Methods for sequences
                  ),

                  active = list()
)
