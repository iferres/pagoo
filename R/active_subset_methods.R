.active_subset <- list(

  pan_matrix = function(){
    ii <- as.integer(names(private$.i))
    jj <- as.integer(names(private$.j))
    private$
      .x$
      .__enclos_env__$
      private$
      .panmatrix[ii, jj, drop=FALSE]
  },

  organisms = function(){
    ii <- private$.i
    jj <- private$.j
    df <- private$.x$.__enclos_env__$private$.data
    orgs <- private$.x$organisms
    wh <- which(df$org%in%ii & df$cluster%in%jj)
    un <- df[wh, 'org', drop = TRUE]
    orgs[orgs$org %in% un, ]
  },

  #' @importFrom S4Vectors split
  genes = function(){
    jj <- private$.j
    ii <- private$.i
    df <- private$.x$.__enclos_env__$private$.data
    wh <- which(df$org%in%ii & df$cluster%in%jj)
    df <- df[wh, ]
    split(df, df$cluster, drop = TRUE)
  },

  clusters = function(){
    jj <- private$.j
    ii <- private$.i
    df <- private$.x$.__enclos_env__$private$.data
    clu <- private$.x$clusters
    wh <- which(df$org%in%ii & df$cluster%in%jj)
    un <- df[wh, 'cluster', drop = TRUE]
    clu[clu$cluster %in% un, ]
  },

  #' @importFrom S4Vectors split
  sequences = function(){
    jj <- private$.j
    ii <- private$.i
    sqs <- private$
      .x$
      .__enclos_env__$
      private$
      .sequences
    sset <- which(mcols(sqs)$cluster %in% jj &
                    mcols(sqs)$org %in% ii)
    split(sqs[sset], mcols(sqs[sset])$cluster)
  }

)
