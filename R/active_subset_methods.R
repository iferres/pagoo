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
    df <- private$.x$.__enclos_env__$private$.DF
    orgs <- private$.x$organisms
    wh <- which(df$org%in%ii & df$group%in%jj)
    un <- df[wh, 'org', drop = TRUE]
    orgs[orgs$org %in% un, ]
  },

  #' @importFrom S4Vectors split
  genes = function(){
    jj <- private$.j
    ii <- private$.i
    df <- private$.x$.__enclos_env__$private$.DF
    wh <- which(df$org%in%ii & df$group%in%jj)
    df <- df[wh, ]
    split(df, df$group, drop = TRUE)
  },

  clusters = function(){
    jj <- private$.j
    ii <- private$.i
    df <- private$.x$.__enclos_env__$private$.DF
    clu <- private$.x$clusters
    wh <- which(df$org%in%ii & df$group%in%jj)
    un <- df[wh, 'group', drop = TRUE]
    clu[clu$group %in% un, ]
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
    sset <- which(mcols(sqs)$group %in% jj &
                    mcols(sqs)$org %in% ii)
    split(sqs[sset], mcols(sqs[sset])$group)
  }

)
