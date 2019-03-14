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
    # rn <- rownames(self$pan_matrix)
    # private$.x$organisms[private$.x$organisms %in% rn]
    ii <- private$.i
    jj <- private$.j
    fr <- private$
      .x$
      .__enclos_env__$
      private$
      .dt[group%in%jj & org%in%ii, org, ]
    un <- unique(as.character(fr))
    private$.x$organisms[private$.x$organisms %in% un]

  },

  clusters = function(){

    jj <- private$.j
    ii <- private$.i
    split(private$
            .x$
            .__enclos_env__$
            private$
            .dt[group%in%jj & org%in%ii,,],
          by='group',
          keep.by=F,
          drop = TRUE)
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
