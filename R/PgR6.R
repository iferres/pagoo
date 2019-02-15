# Basic pgr6 class
# This is a basic pgr6 class constructor from which inherit more complex ones.

# Test: bare input
# cluster_list <- list(structure(c("gene1", "gene1", "gene1"),
#                                .Names = c("org1","org2", "org3")),
#                      structure(c("gene2", "gene3", "gene2"),
#                                .Names = c("org1","org1", "org3")),
#                      structure(c("gene4", "gene2"),
#                                .Names = c("org1","org2")))
#' @name PgR6
#' @title PgR6 basic class
#' @description A basic \code{pgr6} class constructor.
#' @importFrom R6 R6Class
#' @export
PgR6 <- R6Class('PgR6',

                # Private fields #
                private = list(
                  version = NULL,
                  p_clusters = NULL,
                  p_panmatrix = NULL,
                  p_organisms = NULL,
                  p_dropped = NULL,
                  p_level = NULL
                ),

                # Public functions #
                public = list(

                  initialize = function(cluster_list,
                                        prefix = 'group',
                                        sep = '__'){

                    # Check cluster_list input #
                    # 1. Is list. length > 0?
                    ln_cluster_list <- length(cluster_list)
                    if (any(c(!is.list(cluster_list), ln_cluster_list<=1)))
                      stop('please, check input. Must provide a list of length > 1')
                    # 2. Check cluster_list names. Create if !exist.
                    if (!length(names(cluster_list))){
                      nch <- nchar(ln_cluster_list)
                      np <- paste0('%0', nch, 'd')
                      names(cluster_list) <- paste0(prefix, sprintf(np, 1:ln_cluster_list))
                    }
                    # 3. Check if organism names are given.
                    hasOrgNames <- sapply(cluster_list, function(x) all(length(names(x))) )
                    if (!all(hasOrgNames))
                      stop('you must provide organism names for all genes')
                    # 4. Modify gene names
                    cluster_list <- lapply(cluster_list, function(x){
                      nn <- names(x)
                      ng <- paste(nn, sep, x, sep = '')
                      names(ng) <- nn
                      ng
                    })


                    # Create panmatrix #
                    # 1. Which organisms are?
                    organisms <- unique(unlist(lapply(cluster_list, names)))
                    names(organisms) <- 1:length(organisms)
                    # 2. Compute panmatrix
                    panmatrix <- sapply(cluster_list, function(x){
                      table(factor(names(x), levels = organisms))
                    }, simplify = TRUE)

                    # Populate private$ #
                    private$version <- packageVersion('pgr6')
                    private$p_clusters <- list2env(cluster_list, hash = TRUE)
                    private$p_organisms <- organisms
                    private$p_panmatrix <- panmatrix
                    private$p_level <- 95 #default
                  },

                  # # Print method #
                  # print = function(){
                  #   clss <- class(self)[1]
                  #   cat(paste0('<',clss,'> object. Package pgr6, version: ',
                  #              private$version,
                  #              '\n'))
                  # },

                  # Basic Subset Methods #
                  # Drop organisms from dataset
                  drop = function(x){
                    orgs <- private$p_organisms
                    if (is.numeric(x)){
                      vec <- c(private$p_dropped, orgs[x])
                      un <- vec[unique(names(vec))]
                      dp <- un[!is.na(un)]
                    }else if (is.character(x)){
                      vec <- c(private$dropped, orgs[orgs%in%x])
                      un <- vec[unique(names(vec))]
                      dp <- un[!is.na(un)]
                    }else{
                      stop('x must be numeric or character')
                    }
                    if (length(dp)) private$p_dropped <- dp else private$p_dropped <- NULL
                    invisible(self)
                  },

                  # Recover from trash previously dropped organisms
                  recover = function(x){
                    orgs <- private$p_organisms
                    if (is.numeric(x)){
                      dp <- private$p_dropped[!names(private$p_dropped)%in%x]
                    } else if (is.character(x)){
                      dp <- private$p_dropped[!private$p_dropped%in%x]
                    }
                    if (length(dp)) private$p_dropped <- dp else private$p_dropped <- NULL
                    invisible(self)
                  }#,

                  # # Get clusters by name
                  # get = function(x){
                  #   aa <- lapply(x, function(y){
                  #     self$clusters[[y]]
                  #   })
                  #   names(aa) <- x
                  #   aa
                  # }

                ),

                # Active binding variables #
                active = list(
                  pan_matrix = function(){
                    pdp <- private$p_dropped
                    if (length(pdp)){
                      pm2 <- private$p_panmatrix[ -as.integer(names(pdp)), , drop = FALSE]
                      cs0 <- colSums(pm2)==0
                      if (any(cs0)){
                        pm2[, -which(cs0), drop = FALSE]
                      }else{
                        pm2
                      }
                    }else{
                      private$p_panmatrix
                    }
                  },

                  organisms = function(){
                    nms <- dimnames(self$pan_matrix)[[1]]
                    pnms <- dimnames(private$p_panmatrix)[[1]]
                    idx <- which(pnms%in%nms)
                    names(nms) <- idx
                    nms
                  },

                  clusters = function(){
                    ogs <- dimnames(self$pan_matrix)[[2]]
                    names(ogs) <- ogs
                    rr <- lapply(ogs, function(x){
                      sset <- private$p_clusters[[x]]#[self$organisms]
                      onlyPublic <- sset[names(sset)%in%self$organisms]
                      onlyPublic[!is.na(onlyPublic)]
                    })
                    class(rr) <- 'ClusterList'
                    attr(rr, 'organisms') <- self$organisms
                    rr
                  },

                  core_level = function(value){
                    if (missing(value)) return(private$p_level)
                    if (value>100) stop("can't set 'core_level' > 100%.")
                    if (value<85) warning("setting 'core_level' under 85%")
                    private$p_level <- value
                    private$p_level
                  },

                  core_clusters = function(){
                    ln <- length(self$organisms)
                    co <- round(private$p_level * ln / 100)
                    pm <- self$pan_matrix
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    wh <- which(colSums(pm) >= co)
                    dimnames(self$pan_matrix)[[2]][wh]
                  },

                  random_seed = function(){
                    if(!exists('.Random.seed')) set.seed(NULL)
                    .Random.seed
                  },

                  dropped = function(){
                    private$p_dropped
                  }

                )
)



## S3 Subset Methods for some active bindings

# Method for subseting the $clusters field as if were $pan_matrix.
`[.ClusterList` <- function(x, i, j){

  Narg <- nargs()
  orgs <- attr(x, 'organisms')

  #simple subset, single arg
  if(Narg<3){
    if (missing(i)) {
      return(x)
    }else{
      rr <- unclass(x)[i]
    }
  }

  #double arg
  #Take all orgs, subset clusters
  if (missing(i)){
    rr <- unclass(x)[j]
  }

  #Take all clusters, subset orgs
  if (missing(j)){
    if (is.character(i)){
      i <- which(orgs%in%i)
    }
    sorgs <- orgs[i]
    rr <- lapply(x, function(y){
      ss <- y[which(names(y)%in%sorgs)]
      ss[!is.na(ss)]
    })

  }


  #Subset orgs and clusters
  if(!missing(i) & !missing(j)){
    rr <- unclass(x)[j]
    if (is.character(i)){
      i <- which(orgs%in%i)
    }
    sorgs <- orgs[i]
    rr <- lapply(rr, function(y){
      ss <- y[which(names(y)%in%sorgs)]
      ss[!is.na(ss)]
    })
  }

  ln <- sapply(rr, length)
  wh <- which(ln==0)
  if (length(wh)) rr[wh] <- NA
  class(rr) <- 'ClusterList'
  rr
}


