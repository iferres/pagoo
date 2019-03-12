# Basic pgr6 class
# This is a basic pgr6 class constructor from which inherit more complex ones.

# Test: bare input
# cluster_list <- list(structure(c("gene1", "gene1", "gene1"),
#                                .Names = c("org1","org2", "org3")),
#                      structure(c("gene2", "gene3", "gene2"),
#                                .Names = c("org1","org1", "org3")),
#                      structure(c("gene4", "gene2"),
#                                .Names = c("org1","org2")))

# cluster_df <- structure(list(group = c("group1", "group1", "group1", "group2",
#                                        "group2", "group2", "group3", "group3"),
#                              org = c("org1", "org2", "org3", "org1",
#                                      "org1", "org3", "org1", "org2"),
#                              gene = c("gene1","gene1", "gene1", "gene2",
#                                       "gene3", "gene2", "gene4", "gene2")),
#                         .Names = c("group", "org", "gene"),
#                         row.names = c(NA, -8L),class = "data.frame")

#' @name PgR6
#' @title PgR6 basic class
#' @description A basic \code{PgR6} class constructor. It contains basic fields
#' and subset functions to handle a pangenome.
#' @param cluster_df A \code{data.frame} or \code{data.table} cointaining at
#' least the following 3 columns: \code{"gene"}, with gene name; \code{"org"},
#' organism name; and \code{"group"} with the name of the cluster it belongs
#' to.
#' @param sep A separator for creating a \code{"gid"} (gene id) column. By
#' default is "__" (two underscores). \code{gid} are created by pasting the
#' organism name to the gene name, in order to have unique gene identifiers.
#'
#' @importFrom R6 R6Class
#' @importFrom data.table as.data.table setcolorder dcast
#' @export
PgR6 <- R6Class('PgR6',

                # Private fields #
                private = list(
                  version = NULL,
                  .dt = NULL,
                  .panmatrix = NULL,
                  .organisms = NULL,
                  .dropped = NULL,
                  .level = NULL,
                  .sep = NULL
                ),

                # Public functions #
                public = list(

                  initialize = function(cluster_df,
                                        sep = '__'){

                    # Check cluster_df input
                    # 1. Is data.frame or data.table?
                    cl <- class(cluster_df)
                    if (! cl %in% c('data.frame', 'data.table'))
                      stop('"cluster_df" is not a data.frame nor a data.table.')

                    # 2. Are colnames correct?
                    dn <- dimnames(cluster_df)[[2]]
                    if (!all(c('group', 'org', 'gene') %in% dn))
                      stop('colnames must have "group", "org", and "gene".')

                    # 3. To data.table if not yet.
                    if(cl[1]=='data.frame'){
                      cluster_df <- as.data.table(cluster_df)
                    }

                    # 4. Create gid
                    cluster_df[, gid := paste(org, gene, sep = sep)]
                    if (any(duplicated(cluster_df[, gid])))
                      stop('Duplicated gene names in cluster_df.')

                    # 5. Column order, just for format
                    setcolorder(cluster_df, c('group','gid', 'org', 'gene'))

                    # # 6. Change class
                    cluster_df[, group := factor(group)]
                    orgs <- unique(cluster_df$org)
                    cluster_df[, org := factor(org, levels = orgs)]

                    ############# OLD ###############

                    # # Check cluster_list input #
                    # # 1. Is list. length > 0?
                    # ln_cluster_list <- length(cluster_list)
                    # if (any(c(!is.list(cluster_list), ln_cluster_list<=1)))
                    #   stop('please, check input. Must provide a list of length > 1')
                    # # 2. Check cluster_list names. Create if !exist.
                    # if (!length(names(cluster_list))){
                    #   nch <- nchar(ln_cluster_list)
                    #   np <- paste0('%0', nch, 'd')
                    #   names(cluster_list) <- paste0(prefix, sprintf(np, 1:ln_cluster_list))
                    # }
                    # # 3. Check if organism names are given.
                    # hasOrgNames <- sapply(cluster_list, function(x) all(length(names(x))) )
                    # if (!all(hasOrgNames))
                    #   stop('you must provide organism names for all genes')
                    # # 4. Modify gene names
                    # cluster_list <- lapply(cluster_list, function(x){
                    #   nn <- names(x)
                    #   ng <- paste(nn, sep, x, sep = '')
                    #   names(ng) <- nn
                    #   ng
                    # })
                    # private$.sep <- sep
                    #################################

                    # Create panmatrix #
                    panmatrix <- dcast(cluster_df,
                                       org~group,
                                       value.var = 'gene',
                                       fun.aggregate = length)
                    rn <- panmatrix$org
                    panmatrix$org <- NULL
                    panmatrix <- as.matrix(panmatrix)
                    rownames(panmatrix) <- rn


                    ############# OLD ###############
                    # # 1. Which organisms are?
                    # organisms <- unique(unlist(lapply(cluster_list, names)))
                    # names(organisms) <- 1:length(organisms)
                    # # 2. Compute panmatrix
                    # panmatrix <- sapply(cluster_list, function(x){
                    #   table(factor(names(x), levels = organisms))
                    # }, simplify = TRUE)
                    ################################

                    # Populate private$ #
                    private$version <- packageVersion('pgr6')
                    # private$.clusters <- list2env(cluster_list, hash = TRUE)
                    private$.dt <- cluster_df
                    # private$.organisms <- organisms
                    names(orgs) <- seq_along(orgs)
                    private$.organisms <- orgs
                    private$.panmatrix <- panmatrix
                    private$.level <- 95 #default
                    private$.sep <- sep
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
                    orgs <- private$.organisms
                    if (is.numeric(x)){
                      vec <- c(private$.dropped, orgs[x])
                      un <- vec[unique(names(vec))]
                      dp <- un[!is.na(un)]
                    }else if (is.character(x)){
                      vec <- c(private$dropped, orgs[orgs%in%x])
                      un <- vec[unique(names(vec))]
                      dp <- un[!is.na(un)]
                    }else{
                      stop('x must be numeric or character')
                    }
                    if (length(dp)) private$.dropped <- dp else private$.dropped <- NULL
                    invisible(self)
                  },

                  # Recover from trash previously dropped organisms
                  recover = function(x){
                    orgs <- private$.organisms
                    if (is.numeric(x)){
                      dp <- private$.dropped[!names(private$.dropped)%in%x]
                    } else if (is.character(x)){
                      dp <- private$.dropped[!private$.dropped%in%x]
                    }
                    if (length(dp)) private$.dropped <- dp else private$.dropped <- NULL
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
                    pdp <- private$.dropped
                    if (length(pdp)){
                      pm2 <- private$.panmatrix[ -as.integer(names(pdp)), , drop = FALSE]
                      cs0 <- colSums(pm2)==0
                      if (any(cs0)){
                        pm2[, -which(cs0), drop = FALSE]
                      }else{
                        pm2
                      }
                    }else{
                      private$.panmatrix
                    }
                  },

                  organisms = function(){
                    nms <- dimnames(self$pan_matrix)[[1]]
                    pnms <- dimnames(private$.panmatrix)[[1]]
                    idx <- which(pnms%in%nms)
                    names(nms) <- idx
                    nms
                  },

                  # clusters = function(){
                  #   ogs <- dimnames(self$pan_matrix)[[2]]
                  #   names(ogs) <- ogs
                  #   rr <- lapply(ogs, function(x){
                  #     sset <- private$.clusters[[x]]#[self$organisms]
                  #     onlyPublic <- sset[names(sset)%in%self$organisms]
                  #     onlyPublic[!is.na(onlyPublic)]
                  #   })
                  #   class(rr) <- 'ClusterList'
                  #   attr(rr, 'organisms') <- self$organisms
                  #   attr(rr, 'separator') <- private$.sep
                  #   rr
                  # },

                  clusters = function(){

                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    orgs <- dn[[1]]
                    rr <- split(private$.dt[group%in%ogs & org%in%orgs, ],
                                by = 'group', keep.by = F)
                    # class(rr) <- 'ClusterList'
                    # attr(rr, 'organisms') <- self$organisms
                    # attr(rr, 'separator') <- private$.sep
                    # attr(rr, 'pgr6env') <- parent.frame()
                    rr
                  },

                  core_level = function(value){
                    if (missing(value)) return(private$.level)
                    if (value>100) stop("can't set 'core_level' > 100%.")
                    if (value<85) warning("setting 'core_level' under 85%")
                    private$.level <- value
                    private$.level
                  },

                  core_clusters = function(){
                    ln <- length(self$organisms)
                    co <- round(private$.level * ln / 100)
                    pm <- self$pan_matrix
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    wh <- which(colSums(pm) >= co)
                    dimnames(self$pan_matrix)[[2]][wh]
                  },

                  cloud_clusters = function(){
                    pm <- self$pan_matrix
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    dups <- duplicated(pm, MARGIN = 1)
                    wdups <- which(dups)
                    if (length(wdups)){
                      pm <- pm[-wdups, ]
                    }
                    dimnames(pm)[[2]][which(colSums(pm) == 1L)]
                  },

                  shell_clusters = function(){
                    dn <- dimnames(self$pan_matrix)[[2]]
                    dn[which(!dn %in% c(self$core_clusters, self$cloud_clusters))]
                  },

                  summary_stats = function(){
                    total <- dim(self$pan_matrix)[2]
                    core <- length(self$core_clusters)
                    cloud <- length(self$cloud_clusters)
                    shell <- total - core - cloud
                    data.frame(Category = c('Total',
                                            'Core',
                                            'Shell',
                                            'Cloud'),
                               Number = c(total,
                                          core,
                                          shell,
                                          cloud))
                  },

                  random_seed = function(){
                    if(!exists('.Random.seed')) set.seed(NULL)
                    .Random.seed
                  },

                  dropped = function(){
                    private$.dropped
                  }

                )
)


# #' @importFrom data.table `[.data.table`
PgR6_subset <- R6::R6Class('PgR6_subset',

                           private = list(.x = NULL,
                                          .i = NULL,
                                          .j = NULL,
                                          .simple = NULL
                                          ),


                           public = list(
                             initialize = function(x, i, j, .simple){

                               private$.simple <- .simple

                               nnclu <- dimnames(x$pan_matrix)[[2]]
                               names(nnclu) <- seq_along(nnclu)
                               nnorg <- x$organisms

                               # Simple subset
                               if(.simple){
                                 j <- i
                                 private$.i <- nnorg
                                 if (is.numeric(j)){
                                   private$.j <- nnclu[j]
                                 }else{
                                   wh <- which(nnclu%in%j)
                                   private$.j <- nnclu[wh]
                                 }


                               # Double arg
                               }else{

                                 ## Take all orgs, subset clusters
                                 if (missing(i)){
                                   private$.i <- nnorg
                                   if (is.numeric(j)){
                                     private$.j <- nnclu[j]
                                   }else{
                                     wh <- which(nnclu%in%j)
                                     private$.j <- nnclu[wh]
                                   }

                                 }

                                 ## Take all cluster, subset orgs
                                 if (missing(j)){
                                   private$.j <- nnclu
                                   if (is.numeric(i)){
                                     private$.i <- nnorg[i]
                                   }else{
                                     wh <- which(nnorg%in%i)
                                     private$.i <- nnorg[wh]
                                   }
                                 }

                                 ## Subset cluster and orgs
                                 if (!missing(i) & !missing(j)){
                                   if (is.numeric(i)){
                                     private$.i <- nnorg[i]
                                   }else{
                                     wh <- which(nnorg%in%i)
                                     private$.i <- nnorg[wh]
                                   }
                                   if (is.numeric(j)){
                                     private$.j <- nnclu[j]
                                   }else{
                                     wh <- which(nnclu%in%j)
                                     private$.j <- nnclu[wh]
                                   }
                                 }
                               }

                               private$.x <- x

                             }
                           ),

                           cloneable = FALSE, class = FALSE, portable = FALSE

)

## S3 Subset Methods for some active bindings
#' @export
`[.PgR6` <- function(self, i, j){
  Nargs <- nargs()
  mj <- missing(j)
  isSimple <- (Nargs + !(mj))<3
  asst <- names(active_subset)
  cl_abind <- which(sapply(names(self), bindingIsActive, env = self))
  avail <- which(names(cl_abind) %in% asst)
  addm <- names(cl_abind)[avail]
  lapply(addm, function(x) PgR6_subset$set('active', x, active_subset[[x]], overwrite = TRUE))
  # on.exit(lapply(addm, function(x) PgR6_subset$set('active', x, NULL, overwrite = TRUE)))
  rr <- PgR6_subset$new(x = self, i, j, .simple = isSimple)
  invisible(rr)
}

# # Method for subseting the $clusters field as if were $pan_matrix.
# `[.ClusterList` <- function(x, i, j){
#
#   Narg <- nargs()
#   # orgs <- attr(x, 'organisms')
#   orgs <- self$organisms
#
#   #simple subset, single arg
#   if(Narg<3){
#     if (missing(i)) {
#       return(x)
#     }else{
#       rr <- unclass(x)[i]
#     }
#   }
#
#   #double arg
#   #Take all orgs, subset clusters
#   if (missing(i)){
#     rr <- unclass(x)[j]
#   }
#
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
#
#   }
#
#
#   #Subset orgs and clusters
#   if(!missing(i) & !missing(j)){
#     rr <- unclass(x)[j]
#     if (is.character(i)){
#       i <- which(orgs%in%i)
#     }
#     sorgs <- orgs[i]
#     rr <- lapply(rr, function(y){
#       ss <- y[which(names(y)%in%sorgs)]
#       ss[!is.na(ss)]
#     })
#   }
#
#   ln <- sapply(rr, length)
#   wh <- which(ln==0)
#   if (length(wh)) rr[wh] <- NA
#   class(rr) <- 'ClusterList'
#   rr
# }
#
#
