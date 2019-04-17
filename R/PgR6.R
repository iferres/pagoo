# Basic pgr6 class
# This is a basic pgr6 class constructor from which inherit more complex ones.

# Test: bare input
# cluster_list <- list(structure(c("gene1", "gene1", "gene1"),
#                                .Names = c("org1","org2", "org3")),
#                      structure(c("gene2", "gene3", "gene2"),
#                                .Names = c("org1","org1", "org3")),
#                      structure(c("gene4", "gene2"),
#                                .Names = c("org1","org2")))

# DF <- structure(list(group = c("group1", "group1", "group1", "group2",
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
#' @section Class Constructor:
#' \describe{
#'     \item{\code{new(DF, sep = "__")}}{
#'         \itemize{
#'             \item{Create a \code{PgR6} object.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{DF}}: A \code{data.frame} or \code{data.table} containing at least the
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
#' @importFrom R6 R6Class
#' @import S4Vectors
#' @importFrom reshape2 dcast
# #' @importFrom data.table as.data.table setcolorder dcast
#' @export
PgR6 <- R6Class('PgR6',

                # Private fields #
                private = list(
                  version = NULL,
                  .DF = NULL,
                  .panmatrix = NULL,
                  .organisms = NULL,
                  .dropped = NULL,
                  .level = NULL,
                  .sep = NULL
                ),

                # Public functions #
                public = list(

                  initialize = function(DF,
                                        org_meta,
                                        group_meta,
                                        sep = '__'){

                    # Check DF input
                    # 1. Is data.frame or DataFrame?
                    cl <- class(DF)
                    if (! cl %in% c('data.frame', 'DataFrame'))
                      stop('"DF" is not a data.frame.')

                    # 2. Are colnames correct?
                    dn <- dimnames(DF)[[2]]
                    if (!all(c('group', 'org', 'gene') %in% dn))
                      stop('colnames must have "group", "org", and "gene".')

                    # 3. To DataFrame if not yet.
                    if(cl[1]=='data.frame'){
                      DF <- DataFrame(DF)
                    }

                    # 4. Create gid
                    DF$gid <- paste(DF$org, DF$gene, sep = sep)
                    if (any(duplicated(DF$gid)))
                      stop('Duplicated gene names in DF.')

                    # 5. Column order, just for format
                    dn <- dimnames(DF)[[2]]
                    req <- c('group', 'org', 'gene', 'gid')
                    extra <- dn[which(!dn%in%req)]
                    DF <- DF[, c(req, extra)]

                    # # 6. Change class
                    DF$group <- factor(DF$group)
                    # DF[, group := factor(group)]
                    orgs <- levels(DF$org)
                    DF$org <- factor(DF$org, levels = orgs)
                    # DF[, org := factor(org, levels = orgs)]

                    # Create organism field. Add metadata (if provided)
                    names(orgs) <- seq_along(orgs)
                    orgs_DF <- DataFrame(org=orgs, row.names = seq_along(orgs))
                    if (!missing(org_meta)){
                      if (!class(org_meta)%in%c('data.frame', 'DataFrame'))
                        stop('"org_meta" should be a data.frame.')
                      if('org'%in%colnames(org_meta)){
                        ma <- match(org_meta$org, orgs_DF$org)
                        if (any(is.na(ma))) stop('org_meta$org do not match with DF$org')
                        org_meta <- org_meta[ma, ]
                        oc <- which(colnames(org_meta)=='org')
                        orgs <- cbind(orgs_DF, DataFrame(org_meta[, -oc, drop=F]))
                      }else{
                        warning('"org_meta" should contain an "org" column. Ignoring this parameter.')
                      }
                    }

                    # Create panmatrix #
                    panmatrix <- dcast(as.data.frame(DF),
                                       org~group,
                                       value.var = 'gene',
                                       fun.aggregate = length)
                    rn <- panmatrix$org
                    panmatrix$org <- NULL
                    panmatrix <- as.matrix(panmatrix)
                    rownames(panmatrix) <- rn

                    # Populate private$ #
                    private$version <- packageVersion('pgr6')
                    private$.DF <- DF
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
                    orgs <- private$.organisms[, 1, drop=TRUE]
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
                    orgs <- private$.organisms[, 1, drop=TRUE]
                    if (is.numeric(x)){
                      dp <- private$.dropped[!names(private$.dropped)%in%x]
                    } else if (is.character(x)){
                      dp <- private$.dropped[!private$.dropped%in%x]
                    }
                    if (length(dp)) private$.dropped <- dp else private$.dropped <- NULL
                    invisible(self)
                  }#,


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

                  organisms = function(value){
                    drp <- private$.dropped
                    if (length(drp)>0){
                      idx <- as.integer(names(drp))
                      orgs[-idx, , drop=FALSE]
                    }else{
                      orgs
                    }
                  },

                  clusters = function(){

                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    orgs <- dn[[1]]
                    df <- private$.DF
                    act <- which(df$group%in%ogs & df$org%in%orgs)
                    df <- df[act, ]
                    rr <- split(df[, -c(1,2,3)], f = df$group , drop = TRUE)
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
                    ln <- length(self$organisms[, 1, drop=TRUE])
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

                           cloneable = FALSE, class = FALSE, portable = TRUE

)

## S3 Subset Methods for some active bindings
#' @export
`[.PgR6` <- function(self, i, j){
  Nargs <- nargs()
  mj <- missing(j)
  isSimple <- (Nargs + !(mj))<3
  asst <- names(.active_subset)
  cl_abind <- which(sapply(names(self), bindingIsActive, env = self))
  avail <- which(names(cl_abind) %in% asst)
  addm <- names(cl_abind)[avail]
  lapply(addm, function(x) PgR6_subset$set('active', x, .active_subset[[x]], overwrite = TRUE))
  on.exit(lapply(addm, function(x) PgR6_subset$active[[x]] <- NULL ))
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
