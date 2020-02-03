#' @name PgR6
#' @title PgR6 basic class
#' @description A basic \code{PgR6} class constructor. It contains basic fields
#' and subset functions to handle a pangenome.
#' @section Class Constructor:
#' \describe{
#'     \item{\code{new(DF, org_meta, group_meta, sep = "__")}}{
#'         \itemize{
#'             \item{Create a \code{PgR6} object.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{DF}}: A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing at least the
#'                     following columns: \code{gene} (gene name), \code{org} (organism name to which the gene belongs to),
#'                     and \code{group} (group of orthologous to which the gene belongs to). More columns can be added as metadata
#'                     for each gene.
#'                  }
#'                     \item{\bold{\code{org_meta}}: (optional) A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'                     containging additional metadata for organisms. This \code{data.frame} must have a column named "org" with
#'                     valid organisms names (that is, they should match with those provided in \code{DF}, column \code{org}), and
#'                     additional columns will be used as metadata. Each row should correspond to each organism.
#'
#'                  }
#'                     \item{\bold{\code{group_meta}}: (optional) A \code{data.frame} or \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'                     containging additional metadata for clusters. This \code{data.frame} must have a column named "group" with
#'                     valid organisms names (that is, they should match with those provided in \code{DF}, column \code{group}), and
#'                     additional columns will be used as metadata. Each row should correspond to each cluster.
#'
#'                  }
#'                     \item{\bold{\code{sep}}: A separator. By default is '__'(two underscores). It will be used to
#'                     create a unique \code{gid} (gene identifier) for each gene. \code{gid}s are created by pasting
#'                     \code{org} to \code{gene}, separated by \code{sep}.
#'                  }
#'                     \item{\bold{\code{verbose}}: \code{logical}. Whether to display progress messages when loading class.
#'
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
#'     \item{\code{add_metadata(map = 'org', df)}}{
#'         \itemize{
#'             \item{Add metadata to the object. You can add metadata to each organism, to each
#'             group of orthologous, or to each gene. Elements with missing data should be filled
#'             by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'             data).}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{map}}: \code{character} identifiying the metadata to map. Can
#'                     be one of \code{"org"}, \code{"group"}, or \code{"gid"}.}
#'                     \item{\bold{\code{df}}: \code{data.frame} or \code{DataFrame} with the metadata to
#'                     add. For ethier case, a column named as \code{"map"} must exists, which should
#'                     contain identifiers for each element. In the case of adding gene (\code{gid})
#'                     metadata,each gene should be referenced by the name of the organism and the name
#'                     of the gene as provided in the \code{"DF"} data.frame, separated by the
#'                     \code{"sep"} argument.}
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{\code{self} invisibly, but with additional metadata.}
#'                 }
#'             }
#'         }
#'     }
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
#'     \item{\bold{\code{$pan_matrix}}}{: The panmatrix. Rows are organisms, and
#'     columns are groups of orthologous. Cells indicates the presence (>=1) or
#'     absence (0) of a given gene, in a given organism. Cells can have values
#'     greater than 1 if contain in-paralogs.}
#'     \item{\bold{\code{$organisms}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with available
#'     organism names, and organism number identifier as \code{rownames()}. (Dropped
#'     organisms will not be displayed in this field, see \code{$dropped} below).
#'     Additional metadata will be shown if provided, as additional columns.}
#'     \item{\bold{\code{$clusters}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with the groups
#'     of orthologous (clusters). Additional metadata will be shown as additional columns,
#'     if provided before. Each row corresponds to each cluster.}
#'     \item{\bold{\code{$genes}}}{: A \code{\link[IRanges:DataFrameList-class]{SplitDataFrameList}} object with
#'     one entry per cluster. Each element contains a \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'     with gene ids (\code{<gid>}) and additional metadata, if provided. \code{gid} are
#'     created by \code{paste}ing organism and gene names, so duplication in gene names
#'     are avoided.}
#'     \item{\bold{\code{$core_level}}}{: The percentage of organisms a gene must be in
#'     to be considered as part of the coregenome. \code{core_level = 95} by default.
#'     Can't be set above 100, and below 85 raises a warning.}
#'     \item{\bold{\code{$core_genes}}}{: Like \code{genes}, but only showing core genes.}
#'     \item{\bold{\code{$core_clusters}}}{: Like \code{$clusters}, but only showing core
#'     clusters.}
#'     \item{\bold{\code{$cloud_genes}}}{: Like \code{genes}, but only showing cloud genes.
#'     These are defined as those clusters which contain a single gene (singletons), plus
#'     those which have more than one but its organisms are probably clonal due to identical
#'     general gene content. Colloquially defined as strain-specific genes.}
#'     \item{\bold{\code{$cloud_clusters}}}{: Like \code{$clusters}, but only showing cloud
#'     clusters as defined above.}
#'     \item{\bold{\code{$shell_genes}}}{: Like \code{genes}, but only showing shell genes.
#'     These are defined as those clusters than don't belong nethier to the core genome,
#'     nor to cloud genome. Colloquially defined as genes that are present in some but not
#'     all strains, and that aren't strain-specific.}
#'     \item{\bold{\code{$shell_clusters}}}{: Like \code{$clusters}, but only showing shell
#'     clusters, as defined above.}
#'     \item{\bold{\code{$summary_stats}}}{: A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with
#'     information about the number of core, shell, and cloud clusters, as well as the
#'     total number of clusters.}
#'     \item{\bold{\code{$random_seed}}}{: The last \code{.Random.seed}. Used for
#'     reproducibility purposes only.}
#'     \item{\bold{\code{$dropped}}}{: A \code{character} vector with dropped organism
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
#' @importFrom S4Vectors DataFrame
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
                  .groups = NULL,
                  .dropped = NULL,
                  .level = NULL,
                  .sep = NULL
                ),

                # Public functions #
                public = list(

                  initialize = function(DF,
                                        org_meta,
                                        group_meta,
                                        sep = '__',
                                        verbose = TRUE){

                    # Check DF input
                    # 1. Is data.frame or DataFrame?
                    if (verbose) message('Checking class.')
                    cl <- class(DF)
                    if (! cl %in% c('data.frame', 'DataFrame'))
                      stop('"DF" is not a data.frame.')

                    # 2. Are colnames correct?
                    if (verbose) message('Checking dimnames.')
                    dn <- dimnames(DF)[[2]]
                    if (!all(c('group', 'org', 'gene') %in% dn))
                      stop('colnames must have "group", "org", and "gene".')

                    # 3. To DataFrame if not yet.
                    if(cl[1]=='data.frame'){
                      DF <- DataFrame(DF)
                    }

                    # 4. Create gid
                    if (verbose) message('Creating gid (gene ids).')
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
                    orgs <- unique(DF$org)
                    DF$org <- factor(DF$org, levels = orgs)
                    # DF[, org := factor(org, levels = orgs)]

                    # Create organism field. Add metadata (if provided)
                    names(orgs) <- seq_along(orgs)
                    orgs_DF <- DataFrame(org=orgs, row.names = seq_along(orgs))
                    if (!missing(org_meta)){
                      if (verbose) message('Checking provided organism metadata.')
                      if (!class(org_meta)%in%c('data.frame', 'DataFrame'))
                        stop('"org_meta" should be a data.frame.')
                      if('org'%in%colnames(org_meta)){
                        ma <- match(org_meta$org, orgs_DF$org)
                        if (any(is.na(ma))) stop('org_meta$org do not match with DF$org')
                        org_meta <- org_meta[ma, ]
                        oc <- which(colnames(org_meta)=='org')
                        orgs_DF <- cbind(orgs_DF, DataFrame(org_meta[, -oc, drop=F]))
                      }else{
                        warning('"org_meta" should contain an "org" column. Ignoring this parameter.')
                      }
                    }

                    #create group field. Add metadata (if provided)
                    group_DF <- DataFrame(group=levels(DF$group))
                    if (!missing(group_meta)){
                      if (verbose) message('Checking provided group metadata.')
                      if (!class(group_meta)%in%c('data.frame', 'DataFrame'))
                        stop('"group_meta" should be a data.frame')
                      if ('group'%in%colnames(group_meta)){
                        ma <- match(group_meta$group, group_DF$group)
                        if (any(is.na(ma))) stop('grop_meta$group do not match with DF$group')
                        group_meta <- group_meta[ma, ]
                        oc <- which(colnames(group_meta)=='group')
                        group_DF <- cbind(group_DF, DataFrame(group_meta[, -oc, drop = F]))
                      }else{
                        warning('"group_meta" should contain a "group" column. Ignoring this parameter.')
                      }
                    }

                    # Create panmatrix #
                    if (verbose) message('Creating panmatrix.')
                    panmatrix <- dcast(as.data.frame(DF),
                                       org~group,
                                       value.var = 'gene',
                                       fun.aggregate = length)
                    rn <- panmatrix$org
                    panmatrix$org <- NULL
                    panmatrix <- as.matrix(panmatrix)
                    rownames(panmatrix) <- rn

                    # Populate private$ #
                    if (verbose) message('Populating class.')
                    private$version <- packageVersion('pagoo')
                    private$.DF <- DF
                    private$.organisms <- orgs_DF
                    private$.groups <- group_DF
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

                  add_metadata = function(map = 'org', df){
                    map <- match.arg(map, choices = c('org', 'group', 'gid'))
                    cl <- class(df)
                    if (!cl%in%c('data.frame', 'DataFrame')){
                      stop('df should be a data.frame .')
                    }

                    if (map == 'org'){
                      if('org'%in%colnames(df)){
                        ma <- match(df$org, private$.organisms$org)
                        if (any(is.na(ma))) stop('df$org do not match with object organisms.')
                        if (dim(df)[1]!=dim(private$.organisms)[1]){
                          stp <- paste('df has', dim(df)[1], 'rows while object require', dim(private$.organisms)[1],'.')
                          stop(stp)
                        }
                        df <- df[ma, ]
                        oc <- which(colnames(df)=='org')
                        private$.organisms <- cbind(private$.organisms, DataFrame(df[, -oc, drop=F]))
                      }else{
                        stop('"df" should contain an "org" column.')
                      }
                    }else if (map == 'group'){
                      if('group'%in%colnames(df)){
                        ma <- match(df$group, private$.groups$group)
                        if (any(is.na(ma))) stop('df$group do not match with object groups.')
                        if (dim(df)[1]!=dim(private$.groups)[1]){
                          stp <- paste('df has', dim(df)[1], 'rows while object require', dim(private$.groups)[1],'.')
                          stop(stp)
                        }
                        df <- df[ma, ]
                        oc <- which(colnames(df)=='group')
                        private$.groups <- cbind(private$.groups, DataFrame(df[, -oc, drop=F]))
                      }else{
                        stop('"df" should contain an "group" column.')
                      }
                    }else{
                      if('gid'%in%colnames(df)){
                        ma <- match(df$gid, private$.DF$gid)
                        if (any(is.na(ma))) stop('df$gid do not match with object gid.')
                        if (dim(df)[1]!=dim(private$.DF)[1]){
                          stp <- paste('df has', dim(df)[1], 'rows while object require', dim(private$.DF)[1],'.')
                          stop(stp)
                        }
                        df <- df[ma, ]
                        oc <- which(colnames(df)=='gid')
                        private$.DF <- cbind(private$.DF, DataFrame(df[, -oc, drop=F]))
                      }else{
                        stop('"df" should contain a "gid" column.')
                      }
                    }
                    invisible(self)
                  },

                  # Basic Subset Methods #
                  # Drop organisms from dataset
                  drop = function(x){
                    orgs <- as.character(private$.organisms[['org']])
                    orgs <- setNames(orgs, seq_along(orgs))
                    if (is.numeric(x)){
                      vec <- c(private$.dropped, orgs[x])
                      un <- vec[unique(names(vec))]
                      dp <- un[!is.na(un)]
                    }else if (is.character(x)){
                      vec <- c(private$.dropped, orgs[orgs%in%x])
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
                    orgs <- as.character(private$.organisms[['org']])
                    orgs <- setNames(orgs, seq_along(orgs))
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

                  organisms = function(){
                    orgs <- private$.organisms
                    drp <- private$.dropped
                    if (length(drp)>0){
                      idx <- as.integer(names(drp))
                      orgs[-idx, ,drop = FALSE]
                    }else{
                      orgs
                    }
                  },

                  genes = function(){
                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    orgs <- dn[[1]]
                    df <- private$.DF
                    act <- which(df$group%in%ogs & df$org%in%orgs)
                    df <- df[act, ]
                    split(df[, -c(1,2,3), drop = FALSE], f = df$group, drop = TRUE)
                  },

                  clusters = function(){
                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    df <- private$.groups
                    act <- which(df$group%in%ogs)
                    df <- df[act, , drop = FALSE]
                    df
                  },

                  core_level = function(value){
                    if (missing(value)) return(private$.level)
                    if (value>100) stop("can't set 'core_level' > 100%.")
                    if (value<85) warning("setting 'core_level' under 85%")
                    private$.level <- value
                    private$.level
                  },

                  core_genes = function(){
                    ln <- length(self$organisms[['org']])
                    co <- round(private$.level * ln / 100)
                    pm <- self$pan_matrix
                    orgs <- dimnames(pm)[[1]]
                    ogs <- dimnames(pm)[[2]]
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    ogs <- ogs[which(colSums(pm) >= co)]
                    df <- private$.DF
                    act <- which(df$group%in%ogs & df$org%in%orgs)
                    df <- df[act, ]
                    split(df[, -c(1,2,3), drop = FALSE], f = df$group, drop = TRUE)

                  },

                  core_clusters = function(){
                    ln <- length(self$organisms[['org']])
                    co <- round(private$.level * ln / 100)
                    pm <- self$pan_matrix
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    wh <- which(colSums(pm) >= co)
                    self$clusters[wh, ,drop = FALSE]
                    # dimnames(self$pan_matrix)[[2]][wh]
                  },

                  cloud_genes = function(){
                    pm <- self$pan_matrix
                    orgs <- dimnames(pm)[[1]]
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    dups <- duplicated(pm, MARGIN = 1)
                    wdups <- which(dups)
                    if (length(wdups)){
                      pm <- pm[-wdups, ]
                    }
                    ogs <- dimnames(pm)[[2]]
                    ogs <- ogs[which(colSums(pm) == 1L)]
                    df <- private$.DF
                    act <- which(df$group%in%ogs & df$org%in%orgs)
                    df <- df[act, ]
                    split(df[, -c(1,2,3), drop = FALSE], f = df$group, drop = TRUE)
                  },

                  cloud_clusters = function(){
                    pm <- self$pan_matrix
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    dups <- duplicated(pm, MARGIN = 1)
                    wdups <- which(dups)
                    if (length(wdups)){
                      pm <- pm[-wdups, ]
                    }
                    wh <- which(colSums(pm) == 1L)
                    self$clusters[wh, , drop = FALSE]
                  },

                  shell_genes = function(){
                    df <- private$.DF
                    dn <- dimnames(self$pan_matrix)
                    orgs <- dn[[1]]
                    ogs <- dn[[2]]
                    core <- self$core_clusters$group
                    cloud <- self$cloud_clusters$group
                    shell <- ogs[which(!ogs %in% c(core,cloud))]
                    act <- which(df$group%in%shell & df$org%in%orgs)
                    df <- df[act, ]
                    split(df[, -c(1,2,3), drop = FALSE], f = df$group, drop = TRUE)
                  },

                  shell_clusters = function(){
                    ogs <- dimnames(self$pan_matrix)[[2]]
                    core <- self$core_clusters$group
                    cloud <- self$cloud_clusters$group
                    wh <- which(!ogs %in% c(core,cloud))
                    self$clusters[wh, , drop = FALSE]
                  },

                  summary_stats = function(){
                    total <- dim(self$pan_matrix)[2]
                    core <- dim(self$core_clusters)[1]
                    cloud <- dim(self$cloud_clusters)[1]
                    shell <- total - core - cloud
                    DataFrame(Category = c('Total',
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
                               nnorg <- x$organisms[['org']]

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
