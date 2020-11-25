
#' @name PgR6
#' @title PgR6 basic class
#' @description A basic \code{PgR6} class constructor. It contains basic fields
#' and subset functions to handle a pangenome. Final users should use \code{\link{pagoo}}
#' instead of this, since is more easy to understand.
#' @importFrom R6 R6Class
#' @importFrom S4Vectors DataFrame
#' @importFrom reshape2 dcast
# #' @importFrom data.table as.data.table setcolorder dcast
#' @export
PgR6 <- R6Class('PgR6',

                # Private fields #
                private = list(
                  version = NULL,
                  .data = NULL,
                  .panmatrix = NULL,
                  .organisms = NULL,
                  .clusters = NULL,
                  .dropped = NULL,
                  .level = NULL,
                  .sep = NULL
                ),

                # Public functions #
                public = list(
                  #' @description A basic \code{PgR6} class constructor. It contains basic fields
                  #' and subset functions to handle a pangenome.
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
                  #' @param verbose \code{logical}. Whether to display progress messages when loading class.
                  #' @param DF Deprecated. Use \code{data} instead.
                  #' @param group_meta Deprecated. Use \code{cluster_meta} instead.
                  #' @return An R6 object of class PgR6. It contains basic fields and methods for
                  #' analyzing a pangenome.
                  initialize = function(data,
                                        org_meta,
                                        cluster_meta,
                                        core_level = 95,
                                        sep = '__',
                                        verbose = TRUE,
                                        DF,
                                        group_meta){

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

                    # Check data input
                    # 1. Is data.frame or DataFrame?
                    if (verbose) message('Checking class.')
                    cl <- class(data)
                    if (! cl %in% c('data.frame', 'DataFrame', 'DFrame'))
                      stop('"data" is not a data.frame.')

                    # 2. Are colnames correct?
                    if (verbose) message('Checking dimnames.')
                    dn <- dimnames(data)[[2]]
                    if (!all(c('cluster', 'org', 'gene') %in% dn))
                      stop('colnames must have "cluster", "org", and "gene".')

                    # 3. To DataFrame if not yet.
                    if(cl[1]=='data.frame'){
                      data <- DataFrame(data)
                    }

                    # 4. Create gid
                    if (verbose) message('Creating gid (gene ids).')
                    data$gid <- paste(data$org, data$gene, sep = sep)
                    if (any(duplicated(data$gid)))
                      stop('Duplicated gene names in data.')

                    # 5. Column order, just for format
                    dn <- dimnames(data)[[2]]
                    req <- c('cluster', 'org', 'gene', 'gid')
                    extra <- dn[which(!dn%in%req)]
                    data <- data[, c(req, extra)]

                    # # 6. Change class
                    data$cluster <- factor(data$cluster)
                    # data[, group := factor(group)]
                    orgs <- unique(data$org)
                    data$org <- factor(data$org, levels = orgs)
                    data$gene <- factor(data$gene)
                    # data[, org := factor(org, levels = orgs)]

                    # Create organism field. Add metadata (if provided)
                    names(orgs) <- seq_along(orgs)
                    orgs_data <- DataFrame(org=factor(orgs, levels = orgs), row.names = seq_along(orgs))
                    if (!missing(org_meta)){
                      if (verbose) message('Checking provided organism metadata.')
                      if (!class(org_meta)%in%c('data.frame', 'DataFrame', 'DFrame'))
                        stop('"org_meta" should be a data.frame.')
                      if('org'%in%colnames(org_meta)){
                        ma <- match(org_meta$org, orgs_data$org)
                        if (any(is.na(ma))) stop('org_meta$org do not match with data$org')
                        org_meta <- org_meta[ma, ]
                        oc <- which(colnames(org_meta)=='org')
                        orgs_data <- cbind(orgs_data, DataFrame(org_meta[, -oc, drop=F]))
                      }else{
                        warning('"org_meta" should contain an "org" column. Ignoring this parameter.')
                      }
                    }

                    #create group field. Add metadata (if provided)
                    clus_lvls <- levels(data$cluster)
                    cluster_data <- DataFrame(cluster=factor(clus_lvls, levels = clus_lvls))
                    if (!missing(cluster_meta)){
                      if (verbose) message('Checking provided cluster metadata.')
                      if (!class(cluster_meta)%in%c('data.frame', 'DataFrame', 'DFrame'))
                        stop('"cluster_meta" should be a data.frame')
                      if ('cluster'%in%colnames(cluster_meta)){
                        ma <- match(cluster_meta$cluster, cluster_data$cluster)
                        if (any(is.na(ma))) stop('cluster_meta$cluster do not match with data$cluster')
                        cluster_meta <- cluster_meta[ma, ]
                        oc <- which(colnames(cluster_meta)=='cluster')
                        cluster_data <- cbind(cluster_data, DataFrame(cluster_meta[, -oc, drop = F]))
                      }else{
                        warning('"cluster_meta" should contain a "cluster" column. Ignoring this parameter.')
                      }
                    }

                    # Create panmatrix #
                    if (verbose) message('Creating panmatrix.')
                    panmatrix <- dcast(as.data.frame(data),
                                       org~cluster,
                                       value.var = 'gene',
                                       fun.aggregate = length)
                    rn <- panmatrix$org
                    panmatrix$org <- NULL
                    panmatrix <- as.matrix(panmatrix)
                    rownames(panmatrix) <- rn

                    # Populate private$ #
                    if (verbose) message('Populating class.')
                    private$version <- packageVersion('pagoo')
                    private$.data <- data
                    private$.organisms <- orgs_data
                    private$.clusters <- cluster_data
                    private$.panmatrix <- panmatrix
                    private$.level <- core_level
                    private$.sep <- sep
                  },

                  # # Print method #
                  # print = function(){
                  #   clss <- class(self)[1]
                  #   cat(paste0('<',clss,'> object. Package pgr6, version: ',
                  #              private$version,
                  #              '\n'))
                  # },
                  #' @description
                  #' Add metadata to the object. You can add metadata to each organism, to each
                  #' group of orthologous, or to each gene. Elements with missing data should be filled
                  #' by \code{NA} (dimensions of the provided data.frame must be coherent with object
                  #' data).
                  #' @param map \code{character} identifiying the metadata to map. Can
                  #' be one of \code{"org"}, \code{"group"}, or \code{"gid"}.
                  #' @param data \code{data.frame} or \code{DataFrame} with the metadata to
                  #' add. For ethier case, a column named as \code{"map"} must exists, which should
                  #' contain identifiers for each element. In the case of adding gene (\code{gid})
                  #' metadata,each gene should be referenced by the name of the organism and the name
                  #' of the gene as provided in the \code{"data"} data.frame, separated by the
                  #'  \code{"sep"} argument.
                  #' @return \code{self} invisibly, but with additional metadata.
                  add_metadata = function(map = 'org', data){
                    map <- match.arg(map, choices = c('org', 'cluster', 'gid'))
                    cl <- class(data)
                    if (!cl%in%c('data.frame', 'DataFrame')){
                      stop('data should be a data.frame .')
                    }

                    if (map == 'org'){
                      if('org'%in%colnames(data)){
                        ma <- match(private$.organisms$org, data$org)
                        if (any(is.na(ma))) stop('data$org do not match with object organisms.')
                        if (dim(data)[1]!=dim(private$.organisms)[1]){
                          stp <- paste('data has', dim(data)[1], 'rows while object require', dim(private$.organisms)[1],'.')
                          stop(stp)
                        }
                        data <- data[ma, ]
                        oc <- which(colnames(data)=='org')
                        private$.organisms <- cbind(private$.organisms, DataFrame(data[, -oc, drop=F]))
                      }else{
                        stop('"data" should contain an "org" column.')
                      }
                    }else if (map == 'cluster'){
                      if('cluster'%in%colnames(data)){
                        ma <- match(private$.clusters$cluster, data$cluster)
                        if (any(is.na(ma))) stop('data$cluster do not match with object "cluster" column.')
                        if (dim(data)[1]!=dim(private$.clusters)[1]){
                          stp <- paste('data has', dim(data)[1], 'rows while object require', dim(private$.clusters)[1],'.')
                          stop(stp)
                        }
                        data <- data[ma, ]
                        oc <- which(colnames(data)=='cluster')
                        private$.clusters <- cbind(private$.clusters, DataFrame(data[, -oc, drop=F]))
                      }else{
                        stop('"data" should contain an "cluster" column.')
                      }
                    }else{
                      if('gid'%in%colnames(data)){
                        ma <- match(private$.data$gid, data$gid)
                        if (any(is.na(ma))) stop('data$gid do not match with object gid.')
                        if (dim(data)[1]!=dim(private$.data)[1]){
                          stp <- paste('data has', dim(data)[1], 'rows while object require', dim(private$.data)[1],'.')
                          stop(stp)
                        }
                        data <- data[ma, ]
                        oc <- which(colnames(data)=='gid')
                        private$.data <- cbind(private$.data, DataFrame(data[, -oc, drop=F]))
                      }else{
                        stop('"data" should contain a "gid" column.')
                      }
                    }
                    invisible(self)
                  },

                  # Basic Subset Methods #
                  # Drop organisms from dataset
                  #' @description
                  #' Drop an organism from the dataset. This method allows to hide an organism
                  #' from the real dataset, ignoring it in downstream analyses. All the fields and
                  #' methods will behave as it doesn't exist. For instance, if you decide to drop
                  #' organism 1, the \code{$pan_matrix} field (see below) would not show it when
                  #' called.
                  #' @param x  \code{character} or \code{numeric}. The name of the
                  #' organism wanted to be dropped, or its numeric id as returned in
                  #' \code{$organism} field (see below).
                  #' @return \code{self} invisibly, but with \code{x} dropped. It isn't necessary
                  #' to assign the function call to a new object, nor to re-write it as R6 objects
                  #' are mutable.
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
                  #' @description
                  #' Recover a previously \code{$drop()}ped organism (see above). All fields
                  #' and methods will start to behave considering this organism again.
                  #' @param x \code{character} or \code{numeric}. The name of the
                  #' organism wanted to be recover, or its numeric id as returned in
                  #' \code{$dropped} field (see below).
                  #' @return \code{self} invisibly, but with \code{x} recovered. It isn't necessary
                  #' to assign the function call to a new object, nor to re-write it as R6 objects
                  #' are mutable.
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
                  },

                  # Method for saving the pangenome as flat files (text)
                  #' @description
                  #' Write the pangenome data as flat tables (text). Is not the most recomended way
                  #' to save a pangenome, since you can loose information as numeric precision,
                  #' column classes (factor, numeric, integer), and the state of the object itself
                  #' (i.e. dropped organisms, or core_level), loosing reproducibility. Use
                  #' \code{$save_pangenomeRDS} for a more precise way of saveing a pagoo object.
                  #' Still, it is useful if you want to work with the data outside R, just keep
                  #' the above in mind.
                  #' @param dir The unexisting directory name where to put the data files. Default
                  #' is "pangenome".
                  #' @param force \code{logical}. Whether to overwrite the directory if it already
                  #' exists. Default: \code{FALSE}.
                  #' @return A directory with at least 3 files. "data.tsv" contain the basic
                  #' pangenome data as it is provided to the \code{data} argumnet in the
                  #' initialization method (\code{$new(...)}). "clusters.tsv" contain any metadata
                  #' asociated to the clusters. "organisms.tsv" contain any metadata asociated to
                  #' the organisms. The latter 2 files will contain a single column if no metadata
                  #' was provided.
                  write_pangenome = function(dir = "pangenome", force = FALSE){

                    dexist <- dir.exists(dir)
                    if (dexist) {
                      if (force) {
                        unlink(dir, recursive = TRUE)
                        dir.create(dir)
                      } else {
                        stop("dir already exists. Use force = TRUE to overwrite.")
                      }
                    } else {
                      dir.create(dir)
                    }
                    dir <- paste0(normalizePath(dir), "/")
                    # Genes
                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    orgs <- dn[[1]]
                    data <- private$.data
                    act <- which(data$cluster%in%ogs & data$org%in%orgs)
                    x <- data[act, -4L, drop=FALSE] #Don't include "gid" column (4) as it's generated by pagoo
                    write.table(x = x,
                                file = paste0(dir, "data.tsv"),
                                quote = FALSE,
                                row.names = FALSE,
                                col.names = TRUE,
                                sep = "\t")
                    rm(x)
                    # Clusters
                    write.table(self$clusters,
                                file = paste0(dir, "clusters.tsv"),
                                quote = FALSE,
                                row.names = FALSE,
                                col.names = TRUE,
                                sep = "\t")
                    # Organisms
                    write.table(self$organisms,
                              file = paste0(dir, "organisms.tsv"),
                              quote = FALSE,
                              row.names = FALSE,
                              col.names = TRUE,
                              sep = "\t")

                  },

                  #' @description
                  #' Save a pagoo pangenome object. This function provides a method for saving a pagoo
                  #' object and its state into a "RDS" file. To load the pangenome, use the
                  #' \code{load_pangenomeRDS} function in this package. It *should* be compatible between
                  #' pagoo versions, so you could update pagoo and still recover the same pangenome. Even
                  #' \code{sep} and \code{core_level} are restored unless the user provides those
                  #' arguments in \code{load_pangenomeRDS}. \code{dropped} organisms also kept hidden, as
                  #' you where working with the original object.
                  #' @param file The name of the file to save. Default: "pangenome.rds".
                  #' @return Writes a list with all the information needed to restore the object by
                  #' using the load_pangenomeRDS function, into an RDS (binary) file.
                  save_pangenomeRDS = function(file = "pangenome.rds"){
                    clss <- class(self)
                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    orgs <- dn[[1]]
                    data <- private$.data[, -4, drop=FALSE] #Don't include "gid" column (4) as it's generated by pagoo
                    data <- as.data.frame(data)
                    cluster_meta <- as.data.frame(private$.clusters)
                    org_meta <- as.data.frame(private$.organisms)
                    core_level <- self$core_level
                    sep <- private$.sep
                    dropped <- self$dropped
                    version <- private$version

                    pg_data <- list()
                    pg_data$data <- data
                    pg_data$sep <- sep
                    pg_data$dropped <- dropped
                    pg_data$core_level <- core_level
                    if (dim(cluster_meta)[2] > 1){
                      pg_data$cluster_meta <- cluster_meta
                    }
                    if (dim(org_meta)[2] > 1){
                      pg_data$org_meta <- org_meta
                    }
                    if (!is.null(private$.sequences)){
                      sqs <- private$.sequences
                      spl <- split(sqs, mcols(sqs)$org)
                      pg_data$sequences <- lapply(spl, function(x) {
                        patt <- paste0("^.+", sep, collapse = "")
                        names(x) <- sub(patt, "", names(x))
                        as.character(x)
                      })
                      has_seqs <- TRUE
                    } else {
                      has_seqs <- FALSE
                    }
                    # attr(pg_data, "package") <- eval(parse(text=paste0(clss[1], "$parent_env$.packageName")))
                    attr(pg_data, "package") <- environmentName(parent.env(self$.__enclos_env__))
                    attr(pg_data, "parent_package") <- "pagoo"
                    attr(pg_data, "has_seqs") <- has_seqs
                    attr(pg_data, "class") <- clss
                    attr(pg_data, "version") <- version
                    saveRDS(pg_data, file = file)
                  }


                ),


                # Active binding variables #
                active = list(
                  #' @field pan_matrix The panmatrix. Rows are organisms, and
                  #' columns are groups of orthologous. Cells indicates the presence (>=1) or
                  #' absence (0) of a given gene, in a given organism. Cells can have values
                  #' greater than 1 if contain in-paralogs.
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

                  #' @field organisms A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with available
                  #' organism names, and organism number identifier as \code{rownames()}. (Dropped
                  #' organisms will not be displayed in this field, see \code{$dropped} below).
                  #' Additional metadata will be shown if provided, as additional columns.
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

                  #' @field genes A \code{\link[IRanges:DataFrameList-class]{SplitDataFrameList}} object with
                  #' one entry per cluster. Each element contains a \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
                  #' with gene ids (\code{<gid>}) and additional metadata, if provided. \code{gid} are
                  #' created by \code{paste}ing organism and gene names, so duplication in gene names
                  #' are avoided.
                  genes = function(){
                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    orgs <- dn[[1]]
                    data <- private$.data
                    act <- which(data$cluster%in%ogs & data$org%in%orgs)
                    data <- data[act, ]
                    split(data[, , drop = FALSE], f = data$cluster, drop = TRUE)
                  },

                  #' @field clusters A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with the groups
                  #' of orthologous (clusters). Additional metadata will be shown as additional columns,
                  #' if provided before. Each row corresponds to each cluster.
                  clusters = function(){
                    dn <- dimnames(self$pan_matrix)
                    ogs <- dn[[2]]
                    data <- private$.clusters
                    act <- which(data$cluster%in%ogs)
                    data <- data[act, , drop = FALSE]
                    data
                  },

                  #' @field core_level The percentage of organisms a gene must be in
                  #' to be considered as part of the coregenome. \code{core_level = 95} by default.
                  #' Can't be set above 100, and below 85 raises a warning.
                  core_level = function(value){
                    if (missing(value)) return(private$.level)
                    if (value>100) stop("can't set 'core_level' > 100%.")
                    if (value<85) warning("setting 'core_level' under 85%")
                    private$.level <- value
                    private$.level
                  },

                  #' @field core_genes Like \code{genes}, but only showing core genes.
                  core_genes = function(){
                    ln <- length(self$organisms[['org']])
                    co <- floor(private$.level * ln / 100)
                    pm <- self$pan_matrix
                    orgs <- dimnames(pm)[[1]]
                    ogs <- dimnames(pm)[[2]]
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    ogs <- ogs[which(colSums(pm) >= co)]
                    data <- private$.data
                    act <- which(data$cluster%in%ogs & data$org%in%orgs)
                    data <- data[act, ]
                    split(data[, , drop = FALSE], f = data$cluster, drop = TRUE)

                  },

                  #' @field core_clusters Like \code{$clusters}, but only showing core
                  #' clusters.
                  core_clusters = function(){
                    ln <- length(self$organisms[['org']])
                    co <- floor(private$.level * ln / 100)
                    pm <- self$pan_matrix
                    pm[which(pm>1L, arr.ind = TRUE)] <- 1L
                    wh <- which(colSums(pm) >= co)
                    self$clusters[wh, ,drop = FALSE]
                    # dimnames(self$pan_matrix)[[2]][wh]
                  },

                  #' @field cloud_genes Like \code{genes}, but only showing cloud genes.
                  #' These are defined as those clusters which contain a single gene (singletons), plus
                  #' those which have more than one but its organisms are probably clonal due to identical
                  #' general gene content. Colloquially defined as strain-specific genes.
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
                    data <- private$.data
                    act <- which(data$cluster%in%ogs & data$org%in%orgs)
                    data <- data[act, ]
                    split(data[, , drop = FALSE], f = data$cluster, drop = TRUE)
                  },

                  #' @field cloud_clusters Like \code{$clusters}, but only showing cloud
                  #' clusters as defined above.
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

                  #' @field shell_genes Like \code{genes}, but only showing shell genes.
                  #' These are defined as those clusters than don't belong nethier to the core genome,
                  #' nor to cloud genome. Colloquially defined as genes that are present in some but not
                  #' all strains, and that aren't strain-specific.
                  shell_genes = function(){
                    data <- private$.data
                    dn <- dimnames(self$pan_matrix)
                    orgs <- dn[[1]]
                    ogs <- dn[[2]]
                    core <- self$core_clusters$cluster
                    cloud <- self$cloud_clusters$cluster
                    shell <- ogs[which(!ogs %in% unlist(list(core,cloud)))]
                    act <- which(data$cluster%in%shell & data$org%in%orgs)
                    data <- data[act, ]
                    split(data[, , drop = FALSE], f = data$cluster, drop = TRUE)
                  },

                  #' @field shell_clusters Like \code{$clusters}, but only showing shell
                  #' clusters, as defined above.
                  shell_clusters = function(){
                    ogs <- dimnames(self$pan_matrix)[[2]]
                    core <- self$core_clusters$cluster
                    cloud <- self$cloud_clusters$cluster
                    wh <- which(!ogs %in% unlist(list(core,cloud)))
                    self$clusters[wh, , drop = FALSE]
                  },

                  #' @field summary_stats  A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with
                  #' information about the number of core, shell, and cloud clusters, as well as the
                  #' total number of clusters.
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

                  #' @field random_seed The last \code{.Random.seed}. Used for
                  #' reproducibility purposes only.
                  random_seed = function(){
                    if(!exists('.Random.seed')) set.seed(NULL)
                    .Random.seed
                  },

                  #' @field dropped A \code{character} vector with dropped organism
                  #' names, and organism number identifier as \code{names()}
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


