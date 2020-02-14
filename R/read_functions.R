#' @name pagoo
#' @title Create a Pagoo Object
#' @description This is the main function to load a pagoo object. It's safer and
#' more friendly than using pagoo's class constructors (\link{PgR6},
#' \link{PgR6M}, and \link{PgR6MS}). This function returns either a
#' \code{\link{PgR6M}} class object, or a \code{\link{PgR6MS}} class object,
#' depending on the parameters provided. If sequences are provided, it returns
#' the latter. See below for more details.
#' @param data A \code{data.frame} or
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containing at least the
#'   following columns: \code{gene} (gene name), \code{org} (organism name to
#'   which the gene belongs to), and \code{cluster} (group of orthologous to
#'   which the gene belongs to). More columns can be added as metadata for each
#'   gene.
#' @param org_meta (optional) A \code{data.frame} or
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containging additional
#'   metadata for organisms. This \code{data.frame} must have a column named
#'   "org" with valid organisms names (that is, they should match with those
#'   provided in \code{data}, column \code{org}), and additional columns will be
#'   used as metadata. Each row should correspond to each organism.
#' @param cluster_meta (optional) A \code{data.frame} or
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} containging additional
#'   metadata for clusters. This \code{data.frame} must have a column named
#'   "cluster" with valid organisms names (that is, they should match with those
#'   provided in \code{data}, column \code{cluster}), and additional columns
#'   will be used as metadata. Each row should correspond to each cluster.
#' @param sequences (optional) Can accept: 1) a named \code{list} of named
#'   \code{character} vector. Name of list are names of organisms, names of
#'   character vector are gene names; or 2) a named \code{list} of
#'   \code{\link[Biostrings:XStringSet-class]{DNAStringSetList}} objects (same
#'   requirements as (1), but with BStringSet names as gene names); or 3) a
#'   \code{\link[Biostrings:XStringSetList-class]{DNAStringSetList}} (same
#'   requirements as (2) but \code{DNAStringSetList} names are organisms names).
#'   If this parameter is used, then a \code{\link{PgR6MS}} class object is
#'   returned.
#' @param core_level The initial core_level (thats the percentage of organisms a
#'   core cluster must be in to be considered as part of the core genome). Must
#'   be a number between 100 and 85, (default: 95). You can change it later by
#'   using the \code{$core_level} field once the object was created.
#' @param sep A separator. By default is '__'(two underscores). It will be used
#'   to create a unique \code{gid} (gene identifier) for each gene. \code{gid}s
#'   are created by pasting \code{org} to \code{gene}, separated by \code{sep}.
#' @param verbose \code{logical}. Whether to display progress messages when
#'   loading class.
#' @details This package uses
#' [R6](https://r6.r-lib.org/articles/Introduction.html) classes to provide a
#' unified, comprehensive, standarized, but at the same time flexible, way to
#' analyze a pangenome. The idea is to have a single object which contains both
#' the data and the basic methods to analyze them, as well as manipulate fields,
#' explore, and to use in harmony with the already existing extensive list of R
#' packages available created for comparative genomics and genetics.
#'
#' For more information, tutorials, and resources, please visit XXXXXXXXXXXXXXXX.
#'
#' Above you will find a detailed description of each method/field provided
#' once you load a pagoo object.
#'
#' @section Active bindings:
#' \describe{\itemize{
#'     \item{\bold{\code{$pan_matrix}}}{ The panmatrix. Rows are organisms, and
#'     columns are groups of orthologous. Cells indicates the presence (>=1) or
#'     absence (0) of a given gene, in a given organism. Cells can have values
#'     greater than 1 if contain in-paralogs.}
#'     \item{\bold{\code{$organisms}}}{ A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with available
#'     organism names, and organism number identifier as \code{rownames()}. (Dropped
#'     organisms will not be displayed in this field, see \code{$dropped} below).
#'     Additional metadata will be shown if provided, as additional columns.}
#'     \item{\bold{\code{$clusters}}}{ A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with the groups
#'     of orthologous (clusters). Additional metadata will be shown as additional columns,
#'     if provided before. Each row corresponds to each cluster.}
#'     \item{\bold{\code{$genes}}}{ A \code{\link[IRanges:DataFrameList-class]{SplitDataFrameList}} object with
#'     one entry per cluster. Each element contains a \code{\link[S4Vectors:DataFrame-class]{DataFrame}}
#'     with gene ids (\code{<gid>}) and additional metadata, if provided. \code{gid} are
#'     created by \code{paste}ing organism and gene names, so duplication in gene names
#'     are avoided.}
#'     \item{\bold{\code{$sequences}}}{ A \code{\link[Biostrings:XStringSetList-class]{DNAStringSetList}} with the
#'     set of sequences grouped by cluster. Each group is accessible as were a list. All
#'     \code{Biostrings} methods are available.}
#'     \item{\bold{\code{$core_level}}}{ The percentage of organisms a gene must be in
#'     to be considered as part of the coregenome. \code{core_level = 95} by default.
#'     Can't be set above 100, and below 85 raises a warning.}
#'     \item{\bold{\code{$core_genes}}}{ Like \code{genes}, but only showing core genes.}
#'     \item{\bold{\code{$core_clusters}}}{ Like \code{$clusters}, but only showing core
#'     clusters.}
#'     \item{\bold{\code{$core_sequences}}}{ Like \code{$sequences}, but only showing core
#'     sequences.}
#'     \item{\bold{\code{$cloud_genes}}}{ Like \code{genes}, but only showing cloud genes.
#'     These are defined as those clusters which contain a single gene (singletons), plus
#'     those which have more than one but its organisms are probably clonal due to identical
#'     general gene content. Colloquially defined as strain-specific genes.}
#'     \item{\bold{\code{$cloud_clusters}}}{ Like \code{$clusters}, but only showing cloud
#'     clusters as defined above.}
#'     \item{\bold{\code{$cloud_sequences}}}{ Like \code{$sequences}, but only showing cloud
#'     sequences as defined above.}
#'     \item{\bold{\code{$shell_genes}}}{ Like \code{genes}, but only showing shell genes.
#'     These are defined as those clusters than don't belong nethier to the core genome,
#'     nor to cloud genome. Colloquially defined as genes that are present in some but not
#'     all strains, and that aren't strain-specific.}
#'     \item{\bold{\code{$shell_clusters}}}{ Like \code{$clusters}, but only showing shell
#'     clusters, as defined above.}
#'     \item{\bold{\code{$shell_sequences}}}{ Like \code{$sequences}, but only showing shell
#'     sequences, as defined above.}
#'     \item{\bold{\code{$summary_stats}}}{ A \code{\link[S4Vectors:DataFrame-class]{DataFrame}} with
#'     information about the number of core, shell, and cloud clusters, as well as the
#'     total number of clusters.}
#'     \item{\bold{\code{$random_seed}}}{ The last \code{.Random.seed}. Used for
#'     reproducibility purposes only.}
#'     \item{\bold{\code{$dropped}}}{ A \code{character} vector with dropped organism
#'     names, and organism number identifier as \code{names()}}
#'  }
#' }

#' @section Methods:
#' \describe{
#'     Above is a comprehensive description of all the methods provided by the object.
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{$add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{map}}: \code{character} identifiying the metadata to map. Can
#'                     be one of \code{"org"}, \code{"group"}, or \code{"gid"}.}
#'                     \item{\bold{\code{df}}: \code{data.frame} or \code{DataFrame} with the metadata to
#'                     add. For ethier case, a column named as \code{"map"} must exists, which should
#'                     contain identifiers for each element. In the case of adding gene (\code{gid})
#'                     metadata,each gene should be referenced by the name of the organism and the name
#'                     of the gene as provided in the \code{"data"} data.frame, separated by the
#'                     \code{"sep"} argument.}
#'                 }
#'             }
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Drop an organism}}{
#'             \bold{Description:}{
#'                \itemize{
#'                    Drop an organism from the dataset. This method allows to hide an organism from
#'                    the real dataset, ignoring it in downstream analyses. All the fields and
#'                    methods will behave as it doesn't exist. For instance, if you decide to drop
#'                    organism 1, the \code{$pan_matrix} field (see below) would not show it when
#'                    called.
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{$drop(x)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
#'                 \itemize{
#'                     \item{\bold{\code{x}}:  \code{character} or \code{numeric}. The name of the
#'                     organism wanted to be dropped, or its numeric id as returned in
#'                     \code{$organism} field (see below).}
#'                 }
#'             }
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with \code{x} dropped. It isn't necessary
#'                     to assign the function call to a new object, nor to re-write it as R6 objects
#'                     are mutable.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#'
#'     \item{\bold{Add metadata}}{
#'             \bold{Description:}{
#'                \itemize{
#'                   Add metadata to the object. You can add metadata to each organism, to each
#'                   group of orthologous, or to each gene. Elements with missing data should be filled
#'                   by \code{NA} (dimensions of the provided data.frame must be coherent with object
#'                   data).
#'                }
#'             }
#'
#'             \bold{Usage:}{
#'                \itemize{
#'                   \code{add_metadata(map = 'org', df)}
#'                }
#'             }
#'
#'             \bold{Arguments:}{
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
#'             \bold{Return:}{
#'                 \itemize{
#'                     \code{self} invisibly, but with additional metadata.
#'                 }
#'             }
#'     }
#' }
#'
#' @export
pagoo <- function(data, org_meta, cluster_meta, sequences, core_level = 95, sep = "__", verbose = TRUE){

  if (missing(sequences)){

    PgR6M$new(data = data,
              org_meta = org_meta,
              cluster_meta = cluster_meta,
              core_level = core_level,
              sep = sep,
              verbose = verbose)

  } else{

    PgR6M$new(data = data,
              org_meta = org_meta,
              cluster_meta = cluster_meta,
              sequences = sequences,
              core_level = core_level,
              sep = sep,
              verbose = verbose)
  }
}







#' @export
load_pangenomeRDS = function(file, ...){

  args <- readRDS(file)
  atrs <- attributes(args)
  dots <- list(...)

  if (!"package" %in% names(atrs)) stop("Not recognized rds file.")
  if (atrs$package != "pagoo") stop("Not recognized rds file.")

  dropped <- args$dropped
  args$dropped <- NULL

  sep1 <- args$sep
  sep2 <- dots$sep
  if (!is.null(sep2)){
    if (sep2 != sep1) warning('"sep" argument provided is not the same as the original object. Overwriting.')
    args$sep <- sep2
  }

  core_level1 <- args$core_level
  core_level2 <- dots$core_level
  if (!is.null(core_level2)){
    if (core_level2 != core_level1) warning('"core_level" argument provided is not the same as the original object. Overwriting.')
    args$core_level <- core_level2
  }

  p <- do.call("pagoo", args)

  if (!is.null(dropped)){
    message(paste("Dropping the following organisms:", paste0(dropped, collapse = " "), collapse = " "))
    p$drop(dropped)
  }

  return(p)
}
