
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
