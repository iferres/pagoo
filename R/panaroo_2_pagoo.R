

#' @name panaroo_2_pagoo
#' @title Read panaroos's output into a pagoo's R6 class object
#' @description This function handle conversion of \href{https://gtonkinhill.github.io/panaroo/#/}{panaroo}'s output files into
#' a pagoo R6 class object. It takes the "gene_presence_absence.csv" file and
#' (optionally but recommended) gff input file paths, and returns an object of
#' class \code{\link[pagoo]{PgR6MS}} (or \code{\link[pagoo]{PgR6M}} if left
#' empty the \code{gffs} argument).
#' @param gene_presence_absence_csv \code{character}, path to the
#' "gene_presence_absence.csv" file. (Do not confuse with the file with the
#' same name but with \code{.Rtab} extension).
#' @param gffs A \code{character} vector with paths to original gff files used
#' as roary's input. Typically the return value of \code{list.files()} function.
#' This parameter is optional but highly recommended if you want to manipulate
#' sequences.
#' @param sep \code{character}. Default: \code{"__"} (two underscores). See
#' \link[pagoo]{PgR6MS} for a more detail argument description.
#' @return A pagoo's R6 class object. Ethier \link[pagoo]{PgR6M}, if \code{gffs}
#' argument is left empty, or \link[pagoo]{PgR6MS} if path to gff files is
#' provided.
#' @references
#' Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. Producing polished
#' prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21, 180 (2020).
#' https://doi.org/10.1186/s13059-020-02090-4
#' @examples
#' \dontrun{
#' gffs <- list.files(path = "path/to/gffs/",
#'                    pattern = "[.]gff$",
#'                    full.names = TRUE)
#' gpa_csv <- "path/to/gene_presence_absence.csv"
#'
#' library(pagoo)
#' pg <- panaroo_2_pagoo(gene_presence_absence_csv = gpa_csv,
#'                     gffs = gffs)
#' }
#' @importFrom utils read.csv
#' @importFrom Biostrings DNAStringSetList
#' @importFrom S4Vectors mcols
#' @importFrom stats setNames
#' @export
panaroo_2_pagoo <- function(gene_presence_absence_csv, gffs, sep = '__'){

  message('Reading csv file (panaroo).')
  df <- read.csv(gene_presence_absence_csv,
                 header = TRUE,
                 sep = ',',
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

  message('Processing csv file.')

  df[, -c(1:3)] <- lapply(df[, -c(1:3)], strsplit, ";")

  grstp <- lapply(df, function(x) grep("_stop$", x))
  if (length(unlist(grstp))) {
    warning("Removing refound genes with stop codon (tagged with '_stop')", immediate. = TRUE)
    idx <- mapply(function(X, Y, Z) {
      data.frame(COL = Z, ROW=Y, INDEX=sapply(Y, function(x) grep("_stop$", X[x])))
    },
    X=df,
    Y=grstp,
    Z=mapply(rep, names(grstp), lengths(grstp)))
    idx <- do.call(rbind, idx)
    for (r in seq_len(dim(idx)[1])){
      COL <- idx[r, 1]
      ROW <- idx[r, 2]
      INDEX <- idx[r, 3]
      df[[COL]][[ROW]] <- df[[COL]][[ROW]][-INDEX]
    }
  }

  grlen <- lapply(df, function(x) grep("_len$", x))
  if (length(unlist(grlen))){
    warning("Removing genes with unusual length (tagged with '_len')", immediate. = TRUE)
    idx <- mapply(function(X, Y, Z) {
      data.frame(COL = Z, ROW=Y, INDEX=sapply(Y, function(x) grep("_len$", X[x])))
    },
    X=df,
    Y=grlen,
    Z=mapply(rep, names(grlen), lengths(grlen)))
    idx <- do.call(rbind, idx)
    for (r in seq_len(dim(idx)[1])){
      COL <- idx[r, 1]
      ROW <- idx[r, 2]
      INDEX <- idx[r, 3]
      df[[COL]][[ROW]] <- df[[COL]][[ROW]][-INDEX]
    }
  }


  grref <- lapply(df, function(x) grep("_refound_", x))
  if (length(unlist(grref))){
    warning("Removing refound genes (tagged with '_refound_')", immediate. = TRUE)
    idx <- mapply(function(X, Y, Z) {
      data.frame(COL = Z, ROW=Y, INDEX=sapply(Y, function(x) grep("_refound_", X[x])))
    },
    X=df,
    Y=grref,
    Z=mapply(rep, names(grref), lengths(grref)))
    idx <- do.call(rbind, idx)
    for (r in seq_len(dim(idx)[1])){
      COL <- idx[r, 1]
      ROW <- idx[r, 2]
      INDEX <- idx[r, 3]
      df[[COL]][[ROW]] <- df[[COL]][[ROW]][-INDEX]
    }
  }




  # message("Cleaning merged clusters.")
  # df$Gene <- sapply(strsplit(df$Gene, "~~~"), function(x){
  #   nch <- nchar(x)
  #   x[which.min(nch)]
  # })

  cluster_meta <- df[, c('Gene', 'Annotation')]
  colnames(cluster_meta) <- c('cluster', 'Annotation')

  dims <- dim(df)
  df <- df[, 4:dims[2]]
  lp <- lapply(df, setNames, cluster_meta$cluster)

  lths <- lapply(lp, lengths)
  sums <- sapply(lths, sum)
  mm <- cbind(
    unlist(lp, use.names = FALSE),
    unlist(lapply(lths, function(x){
      unlist(rep(names(x), x), use.names = FALSE)
    }), use.names = FALSE),
    rep(names(sums), sums)
  )
  mm <- as.data.frame(mm)
  colnames(mm) <- c('gene', 'cluster', 'org')


  if(missing(gffs)){

    message('Loading PgR6M class object.')
    pg <- PgR6M$new(data = mm,
                    cluster_meta = cluster_meta,
                    sep = sep)

  }else{

    names(gffs) <- sub('[.]gff$', '', basename(gffs))
    seqs <- lapply(gffs, function(x){
      message(paste('Reading gff file', x))
      sq <- read_gff(x)
      #There are some tRNAs genes that are not listed in
      # the gene_presence_absence.csv file
      sq <- sq[which(names(sq) %in% mm$gene)]
    })

    # Get gene information from gffs
    mcls <- lapply(seqs, mcols)
    mcls <- lapply(names(mcls), function(x){
      mcls[[x]]$org <- x
      mcls[[x]]
    })
    # names(mcls) <- names(seqs)

    # Just common columns
    # cols <- Reduce(intersect, lapply(mcls, colnames))

    # Selected columns
    cols <- c('seqid', 'type', 'start', 'end', 'strand', 'product', 'org', 'locus_tag')
    mcls <- lapply(mcls, function(x) x[ , cols])

    # Match !
    mcls <- do.call(rbind, mcls)
    ma <- match(paste(mm$org, mm$gene, sep=sep), paste(mcls$org, mcls$locus_tag, sep = sep))

    mm <- cbind(mm, mcls[ma, c('seqid', 'type', 'start', 'end', 'strand', 'product')])

    message('Loading PgR6MS class object.')
    pg <- PgR6MS$new(data = mm,
                     cluster_meta = cluster_meta,
                     sequences = DNAStringSetList(seqs),
                     sep = sep)

  }

  message('Done.')
  return(pg)
}
