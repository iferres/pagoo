

#' @name panaroo_2_pagoo
#' @title Read panaroos's output into a pagoo's R6 class object
#' @description This function handle conversion of \href{https://gtonkinhill.github.io/panaroo/#/}{panaroo}'s output files into
#' a pagoo R6 class object. It takes the "gene_presence_absence.csv" file and
#' (optionally but recommended) gff input file paths, and returns an object of
#' class \code{\link[pagoo]{PgR6MS}} (or \code{\link[pagoo]{PgR6M}} if left
#' empty the \code{gffs} argument). Panaroo identifies some genes with unusual
#' lengths tagging them with 'stop', 'length', or 'refound' labels. In the
#' current version, this function discards those genes.
#'
#' @param gene_presence_absence_csv \code{character}, path to the
#' "gene_presence_absence.csv" file. (Do not confuse with the file with the
#' same name but with \code{.Rtab} extension).
#' @param gffs A \code{character} vector with paths to original gff files used
#' as roary's input. Typically the return value of \code{list.files()} function.
#' This parameter is optional but highly recommended if you want to manipulate
#' sequences.
#' @param sep \code{character}. Default: \code{"__"} (two underscores). See
#' \link[pagoo]{PgR6MS} for a more detail argument description.
#' @return A pagoo's R6 class object. Either \link[pagoo]{PgR6M}, if \code{gffs}
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
#' @importFrom magrittr `%>%`
#' @export
panaroo_2_pagoo <- function(gene_presence_absence_csv, gffs, sep = '__'){

  message('Reading csv file (panaroo).')
  df <- read.csv(gene_presence_absence_csv,
                 header = TRUE,
                 sep = ',',
                 stringsAsFactors = FALSE,
                 check.names = FALSE,
                 colClasses = "character")

  message('Processing csv file.')

  infocols <- match(c("Gene", "Non-unique Gene name", "Annotation"), colnames(df), nomatch = 0)

  df[, -infocols] <- lapply(df[, -infocols], strsplit, ";")


  # Identify weird genes:
  dfrm <- apply(df[, -infocols], 2,
                # Look for weird genes:
                function(x){
                  grep("(_refound|_stop|_len)_*", x)
                }
  ) %>% {
    # Create a data.frame with the coordinates of the match in the original
    # data.frame:
    data.frame(
      columns = mapply(rep, names(.), lengths(.)) %>% unlist(use.names = FALSE),
      rows = unlist(., use.names = FALSE)
    )
  } %>%
    # Add the index of the match (It could be a paralogue):
    {
      idx <- mapply(function(row, column, x){
        `[[`(x, row, column) %>% grep("(_refound|_stop|_len)_*", .)
      }, row = .$rows, column = .$columns,
      MoreArgs = list(x=df[, -infocols]))
      .$index <- idx
      .
    }

  # Remove those weird genes from the data.frame:
  if (nrow(dfrm)) paste0("Removing ", nrow(dfrm), " genes tagged as 'refound', 'stop', and/or 'length' by panaroo.") %>%
    message()
  for (r in seq_len(nrow(dfrm))){
    COL <- dfrm$columns[[r]]
    ROW <- dfrm$rows[[r]]
    INDEX <- dfrm$index[[r]]
    df[[COL]][[ROW]] <- df[[COL]][[ROW]][-INDEX]
  }

  # remove clusters without genes (seems to be a panaroo's bug)
  emptyclus <- which(rowSums(do.call(cbind, lapply(df[,- infocols], lengths))) == 0)
  if (length(emptyclus)){
    df <- df[-emptyclus, ,drop=FALSE]
  }

  cluster_meta <- df[, c('Gene', 'Annotation')]
  colnames(cluster_meta) <- c('cluster', 'Annotation')

  dims <- dim(df)
  df <- df[,-infocols]
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
      sq
    })

    # Get gene information from gffs
    mcls <- lapply(seqs, mcols)
    mcls <- lapply(names(mcls), function(x){
      mcls[[x]]$org <- x
      mcls[[x]]
    })

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


