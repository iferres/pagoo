

#' @name roary_2_pagoo
#' @title Read roary's output into a pagoo's R6 class object
#' @description This function handle conversion of \href{https://sanger-pathogens.github.io/Roary/}{roary}'s output files into
#' a pagoo R6 class object. It takes the "gene_presence_absence.csv" file and
#' (optionally but recommended) gff input file paths, and returns an object of
#' class \code{\link[pagoo]{PgR6MS}} (or \code{\link[pagoo]{PgR6M}} if left
#' empthy the \code{gffs} argument).
#' @param gene_presence_absence.csv \code{character}, path to the
#' "gene_presence_absence.csv" file. (Do not confuse with the file with the
#' same name but with \code{.Rtab} extension).
#' @param gffs A \code{character} vector with paths to original gff files used
#' as roary's input. Typically the return value of \code{list.files()} function.
#' This parameter is optional but highly recommended if you want to manipulate
#' sequences.
#' @param sep \code{character}. Default: \code{"__"} (two underscores). See
#' \link[pagoo]{PgR6MS} for a more detail argument description.
#' @return A pagoo's R6 class object. Ethier \link[pagoo]{PgR6M}, if \code{gffs}
#' argument is left empthy, or \link[pagoo]{PgR6MS} if path to gff files is
#' provided.
#' @references
#' Andrew J. Page, Carla A. Cummins, Martin Hunt, Vanessa K. Wong, Sandra Reuter,
#'  Matthew T. G. Holden, Maria Fookes, Daniel Falush, Jacqueline A. Keane, Julian
#'  Parkhill, "Roary: Rapid large-scale prokaryote pan genome analysis",
#'  Bioinformatics, 2015;31(22):3691-3693
#' @examples
#' \dontrun{
#' gffs <- list.files(path = "path/to/gffs/",
#'                    pattern = "[.]gff$",
#'                    full.names = TRUE)
#' gpa_csv <- "path/to/gene_presence_absence.csv"
#'
#' library(pagoo)
#' pg <- roary_2_pagoo(gene_presence_absence_csv = gpa_csv,
#'                     gffs = gffs)
#' }
#' @importFrom utils read.csv
#' @importFrom reshape2 melt
#' @importFrom Biostrings DNAStringSetList
#' @importFrom S4Vectors mcols
#' @export
roary_2_pagoo <- function(gene_presence_absence_csv, gffs, sep = '__'){

  message('Reading csv file (roary).')
  df <- read.csv(gene_presence_absence_csv,
                 header = TRUE,
                 sep = ',',
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

  message('Processing csv file.')
  group_meta <- df[, c('Gene', 'Annotation')]
  colnames(group_meta) <- c('group', 'Annotation')

  dims <- dim(df)
  df <- df[, 15:dims[2]]
  lp <- lapply(df, strsplit, ';')
  lp <- lapply(lp, setNames, group_meta$group)

  mm <- melt(lp)
  colnames(mm) <- c('gene', 'group', 'org')


  if(missing(gffs)){

    message('Loading PgR6M class object.')
    pg <- PgR6M$new(DF = mm,
                    group_meta = group_meta,
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
    pg <- PgR6MS$new(DF = mm,
                     group_meta = group_meta,
                     sequences = DNAStringSetList(seqs),
                     sep = sep)

  }

  message('Done.')
  return(pg)
}


# gffs <- list.files('../../pewit2_dataset/', pattern = 'gff$', full.names = T)


#' @importFrom Biostrings DNAStringSet subseq reverseComplement
#' @importFrom S4Vectors mcols<-
read_gff <- function(in_gff){

  rl <- readLines(in_gff)
  gp_cm <- grep('^##', rl)
  rl_fa <- rl[(rev(gp_cm)[1]+1L):length(rl)]
  rl_tb <- rl[(rev(gp_cm)[2]+1L):(rev(gp_cm)[1]-1L)]

  # table
  m <- strsplit(rl_tb, '\t')
  df <- as.data.frame(do.call(rbind, m), stringsAsFactors = FALSE)
  m2 <- strsplit(df$V9, ';')
  lp <- lapply(m2, function(x){
    spl <- strsplit(x, '=')
    val <- vapply(spl, '[', 2, FUN.VALUE = NA_character_)
    nam <- vapply(spl, '[', 1, FUN.VALUE = NA_character_)
    setNames(val, nam)
  })
  una <- unique(names(unlist(lp, use.names = T)))
  df2 <- as.data.frame(do.call(rbind, lapply(lp, '[', una)), stringsAsFactors = FALSE)
  mcls <- cbind(df[, -9], df2)
  colnames(mcls)[1:8] <- c('seqid', 'source', 'type', 'start',
                           'end', 'score', 'strand', 'phase')
  mcls$start <- as.integer(mcls$start)
  mcls$end <- as.integer(mcls$end)

  # sequence
  he <- grep('>', rl_fa)
  st <- he + 1L
  lst <- rev(st)[1]
  en <- c((he - 1)[-1], length(rl_fa))
  toapl <- cbind(st, en)
  dna <- DNAStringSet(paste0(rl_fa[ toapl[, 1]:toapl[, 2] ],collapse=''))

  # dna <- DNAStringSet(apply(toapl, 1, function(x){
  #   paste0(rl_fa[x[1]:x[2]], collapse = '')
  # }))
  names(dna) <- sub('>', '', rl_fa[he])

  sequences <- subseq(dna[mcls$seqid], start = mcls$start, end = mcls$end)
  sequences[mcls$strand=='-'] <- reverseComplement(sequences[mcls$strand=='-'])
  names(sequences) <- mcls$locus_tag
  mcols(sequences) <- mcls

  sequences
}


