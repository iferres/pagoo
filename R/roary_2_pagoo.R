

#' @name roary_2_pagoo
#' @title Read roary's output into a pagoo's R6 class object
#' @description This function handle conversion of \href{https://sanger-pathogens.github.io/Roary/}{roary}'s output files into
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
#' @param paralog_sep \code{character}. A gene separator for cases where the
#' clusters have in-paralogs. (Default: "\\t" - tab).
#' @return A pagoo's R6 class object. Either \link[pagoo]{PgR6M}, if \code{gffs}
#' argument is left empty, or \link[pagoo]{PgR6MS} if path to gff files is
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
#' @importFrom Biostrings DNAStringSetList
#' @importFrom S4Vectors mcols
#' @importFrom stats setNames
#' @export
roary_2_pagoo <- function(gene_presence_absence_csv, gffs, sep = '__', paralog_sep = "\t"){

  message('Reading csv file (roary).')
  df <- read.csv(gene_presence_absence_csv,
                 header = TRUE,
                 sep = ',',
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

  message('Processing csv file.')
  cluster_meta <- df[, c('Gene', 'Annotation')]
  colnames(cluster_meta) <- c('cluster', 'Annotation')

  dims <- dim(df)
  df <- df[, 15:dims[2]]
  lp <- lapply(df, strsplit, paralog_sep)
  lp <- lapply(lp, setNames, cluster_meta$cluster)

  lths <- lapply(lp, sapply, length)
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

    names(gffs) <- sub('[.]gff3*$', '', basename(gffs))
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


# gffs <- list.files('../../pewit2_dataset/', pattern = 'gff$', full.names = T)


#' @importFrom Biostrings DNAStringSet subseq reverseComplement width xscat
#' @importFrom S4Vectors mcols<- do.call
read_gff <- function(in_gff){

  rl <- readLines(in_gff)
  gp_cm <- grep('^#', rl)
  gp_fa <- grep("^##FASTA", rl)

  co_fa <- (gp_fa + 1L):length(rl)
  rl_fa <- rl[co_fa]
  rl_tb <- rl[-c(co_fa, gp_cm)]

  # table
  m <- strsplit(rl_tb, '\t')
  df <- as.data.frame(do.call(rbind, m), stringsAsFactors = FALSE)
  df <- df[df$V2 != "prokka",]
  m2 <- strsplit(df$V9, ';')
  lp <- lapply(m2, function(x){
    spl <- strsplit(x, '=')
    val <- vapply(spl, '[', 2, FUN.VALUE = NA_character_)
    nam <- vapply(spl, '[', 1, FUN.VALUE = NA_character_)
    setNames(val, nam)
  })
  retain <- vapply(lp, function(x) "ID" %in% names(x), FUN.VALUE = NA)
  lp <- lp[retain]
  una <- unique(names(unlist(lp, use.names = T)))
  df2 <- as.data.frame(do.call(rbind, lapply(lp, function(x) {setNames(x[una], una)})), stringsAsFactors = FALSE)
  mcls <- cbind(df[retain, -9], df2)
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
  # dna <- DNAStringSet(paste0(rl_fa[ toapl[, 1]:toapl[, 2] ],collapse=''))

  dna <- DNAStringSet(apply(toapl, 1, function(x){
    paste0(rl_fa[x[1]:x[2]], collapse = '')
  }))
  names(dna) <- sub('>', '', rl_fa[he])

  # The following is to handle some special cases with bakta annotations:
  # Sometimes bakta annotates cds which starts at the end of the contigs, but
  # finish at the beginning (contigs can be circular). In these cases, bakta
  # reports "end" coordinates as a number larger than the length of the contig,
  # so if not correctly handled, the parser throws an error.
  brks <- mcls$end > width(dna[mcls$seqid])
  if (any(brks)) {
    st <- mcls$start
    en <- mcls$end
    st[brks] <- 1L
    en[brks] <- 0L
    sequences <- subseq(dna[mcls$seqid], start = st, end = en)
    sequences[brks] <- mapply(function(st, en, seqid){
      st1 <- st
      en1 <- width(dna[seqid])
      st2 <- 1L
      en2 <- en - en1
      sseq <- subseq(rep(dna[seqid], 2), c(st1, st2), c(en1, en2))
      do.call(xscat, sseq)
    },
    st = mcls$start[brks],
    en = mcls$end[brks],
    seqid = mcls$seqid[brks],
    SIMPLIFY = FALSE)
  }else{
    sequences <- subseq(dna[mcls$seqid], start = mcls$start, end = mcls$end)
  }

  sequences[mcls$strand=='-'] <- reverseComplement(sequences[mcls$strand=='-'])
  names(sequences) <- mcls$ID
  mcols(sequences) <- mcls

  sequences
}


