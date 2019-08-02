
# df <- list(
#   Gene=c('asdf', 'sdfg', 'dfgh', 'fghj'),
#   Annotation=c('anot1', 'anot2', 'anot1', 'anot3'),
#   org1=c('org1_01', 'org1_02;org1_03', 'org1_04', 'org1_05'),
#   org2=c('org2_01', 'org2_02', 'org2_03;org2_04', 'org2_05'),
#   org3=c('org3_01', 'org3_02;org3_03', 'org3_04', 'org3_05'),
#   org4=c('org4_01;org4_02', 'org4_03', 'org4_04', 'org4_05')
# )
# df <- as.data.frame(df)
# df[] <- lapply(df, as.character)
#
# gnam <- df$Gene
# gann <- df$Annotation
# dims <- dim(df)
# df <- df[, 3:dims[2]]
# lp <- lapply(df, strsplit, ';')
# lp <- lapply(lp, setNames, gnam)
#
# mm <- melt(lp)
# colnames(mm) <- c('gene', 'group', 'org')
# mm$Annotation <- gann[factor(mm$group, levels = gnam)]

#' @importFrom utils read.csv
#' @importFrom reshape2 melt
#' @importFrom Biostrings DNAStringSetList
#' @importFrom S4Vectors mcols
#' @export
roary_2_pagoo <- function(gene_presence_absence_csv, gffs, sep = '___'){

  df <- read.csv(gene_presence_absence_csv,
                 header = TRUE,
                 sep = ',',
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

  group_meta <- df[, c('Gene', 'Annotation')]
  colnames(group_meta) <- c('group', 'Annotation')

  dims <- dim(df)
  df <- df[, 15:dims[2]]
  lp <- lapply(df, strsplit, ';')
  lp <- lapply(lp, setNames, group_meta$group)

  mm <- melt(lp)
  colnames(mm) <- c('gene', 'group', 'org')

  if(missing(gffs)){

    pg <- PgR6M$new(DF = mm,
                    group_meta = group_meta,
                    sep = sep)

  }else{

    names(gffs) <- sub('[.]gff$', '', basename(gffs))
    seqs <- lapply(gffs, function(x){
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



    pg <- PgR6MS$new(DF = mm,
                     group_meta = group_meta,
                     sequences = DNAStringSetList(seqs),
                     sep = sep)

  }

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


