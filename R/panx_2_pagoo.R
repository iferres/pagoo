# @importFrom reshape2 melt
# @export
# panx_2_pagoo <- function(allclusters_final_tsv, nucleotide_fnas, group_prefix = 'group', sep = '__'){
#
#   x <- readLines(allclusters_final_tsv)
#   spl <- strsplit(x, '\t')
#   spl <- setNames(spl, paste0(group_prefix, seq_len(length(x))))
#   sb <- lapply(spl, function(i) sub('[|]', sep, i))
#   df <- melt(sb)
#   df$org <- sub(paste0(sep, '.+$'), '', df$value)
#   colnames(df) <- c('gene', 'cluster', 'org')
#
#   if (missing(nucleotide_fnas)){
#
#     pg <- PgR6M$new(data = df, sep = sep)
#
#   }else{
#
#     message('Not implemented yet')
#     NULL
#
#   }
#
#
# }
