context('Object creation, Data consistency, and Input/Output.')

## Raw Input

tgz <- system.file('extdata', 'toy_data.tar.gz', package = 'pagoo')
untar(tarfile = tgz, exdir = tempdir()) # Decompress example dataset
files <- list.files(path = tempdir(), full.names = TRUE, pattern = 'tsv$|fasta$')

data_file <- grep("case_df.tsv", files, value = TRUE)
data <- read.table(data_file, header = TRUE, sep = '\t', quote = '')

orgs_file <- grep("case_orgs_meta.tsv", files, value = TRUE)
orgs_meta <- read.table(orgs_file, header = TRUE, sep = '\t', quote = '')

clust_file <- grep("case_clusters_meta.tsv", files, value = TRUE)
clust_meta <- read.table(clust_file, header = TRUE, sep = '\t', quote = '')

fasta_files <- grep("[.]fasta", files, value = TRUE)
names(fasta_files) <- sub('[.]fasta', '', basename(fasta_files))


test_that('input data as "data.frame" works', {
  p <- pagoo(data = data)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
})

library(S4Vectors)
test_that('input data as "DataFrame" (S4Vectors) works', {
  DataF <- DataFrame(data)
  p <- pagoo(data = DataF)
  expect_type(p, "environment")
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
})

test_that('input organism metadata as "data.frame" works', {
  p <- pagoo(data, org_meta = orgs_meta)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
})

test_that('input organism metadata as "DataFrame" (S4Vectors) works', {
  Orgs_Meta <- DataFrame(orgs_meta)
  p <- pagoo(data, org_meta = Orgs_Meta)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
})

test_that('input cluster metadata as "data.frame" works', {
  p <- pagoo(data, cluster_meta = clust_meta)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
})

test_that('input cluster metadata as "DataFrame" (S4Vectors) works', {
  Clust_Meta <- DataFrame(clust_meta)
  p <- pagoo(data, cluster_meta = Clust_Meta)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
})

library(Biostrings)
sq <- lapply(fasta_files, readDNAStringSet)

test_that('input sequences as a list of "DNAStringSet" (Biostrings) works', {
  p <- pagoo(data, sequences = sq)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
  expect_is(p, "PgR6MS")
})

test_that('input sequences as a list of "character" (Biostrings) works', {
  sq_char <- lapply(sq, as.character)
  p <- pagoo(data, sequences = sq_char)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
  expect_is(p, "PgR6MS")
})

test_that('input sequences as "DNAStringSetList" (Biostrings) works', {
  sq_DSSL <- DNAStringSetList(sq)
  p <- pagoo(data, sequences = sq_DSSL)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
  expect_is(p, "PgR6MS")
})

# Data class and consistency
test_that('Class of data fields (active bindings) are the expected', {
  p <- pagoo(data, org_meta = orgs_meta, cluster_meta = clust_meta, sequences = sq)
  expect_is(p$genes, "SplitDataFrameList")
  expect_is(p$clusters, "DataFrame")
  expect_is(p$organisms, "DataFrame")
  expect_is(p$sequences, "DNAStringSetList")
  expect_is(p$pan_matrix, "matrix")
})

test_that('Data fields (active bindings) are consistent', {
  p <- pagoo(data, org_meta = orgs_meta, cluster_meta = clust_meta, sequences = sq)
  expect_equal(length(p$genes), nrow(p$clusters))
  expect_equal(nrow(p$clusters), ncol(p$pan_matrix))
  expect_equal(ncol(p$pan_matrix), length(p$sequences))
})

## RDS input/output
test_that('$save_pangenomeRDS() and load_pangenome() works', {
  p <- pagoo(data, org_meta = orgs_meta, cluster_meta = clust_meta, sequences = sq)
  out_rds <- paste0(tempdir(), '/pangenome.RDS')
  p$save_pangenomeRDS(out_rds)
  expect(file.exists(out_rds), failure_message = 'No RDS found.')

  # Check raw rds structure
  x <- readRDS(out_rds)
  expect_type(x, 'list')
  expect_named(x, c("data", "sep", "core_level", "cluster_meta", "org_meta", "sequences"),
               ignore.order = TRUE)
  expect_is(x$data, 'data.frame')
  expect_is(x$cluster_meta, 'data.frame')
  expect_is(x$org_meta, 'data.frame')
  expect_is(x$sequences, 'list')
  expect_length(x$sequences, 5)
  expect(unique(sapply(x$sequences, class)) == 'character',
         failure_message = "$sequences is not a list of character vectors")
  expect_is(x$sep, 'character')
  expect_is(x$core_level, 'numeric')
  rm(x)

  # load rds as pangenome
  p2 <- load_pangenomeRDS(out_rds)
  expect_is(p, "R6")
  expect_is(p, "PgR6")
  expect_is(p, "PgR6M")
  expect_is(p, "PgR6MS")

  # Check consistency between objects
  expect_equal(p$genes, p2$genes)
  expect_equal(p$organisms, p2$organisms)
  expect_equal(p$clusters, p2$clusters)
  expect_equivalent(p$sequences, p2$sequences)
  rm(p2)
  file.remove(out_rds)
})

## Input from third party pangenome reconstruction software
#   MISSING

file.remove(files)