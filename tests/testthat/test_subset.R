context("Subset functions")

rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
p <- load_pangenomeRDS(rds)

# Predefined subsets
test_that('predefined subsets for $genes are consistent', {
  expect_is(p$core_genes, is(p$genes))
  expect_is(p$shell_genes, is(p$genes))
  expect_is(p$cloud_genes, is(p$genes))
  expect(length(p$core_genes) +
           length(p$shell_genes) +
           length(p$cloud_genes) == length(p$genes),
         failure_message = 'Inconsistent length for predefied gene subsets.')
})

test_that('predefined subsets for $clusters are consistent', {
  expect_is(p$core_clusters, is(p$clusters))
  expect_is(p$shell_genes, is(p$clusters))
  expect_is(p$cloud_genes, is(p$clusters))
  expect(nrow(p$core_clusters) +
           nrow(p$shell_clusters) +
           nrow(p$cloud_clusters) == nrow(p$clusters),
         failure_message = 'Inconsistent nrow for predefied clusters subsets.')
})

test_that('predefined subsets for $sequences are consistent', {
  expect_is(p$core_sequences, is(p$sequences))
  expect_is(p$shell_sequences, is(p$sequences))
  expect_is(p$cloud_sequences, is(p$sequences))
  expect(length(p$core_sequences) +
           length(p$shell_sequences) +
           length(p$cloud_sequences) == length(p$sequences),
         failure_message = 'Inconsistent length for predefied sequences subsets.')
})

# Changing core_level affects subsets
test_that('Changing $core_level affects predefined subsets', {
  suma1 <- p$summary_stats
  p$core_level <- 100
  suma2 <- p$summary_stats
  expect_false(all(suma1$Number == suma2$Number))
  p$core_level <- 95 # reset
})

# [ operator: Vector notation
test_that('vector/list notation for subset works', {
  # Check raw (numeric)
  expect_is(p[1], "environment")
  # Check subsetting by cluster name (character)
  clus_1_char <- p[1]$clusters$cluster # character
  expect_is(p[clus_1_char], "environment")
  # Are both equivalent?
  expect_equal(p[1], p[clus_1_char])
  # Expected class of subquery?
  expect_is(p[1]$genes, is(p$genes))
  expect_is(p[1]$clusters, is(p$clusters))
  expect_is(p[1]$organisms, is(p$organisms))
  expect_is(p[1]$sequences, is(p$sequences))
  # Mappings are ok?
  expect_setequal(
    as.character(p[1]$genes[[1]]$cluster), p[1]$clusters$cluster
  )
  expect_setequal(
    as.character(p[1]$genes[[1]]$org), as.character(p[1]$organisms$org)
  )
  expect_setequal(
    as.character(p[1]$genes[[1]]$gid) , names(p[1]$sequences[[1]])
  )
})

# [ operator: Matrix notation
test_that('matrix notation for subset works', {
  expect_is(p[1, 1], "environment")
  clus_1_char <- p[1,1]$clusters$cluster # character
  orgs_1_char <- as.character(p[1,1]$organisms$org)
  expect_is(p[orgs_1_char, clus_1_char], "environment")
  expect_equal(p[1,1], p[orgs_1_char, clus_1_char])
  expect_is(p[1, 1]$genes, is(p$genes))
  expect_is(p[1, 1]$clusters, is(p$clusters))
  expect_is(p[1, 1]$organisms, is(p$organisms))
  expect_is(p[1, 1]$sequences, is(p$sequences))
  expect_setequal(
    as.character(p[1, 1]$genes[[1]]$cluster), p[1, 1]$clusters$cluster
  )
  expect_setequal(
    as.character(p[1, 1]$genes[[1]]$org), as.character(p[1, 1]$organisms$org)
  )
  expect_setequal(
    as.character(p[1, 1]$genes[[1]]$gid) , names(p[1, 1]$sequences[[1]])
  )
})

# Drop/Recover organisms
test_that('dropping/recovering organisms works', {
  norgs <- nrow(p$organisms)
  todrop <- p$organisms$org[c(1, 3, 5)]
  # Dropping by index
  p$drop(c(1, 3, 5))
  expect_equivalent(as.character(todrop), p$dropped)
  expect_equal(nrow(p$organisms), norgs - 3)
  expect_false(any(as.character(todrop) %in%  as.character(p$organisms$org)))
  expect_false(any(as.character(todrop) %in% as.character(unlist(p$genes, use.names = FALSE)$org)))
  p$recover( p$dropped )
  expect_equal(nrow(p$organisms), norgs)
  # Dropping by name
  p$drop( as.character(todrop) )
  expect_equivalent(as.character(todrop), p$dropped)
  expect_equal(nrow(p$organisms), norgs - 3)
  expect_false(any(as.character(todrop) %in%  as.character(p$organisms$org)))
  expect_false(any(as.character(todrop) %in% as.character(unlist(p$genes, use.names = FALSE)$org)))
  p$recover( p$dropped )
  expect_equal(nrow(p$organisms), norgs)
})


