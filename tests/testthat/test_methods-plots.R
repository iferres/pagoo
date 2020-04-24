context("Methods and Plots")

rds <- system.file('extdata', 'campylobacter.RDS', package = 'pagoo')
p <- load_pangenomeRDS(rds)


test_that('methods return expected objects', {
  dist <- p$dist()
  expect_is(dist, 'dist')

  flui <- p$fluidity()
  expect_is(flui, "numeric")
  expect_named(flui, expected = c("Mean", "Std"))

  rare <- p$rarefact()
  expect_is(rare, "matrix")
  expect(all(apply(rare, 1, class) == "numeric"), "Matrix class is not numeric")
  expect_equal(nrow(rare), nrow(p$organisms))

  pow <- p$pg_power_law_fit()
  expect_is(pow, "list")
  expect_length(pow, 2)
  expect_named(pow, c("formula", "params"))

  dec <- p$cg_exp_decay_fit()
  expect_is(dec, "list")
  expect_length(dec, 2)
  expect_named(dec, c("formula", "params"))

  pca <- p$pan_pca()
  expect_is(pca, "prcomp")

  core <- p$core_seqs_4_phylo()
  orgs <- as.character(p$organisms$org)
  expect_is(core, "DNAStringSetList")
  expect_length(core, length(p$core_genes))
  expect(all(sapply(core, length) == nrow(p$organisms)),
         "core_seqs_4_phylo() length elements are not equal to the number of organisms")
  expect_true(all(unlist(lapply(core, function(x) mcols(x)$org == orgs ), use.names = FALSE)),
              "core_seqs_4_phylo() organism order is not the same as .$organisms$org")
})


test_that('Plots returns gg objects', {
  expect_is(p$gg_pie(), "gg")
  expect_is(p$gg_barplot(), "gg")
  expect_is(p$gg_binmap(), "gg")
  expect_is(p$gg_curves(), "gg")
  expect_is(p$gg_dist(), "gg")
  expect_is(p$gg_pca(), "gg")
})
