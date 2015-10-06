context("Summary methods")

test_that("Diploid summaries work", {
  skip_on_cran()
  data("nancycats", package = "adegenet")
  nansum <- summary(nancycats, verbose = FALSE)
  expect_that(nansum, is_a("genindSummary"))
  expect_that(nansum$n, equals(nInd(nancycats)))
  expect_that(length(nansum$n.by.pop), equals(length(popNames(nancycats))))
  expect_that(length(nansum$pop.n.all), equals(length(popNames(nancycats))))
  expect_that(sample(nansum$Hobs, 1), is_more_than(0))
  expect_that(sample(nansum$Hexp, 1), is_more_than(0))
})
