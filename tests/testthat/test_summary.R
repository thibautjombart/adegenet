context("Summary methods")

test_that("Diploid summaries work", {
  skip_on_cran()
  data("nancycats", package = "adegenet")
  nansum <- summary(nancycats)
  expect_that(nansum, is_a("list"))
  expect_that(nansum$N, equals(nInd(nancycats)))
  expect_that(length(nansum$pop.eff), equals(length(popNames(nancycats))))
  expect_that(length(nansum$pop.nall), equals(length(popNames(nancycats))))
  expect_that(sample(nansum$Hobs, 1), is_more_than(0))
  expect_that(sample(nansum$Hexp, 1), is_more_than(0))
})
