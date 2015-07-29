context("subset tests")

data("nancycats", package = "adegenet")

test_that("subsetters work for genind objects", {
  skip_on_cran()
  expect_equivalent(nInd(nancycats[1:10]), 10)
  expect_equivalent(nLoc(nancycats[loc = locNames(nancycats)[1]]), 1)
  expect_equivalent(nLoc(nancycats[loc = 1]), 1)
  expect_equivalent(nLoc(nancycats[loc = -1]), 8)
})

test_that("subsetters work for genind objects", {
  skip_on_cran()
  nanpop <- genind2genpop(nancycats, quiet = TRUE)
  expect_equivalent(nPop(nanpop[1:10]), 10)
  expect_equivalent(nLoc(nanpop[loc = locNames(nanpop)[1]]), 1)
  expect_equivalent(nLoc(nanpop[loc = 1]), 1)
  expect_equivalent(nLoc(nanpop[loc = -1]), 8)
})