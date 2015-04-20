context("subset tests")

data("nancycats", package = "adegenet")

test_that("subsetters work for genind objects", {
  skip_on_cran()
  test_that(nInd(nancycats[1:10]), equals(10))
  test_that(nLoc(nancycats[loc = locNames(nancycats)[1]]), equals(1))
})

test_that("subsetters work for genind objects", {
  skip_on_cran()
  nanpop <- genind2genpop(nancycats, quiet = TRUE)
  test_that(nPop(nanpop[1:10]), equals(1))
  test_that(nLoc(nanpop[loc = locNames(nanpop)[1]]), equals(1))
})