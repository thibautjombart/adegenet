context("Accessor tests")


data("microbov")
strata(microbov) <- data.frame(other(microbov))

test_that("individual accessors work as expected", {
  skip_on_cran()
  expect_that(nInd(microbov), equals(704))
  expect_that(indNames(microbov), is_equivalent_to(microbov@ind.names))
  indNames(microbov)[1] <- "replacement"
  expect_that(indNames(microbov)[1], equals("replacement"))
})

test_that("population accessors work for genind objects", {
  skip_on_cran()
  expect_that(nPop(microbov), equals(15))
  expect_that(popNames(microbov), is_equivalent_to(microbov@pop.names))
  expect_that(popNames(microbov), is_equivalent_to(levels(pop(microbov))))
  popNames(microbov)[1] <- "replacement"
  expect_that(popNames(microbov)[1], equals("replacement"))
  expect_that(unique(head(pop(microbov))), is_equivalent_to(factor("replacement")))
})

test_that("population accessors work for genpop objects", {
  skip_on_cran()
  micpop <- genind2genpop(microbov, quiet = TRUE)
  expect_that(nPop(micpop), equals(15))
  expect_that(popNames(micpop), is_equivalent_to(micpop@pop.names))
  expect_that(popNames(micpop), is_equivalent_to(rownames(micpop@tab)))
  popNames(micpop)[1] <- "replacement"
  expect_that(popNames(micpop)[1], equals("replacement"))
  expect_that(rownames(micpop@tab)[1], equals("replacement"))
})

test_that("locus accessors work for genind objects", {
  skip_on_cran()
  expect_that(nLoc(microbov), equals(30))
  expect_that(locNames(microbov), is_equivalent_to(microbov@loc.names))
  locNames(microbov)[1] <- "replacement"
  expect_that(locNames(microbov)[1], equals("replacement"))
})

test_that("locus accessors work for genpop objects", {
  skip_on_cran()
  micpop <- genind2genpop(microbov, quiet = TRUE)
  expect_that(nLoc(micpop), equals(30))
  expect_that(locNames(micpop), is_equivalent_to(micpop@loc.names))
  locNames(micpop)[1] <- "replacement"
  expect_that(locNames(micpop)[1], equals("replacement"))
})