context("Repool tests")

data("microbov")
strata(microbov) <- data.frame(other(microbov))

test_that("slots are equivalent", {
  skip_on_cran()
  temp   <- seppop(microbov)
  retemp <- repool(temp)
  expect_null(other(retemp))
  expect_failure(expect_identical(microbov@call, retemp@call))
  expect_equal(names(microbov), names(retemp))

  # Alleles in repooled samples are out of order.
  # This makes sure they are ordered.
  retempallnames <- lapply(retemp@all.names, sort)
  retemptab      <- retemp@tab[, colnames(microbov@tab)]

  expect_equivalent(slot(microbov, 'tab'), retemptab)
  expect_equivalent(slot(microbov, 'all.names'), retempallnames)
  expect_equivalent(slot(microbov, 'strata'), slot(retemp, 'strata'))
  expect_equivalent(slot(microbov, 'hierarchy'), slot(retemp, 'hierarchy'))
  expect_equivalent(slot(microbov, 'loc.fac'), slot(retemp, 'loc.fac'))
  expect_equivalent(slot(microbov, 'loc.n.all'), slot(retemp, 'loc.n.all'))
  expect_equivalent(slot(microbov, 'pop'), slot(retemp, 'pop'))
  expect_equivalent(slot(microbov, 'ploidy'), slot(retemp, 'ploidy'))
  expect_equivalent(slot(microbov, 'type'), slot(retemp, 'type'))
  
})
