context("Repool tests")

data("microbov")
strata(microbov) <- data.frame(other(microbov))

test_that("slots are equivalent", {
  skip_on_cran()
  temp   <- seppop(microbov)
  retemp <- repool(temp)
  expect_null(other(retemp))
  expect_that(microbov@call, not(is_equivalent_to(retemp@call)))
  expect_that(names(microbov), equals(names(retemp)))

#   micnames <- names(microbov)
#   micnames <- micnames[!micnames %in% c("other", "call", "all.names")]
#   for (i in micnames){
#     x <- paste0("expect_that(slot(microbov, '", i, "'), is_equivalent_to(slot(retemp, '", i, "')))", "\n")
#     cat(x)
#   }
  # Alleles in repooled samples are out of order.
  # This makes sure they are ordered.
  retempallnames <- lapply(retemp@all.names, sort)
  retemptab      <- retemp@tab[, colnames(microbov@tab)]

  expect_equivalent(slot(microbov, 'tab'), retemptab)
  expect_equivalent(slot(microbov, 'all.names'), retempallnames)
  expect_equivalent(slot(microbov, 'strata'), slot(retemp, 'strata'))
  expect_equivalent(slot(microbov, 'hierarchy'), slot(retemp, 'hierarchy'))
  expect_equivalent(slot(microbov, 'loc.fac'), slot(retemp, 'loc.fac'))
  expect_equivalent(slot(microbov, 'loc.nall'), slot(retemp, 'loc.nall'))
  expect_equivalent(slot(microbov, 'pop'), slot(retemp, 'pop'))
  expect_equivalent(slot(microbov, 'ploidy'), slot(retemp, 'ploidy'))
  expect_equivalent(slot(microbov, 'type'), slot(retemp, 'type'))

})
