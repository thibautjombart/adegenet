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
  retempallnames <- lapply(retemp@all.names, sort)

  expect_that(slot(microbov, 'all.names'), is_equivalent_to(retempallnames))
  expect_that(slot(microbov, 'strata'), is_equivalent_to(slot(retemp, 'strata')))
  expect_that(slot(microbov, 'hierarchy'), is_equivalent_to(slot(retemp, 'hierarchy')))
  expect_that(slot(microbov, 'tab'), is_equivalent_to(slot(retemp, 'tab')))
  expect_that(slot(microbov, 'loc.names'), is_equivalent_to(slot(retemp, 'loc.names')))
  expect_that(slot(microbov, 'loc.fac'), is_equivalent_to(slot(retemp, 'loc.fac')))
  expect_that(slot(microbov, 'loc.nall'), is_equivalent_to(slot(retemp, 'loc.nall')))
  expect_that(slot(microbov, 'ind.names'), is_equivalent_to(slot(retemp, 'ind.names')))
  expect_that(slot(microbov, 'pop'), is_equivalent_to(slot(retemp, 'pop')))
  expect_that(slot(microbov, 'pop.names'), is_equivalent_to(slot(retemp, 'pop.names')))
  expect_that(slot(microbov, 'ploidy'), is_equivalent_to(slot(retemp, 'ploidy')))
  expect_that(slot(microbov, 'type'), is_equivalent_to(slot(retemp, 'type')))
  
})