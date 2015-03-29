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
  micnames <- names(microbov)
  
  micnames <- micnames[!micnames %in% c("other", "call", "all.names")]
  for (i in micnames){
    expect_that(slot(microbov, i), is_equivalent_to(slot(retemp, i)))
  }
  for (i in names(microbov@all.names)){
    expect_that(microbov@all.names[[i]], is_equivalent_to(sort(retemp@all.names[[i]])))
  }
})