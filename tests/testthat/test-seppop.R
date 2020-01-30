context("Test seppop for genind")

data(microbov)

test_that("seppop will use the internal population factor by default", {
  skip_on_cran()
  blist <- seppop(microbov)
  expect_is(blist, "list")
  expect_equal(length(blist), nPop(microbov))
  expect_equivalent(names(blist), popNames(microbov))
})

test_that("seppop will use the external population factor", {
  skip_on_cran()
  coun <- other(microbov)$coun
  clist <- seppop(microbov, pop = coun)
  expect_is(clist, "list")
  expect_equal(length(clist), nlevels(coun))
  expect_equivalent(names(clist), levels(coun))
})

test_that("seppop will use formula input", {
  skip_on_cran()
  strata(microbov) <- data.frame(other(microbov))
  slist <- seppop(microbov, pop = ~spe)
  setPop(microbov) <- ~spe
  expect_is(slist, "list")
  expect_equal(length(slist), nPop(microbov))
  expect_equivalent(names(slist), popNames(microbov))
})

test_that("seppop will throw a warning if there are missing populations", {
  skip_on_cran()
  # Create 10 ambiguous population assignments
  pop(microbov)[sample(nInd(microbov), 10)] <- NA
  
  # Missing data will throw a warning by default
  expect_warning(res <- seppop(microbov), "missing population information")
  # The resulting list will be equal to the number of populations
  expect_length(res, nPop(microbov))
  
  # keepNA does not throw a warning
  expect_failure(expect_warning(res2 <- seppop(microbov, keepNA = TRUE)))
  # The resulting list will have 1 more population
  expect_length(res2, nPop(microbov) + 1L)
  # This population will be equal to the number of missing samples
  expect_equal(nInd(res2[[length(res2)]]), 10L)
})

context("Test seppop for genight")

x <- new("genlight", list(a=rep(1,1e3),b=rep(0,1e3),c=rep(1, 1e3)), parallel = FALSE)
pop(x) <- c("pop1","pop2", "pop1")

test_that("seppop will work for genlight objects", {
  skip_on_cran()
  plist <- seppop(x)
  expect_is(plist, "list")
  expect_equal(length(plist), nPop(x))
  expect_equivalent(names(plist), popNames(x))
  expect_equivalent(vapply(plist, nInd, integer(1)), setNames(c(2, 1), popNames(x)))
})

test_that("seppop will work for genlight objects with external factor", {
  skip_on_cran()
  uniqpop <- rev(LETTERS)[1:3]
  ulist <- seppop(x, pop = uniqpop)
  expect_is(ulist, "list")
  expect_equal(length(ulist), length(uniqpop))
  expect_equivalent(names(ulist), sort(uniqpop))
  expect_equal(vapply(ulist, nInd, integer(1)), rep(1, 3))
})

test_that("seppop will work for genlight objects with formula", {
  skip_on_cran()
  uniqpop <- rev(LETTERS)[1:3]
  strata(x) <- data.frame(pop = pop(x), uni = uniqpop, all = rep("A", 3))
  setPop(x) <- ~all
  alist <- seppop(x, pop = ~all)
  expect_is(alist, "list")
  expect_equal(length(alist), nPop(x))
  expect_equivalent(names(alist), popNames(x))
  expect_equal(nInd(alist[[1]]), nInd(x))
})
