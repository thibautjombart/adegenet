context("Test seppop")

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



x <- new("genlight", list(a=rep(1,1e3),b=rep(0,1e3),c=rep(1, 1e3)), parallel = FALSE)
pop(x) <- c("pop1","pop2", "pop1")

test_that("seppop will work for genlight objects", {
  skip_on_cran()
  plist <- seppop(x)
  expect_is(plist, "list")
  expect_equal(length(plist), nPop(x))
  expect_equivalent(names(plist), popNames(x))
})

test_that("seppop will work for genlight objects with external factor", {
  skip_on_cran()
  uniqpop <- rev(LETTERS)[1:3]
  ulist <- seppop(x, pop = uniqpop)
  expect_is(ulist, "list")
  expect_equal(length(ulist), length(uniqpop))
  expect_equivalent(names(ulist), sort(uniqpop))
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
})