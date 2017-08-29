context("Compoplot Tests")

mat <- t(apply(matrix(sample(10, 30, replace = TRUE), ncol = 3), 1, function(x) x/sum(x)))

test_that("compoplot.matrix works", {
  skip_on_cran()
  expect_silent(compoplot(mat))
})

test_that("compoplot works with a custom palette", {
  skip_on_cran()
  expect_silent(compoplot(mat, col.pal = c("cyan", "magenta", "yellow")))
})

test_that("compoplot works with a named custom palette", {
  skip_on_cran()
  expect_silent(compoplot(mat, col.pal = c(`1`="cyan", `2`="magenta", `3`="yellow")))
})

test_that("compoplot throws a warning if there are too many colors", {
  skip_on_cran()
  expect_warning(compoplot(mat, funky(4)), "populations fewer than")
})

test_that("compoplot throws a warning if there are not enough colors and does its own thing", {
  skip_on_cran()
  expect_warning(compoplot(mat, c("red", "blue")), "Using funky()")
})