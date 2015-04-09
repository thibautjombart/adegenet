context("Test PropShared")

data("microbov", package = "adegenet")
obj <- microbov[1:5, loc = locNames(microbov)[1:2]]
psh <- propShared(obj)
expected_values <- structure(c(1, 0.5, 0.5, 0.75, 0.5, 0.5, 1, 0.75, 0.75, 0.75,
              0.5, 0.75, 1, 0.75, 1, 0.75, 0.75, 0.75, 1, 0.75, 0.5, 0.75, 1,
              0.75, 1), .Dim = c(5L, 5L), .Dimnames = list(c("AFBIBOR9503",
              "AFBIBOR9504", "AFBIBOR9505", "AFBIBOR9506", "AFBIBOR9507"),
              c("AFBIBOR9503", "AFBIBOR9504", "AFBIBOR9505", "AFBIBOR9506",
              "AFBIBOR9507")))
test_that("propShared produces expected results", {
  skip_on_cran()
  expect_that(psh, equals(expected_values))
  expect_true(all(psh <= 1))
  expect_that(ncol(psh), equals(nrow(psh)))
  expect_that(ncol(psh), equals(nInd(obj)))
})