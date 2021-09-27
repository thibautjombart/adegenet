test_that("genind2genpop works with missing loci", {
  data(nancycats)
  p17 <- nancycats[pop = 17, loc = 3:4]
  p   <- genind2genpop(p17)
  expect_s4_class(p, "genpop")
  ptab <- tab(p, freq = TRUE)
  expect_equal(nrow(ptab), 1L)
  expect_equal(ncol(ptab), sum(nAll(p17)))
  expect_false(any(is.na(ptab)))

})
