context("haploGen tests")

test_that("haploGen actually works", {
  expect_is(haploGen(seq.length = 30, geo.sim = TRUE), "haploGen")
})