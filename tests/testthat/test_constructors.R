context("constructor tests")

test_that("genind objects can be constructed with nothing", {
  skip_on_cran()
  gen <- genind()
  expect_that(gen, is_a("genind"))
  expect_true(is.genind(gen))
})

test_that("genpop objects can be constructed with nothing", {
  skip_on_cran()
  gen <- genpop()
  expect_that(gen, is_a("genpop"))
  expect_true(is.genpop(gen))
})