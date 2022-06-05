
test_that("as.matrix.genind treats missing data", {
  # https://github.com/thibautjombart/adegenet/issues/331
  data(nancycats)
  nan1 <- nancycats[pop = 1]
  asis <- as.matrix(nan1)
  mean <- as.matrix(nan1, NA.method = "mean")
  zero <- as.matrix(nan1, NA.method = "zero")

  expect_true(anyNA(asis))
  expect_type(asis, "integer")

  expect_false(anyNA(mean))
  expect_type(mean, "double")

  expect_false(anyNA(zero))
  expect_type(zero, "integer")

})
