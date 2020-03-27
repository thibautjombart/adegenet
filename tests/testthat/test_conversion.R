context("data conversions")

test_that("DNAbin2genind works properly", {

  data("woodmouse", package = "ape")
  wm <- DNAbin2genind(woodmouse)
  expect_is(wm, "genind")
  expect_equal(nInd(wm), nrow(woodmouse))
   
})
