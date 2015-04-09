context("xvalDapc test")

xvalnames <- c("Cross-Validation Results", 
  "Median and Confidence Interval for Random Chance", 
  "Mean Successful Assignment by Number of PCs of PCA", 
  "Number of PCs Achieving Highest Mean Success", 
  "Root Mean Squared Error by Number of PCs of PCA", 
  "Number of PCs Achieving Lowest MSE", 
  "DAPC")

data(sim2pop)


test_that("xvalDapc returns expected results", {
  skip_on_cran()
  xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=3, xval.plot = FALSE)
  cvr  <- xval[["Cross-Validation Results"]]
  msa  <- xval[["Mean Successful Assignment by Number of PCs of PCA"]]
  expect_that(xval, is_a("list"))
  expect_equivalent(names(xval), xvalnames)
  expect_that(xval$DAPC, is_a("dapc"))
  expect_that(nrow(cvr), equals(3 * length(msa)))
})

test_that("xvalDapc throws a warning with populations of 1 sample", {
  skip_on_cran()
  dat_pop <- as.character(pop(sim2pop))
  dat_pop[1] <- "Pop C"
  pop(sim2pop) <- dat_pop
  expect_warning(xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=3, xval.plot = FALSE))
})
