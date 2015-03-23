context("Genlight construction")

test_that("Genlight objects can be created predictably", {
  expect_warning(a <- new("genlight", list(c(1,0,1), c(0,0,1,0)) ))
  expect_warning(b <- new("genlight", list(c(1,0,1,1,1,1), c(1,0)) ))
  locNames(a) <- letters[1:4]
  locNames(b) <- 1:6
  c <- cbind(a, b)
  cbound <- cbind(as.matrix(a), as.matrix(b))
  rbound <- rbind(as.matrix(a),as.matrix(a))
  expect_identical(as.matrix(c), cbound) 
  expect_identical(as.matrix(rbind(a, a)), rbound)
})


test_that("subsetting with/without @other works", {
  x <- new("genlight", list(a=1,b=0,c=1), other=list(1:3, letters, data.frame(2:4)))
  pop(x) <- c("pop1","pop1", "pop2")
  expect_that(x[1:2, ]@other[[1]], equals(1:2))
  expect_that(x[1:2, ]@other[[2]], equals(letters))
  expect_that(x[1:2, ]@other[[3]], equals(x@other[[3]][1:2, , drop = FALSE]))
})