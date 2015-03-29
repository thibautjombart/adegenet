context("Genlight construction")

test_that("Genlight objects can be created predictably", {
  skip_on_cran()
  expect_warning(a <- new("genlight", list(c(1,0,1), c(0,0,1,0)), parallel = FALSE ))
  expect_warning(b <- new("genlight", list(c(1,0,1,1,1,1), c(1,0)), parallel = FALSE ))
  locNames(a) <- letters[1:4]
  locNames(b) <- 1:6
  c <- cbind(a, b)
  cbound <- cbind(as.matrix(a), as.matrix(b))
  rbound <- rbind(as.matrix(a),as.matrix(a))
  expect_identical(as.matrix(c), cbound) 
  expect_identical(as.matrix(rbind(a, a)), rbound)
})


x <- new("genlight", list(a=1,b=0,c=1), other=list(1:3, letters, data.frame(2:4)), parallel = FALSE)
pop(x) <- c("pop1","pop1", "pop2")

test_that("subsetting with/without @other works", {
  skip_on_cran()  
  expect_that(x[1:2, ]@other[[1]], equals(1:2))
  expect_that(x[1:2, ]@other[[2]], equals(letters))
  expect_that(x[1:2, ]@other[[3]], equals(x@other[[3]][1:2, , drop = FALSE]))
})

test_that("population accessors work", {
  skip_on_cran()
  expect_that(nPop(x), equals(2))
  expect_that(pop(x), is_equivalent_to(factor(c("pop1","pop1", "pop2"))))
  expect_that(popNames(x), equals(levels(pop(x))))
  popNames(x)[1] <- "replacement"
  expect_that(popNames(x), equals(c("replacement", "pop2")))
  expect_error(popNames(x) <- NULL)
  expect_error(popNames(x)[2] <- NA)
})