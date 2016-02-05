context("Genlight construction")

test_that("Genlight objects can be created predictably", {
  skip_on_cran()
  expect_warning(a <- new("genlight", list(c(1,0,1), c(0,0,1,0)), parallel = FALSE ))
  expect_warning(b <- new("genlight", list(c(1,0,1,1,1,1), c(1,0)), parallel = FALSE ))
  locNames(a) <- letters[1:4]
  locNames(b) <- 1:6
  c <- cbind(a, b, parallel = FALSE)
  cbound <- cbind(as.matrix(a), as.matrix(b))
  rbound <- rbind(as.matrix(a),as.matrix(a))
  expect_identical(as.matrix(c), cbound) 
  expect_identical(as.matrix(rbind(a, a, parallel = FALSE)), rbound)
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

x <- "
X13049 X13050 X13051 X13052 X13053
AA36881      2     NA      2      2      2
AA36883      2      2      2      2      2
AA36884      2      2      2      2      2
AA36802     NA      2      2      2      2
AA36803      2      2      2      2      2
AA36804      2     NA      2      2      2
AA36181      2     NA      2      2      2
AA36183      2      2      2      2      2"

xxdf <- read.table(text = x)
xx   <- new("genlight", xxdf, parallel = FALSE)
pop(xx) <- rep(LETTERS[1:2], each = 4)

test_that("missing data is properly subset with logical subscripts", {
  skip_on_cran()
  Apop <- pop(xx) == "A"
  Bpop <- pop(xx) == "B"
  expect_identical(NA.posi(xx), NA.posi(xx[]))
  expect_identical(xxdf[Apop, ], as.data.frame(xx[Apop, ]))
  expect_identical(xxdf[Bpop, ], as.data.frame(xx[Bpop, ]))
  keepers <- c(FALSE, rep(TRUE, 4))
  expect_identical(xxdf[keepers], as.data.frame(xx[, keepers]))
})

test_that("missing data is properly subset with positive subscripts", {
  skip_on_cran()
  rl <- sample(5)
  # Can subset single locus
  expect_identical(xxdf[, 1, drop = FALSE], as.data.frame(xx[, 1]))
  # Can subset range of loci
  expect_identical(xxdf[, 1:3, drop = FALSE], as.data.frame(xx[, 1:3]))
  # Can subset by position
  expect_identical(xxdf[, rl, drop = FALSE], as.data.frame(xx[, rl]))
})

test_that("missing data is properly subset with negative subscripts", {
  skip_on_cran()
  expect_identical(xxdf[, -1], as.data.frame(xx[, -1]))
  expect_identical(xxdf[, -c(1, 3)], as.data.frame(xx[, -c(1, 3)]))
})

test_that("missing data is properly subset with a character vector", {
  skip_on_cran()
  lnames <- locNames(xx)
  rl     <- sample(lnames)
  expect_identical(xxdf[, rl], as.data.frame(xx[, rl]))
  expect_identical(xxdf[, lnames[1:2]], as.data.frame(xx[, lnames[1:2]]))
})

test_that("genlight objects do not take a mixture of positive and negative subscripts", {
  skip_on_cran()
  expect_error(xx[, c(2, -1)], "subscripts.")
})


test_that("glSim does not call parallel by default", {
  skip_on_cran()
  detach("package:parallel")
  no_parallel <- sessionInfo()$basePkgs
  expect_false("parallel" %in% no_parallel)
  x <- glSim(2, n.snp.nonstruc = 10, n.snp.struc = 10, parallel = FALSE)
  check_parallel <- sessionInfo()$basePkgs
  expect_false("parallel" %in% check_parallel)
  expect_is(x, "genlight")
})