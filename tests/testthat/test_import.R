context("Import Tests")

test_that("df2genind works with haploids", {
  skip_on_cran()
  x <- matrix(sample(20), nrow = 10, ncol = 2)
  res <- df2genind(x, ploidy = 1)
  expect_that(sum(res@loc.nall), equals(20))
  expect_that(nInd(res), equals(10))
  expect_that(nLoc(res), equals(2))
  resdf <- genind2df(res)
  expect_that(as.matrix(resdf), is_equivalent_to(x))
})

test_that("df2genind makes sense for given example", {
  skip_on_cran()
  df <- data.frame(locusA=c("11","11","12","32"),
                   locusB=c(NA,"34","55","15"),
                   locusC=c("22","22","21","22"))
  row.names(df) <- .genlab("genotype",4)
  obj <- df2genind(df, ploidy=2, ncode=1)
  expect_that(nInd(obj), equals(4))
  expect_that(nLoc(obj), equals(3))
  expect_that(locNames(obj), equals(colnames(df)))
  expect_that(indNames(obj), equals(rownames(df)))
  expect_that(obj@loc.nall, is_equivalent_to(c(3, 4, 2)))
  objdf <- genind2df(obj)
  expect_that(df, is_equivalent_to(df))
})

test_that("df2genind will handle duplicate samples and loci", {
  skip_on_cran()
  x <- 
    "A B
  1/2 3/4
  5/6 4/5
  2/6 3/9"
  xdf <- read.table(text = x, header = TRUE, stringsAsFactors = FALSE)
  inds <- c("one", "one", "two")
  loci <- rep("double", 2)
  expect_warning(xgid <- df2genind(xdf, sep = "/", ind.names = inds, loc.names = loci))
  expect_that(unique(rowSums(tab(xgid))), is_equivalent_to(4))
  expect_that(genind2df(xgid, sep = "/"), is_equivalent_to(xdf))
})