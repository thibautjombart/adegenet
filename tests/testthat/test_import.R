context("Import Tests")

test_that("df2genind works with haploids", {
  skip_on_cran()
  x <- matrix(sample(20), nrow = 10, ncol = 2)
  res <- df2genind(x, ploidy = 1)
  expect_that(sum(res@loc.n.all), equals(20))
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
  expect_that(obj@loc.n.all, is_equivalent_to(c(3, 4, 2)))
  objdf <- genind2df(obj)
  expect_that(df, is_equivalent_to(df))
})

test_that("df2genind handles NAs for 'numerically named' samples correctly", {
  skip_on_cran()
  
  df <- read.table(text = "
AnimalID,Samp,INRA21,AHT137,REN169D01,AHTk253
730,AX.0CCE,092 098,132 132,NA,284 286
498,AP.07P4,092 092,124 142,204 208,280 280
677,AP.088P,092 096,140 146,204 204,280 280
678,AP.088T,096 098,124 148,198 204,280 280
544,AP.07XM,096 098,134 146,198 198,280 286
533,AP.07UM,092 098,134 148,198 204,280 286", 
                   header = TRUE, sep = ",", colClasses = rep("factor", 6))
  
  obj <- df2genind(X = df[, !grepl("AnimalID|Samp", colnames(df))], ind.names = df$AnimalID,
                   sep = " ", ncode = 6)
  g <- tab(obj)
  expect_that(g["730", grepl("REN169D01", colnames(g))], 
              is_equivalent_to(c(REN169D01.204 = as.integer(NA), 
                                 REN169D01.208 = as.integer(NA), 
                                 REN169D01.198 = as.integer(NA)))
  )
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


test_that("read.X functions work as expected", {
  skip_on_cran()
  gpop <- read.genepop(system.file("files/nancycats.gen",package="adegenet"), quiet = TRUE)
  fsta <- read.fstat(system.file("files/nancycats.dat",package="adegenet"), quiet = TRUE)
  gntx <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"), quiet = TRUE)
  stru <- read.structure(system.file("files/nancycats.str",package="adegenet"),
                         onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, 
                         col.pop=2, ask=FALSE, quiet = TRUE)
  data("nancycats", package = "adegenet")
  # Making sure that the populations are all named the same. The order of the
  # isolates are mixed up within these data.
  levels(pop(gpop)) <- levels(pop(nancycats))
  levels(pop(fsta)) <- levels(pop(nancycats))
  levels(pop(gntx)) <- levels(pop(nancycats))
  levels(pop(stru)) <- levels(pop(nancycats))
  
  # Ensuring that the locus and population summaries are equivalent
  summary_stats <- summary(nancycats)
  expect_equivalent(summary(gpop), summary_stats)
  expect_equivalent(summary(fsta), summary_stats)
  expect_equivalent(summary(gntx), summary_stats)
  expect_equivalent(summary(stru), summary_stats)
})