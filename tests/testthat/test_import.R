context("Import Tests")

test_that("df2genind works with haploids", {
  skip_on_cran()
  x <- matrix(as.character(sample(20)), nrow = 10, ncol = 2)
  res <- df2genind(x, ploidy = 1)
  expect_equal(sum(res@loc.n.all), 20)
  expect_equal(nInd(res), 10)
  expect_equal(nLoc(res), 2)
  resdf <- genind2df(res)
  expect_equal(unlist(resdf, use.names=FALSE), as.vector(x))
})

test_that("df2genind makes sense for given example", {
  skip_on_cran()
  df <- data.frame(locusA=c("11","11","12","32"),
                   locusB=c(NA,"34","55","15"),
                   locusC=c("22","22","21","22"))
  row.names(df) <- .genlab("genotype",4)
  obj <- df2genind(df, ploidy=2, ncode=1)
  expect_equal(nInd(obj), 4)
  expect_equal(nLoc(obj), 3)
  expect_equal(locNames(obj), colnames(df))
  expect_equal(indNames(obj), rownames(df))
  expect_equivalent(obj@loc.n.all, c(3, 4, 2))
  objdf <- genind2df(obj)
  expect_equivalent(df, df)
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
  expect_equivalent(g["730", grepl("REN169D01", colnames(g))],
                    c(REN169D01.204 = NA_integer_,
                      REN169D01.208 = NA_integer_,
                      REN169D01.198 = NA_integer_))
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
  expect_equivalent(unique(rowSums(tab(xgid))), 4)
  expect_equivalent(genind2df(xgid, sep = "/"), xdf)
})


test_that("read.X functions work as expected", {
  skip_on_cran()
  suppressWarnings({
    gpop <- read.genepop(system.file("files/nancycats.gen",package="adegenet"), quiet = TRUE)
  })
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

test_that("read.genpop can import duplicate names", {
  skip_on_cran()
  x <- "
  Microsat on Chiracus radioactivus, a pest species
     Loc1, Loc2, Loc3, Y-linked, Loc4
POP
AA8, 0405 0711 0304 0000      0505
AA9, 0405 0609 0208 0000      0505
A10, 0205 0609 0101 0000      0305
A11, 0405 0606 0102 0000      0504
A12, 0202 0609 0105 0000      0507
A13, 0505 0909 0107 0000      0505
A14, 0202 0609 0207 0000      0503
A15, 0405 0609 0101 0000      0505
Pop
AF, 0000 0000 0000 0000      0505
AF, 0205 0307 0102 0000      0505
AF, 0202 0609 0202 0000      0505
AF, 0205 0909 0000 0000      0505
AF, 0205 0307 0202 0000      0505
AF, 0505 0303 0102 0000      0505
AF, 0205 0700 0000 0000      0505
AF, 0505 0900 0000 0000      0405
AF, 0205 0600 0000 0000      0505
AF, 0505 0606 0202 0000      0505
pop
C45, 0505 0606 0202 0000      0505
C45, 0505 0909 0202 0000      0505
C45, 0505 0306 0202 0000      0505
C45, 0505 0909 0102 0000      0405
C45, 0205 0303 0202 0000      0505
C45, 0205 0909 0202 0000      0405
  "
  tmp <- tempfile(fileext = ".gen")
  cat(x, file = tmp)
  expect_warning(gp <- read.genepop(tmp, quiet = TRUE))
  expect_identical(indNames(gp), .genlab("", nInd(gp)))

})



test_that("df2genind can handle periods in input", {
  skip_on_cran()
  dat <-
    data.frame(
      pop = c(1, 1, 2, 2),
      loc1 = c("1/1", "1/2", "1.1/2", "2/2"),
      loc2 = c("1/1", "1/2", "1/2", "2/2")
    )
  expect_warning(datgi <- df2genind(dat[, -1], sep = "/", pop = dat[, 1]))
  expect_equivalent(alleles(datgi)$loc1, c("1", "2", "1_1"))
})

test_that("df2genind can handle periods in input with underscore separator", {
  skip_on_cran()
  dat <-
    data.frame(
      pop = c(1, 1, 2, 2),
      loc1 = c("1/1", "1/2", "1.1/2", "2/2"),
      loc2 = c("1/1", "1/2", "1/2", "2/2")
    )
  dat <- apply(dat, 2, function(i) gsub("/", "_", i))
  expect_warning(datgi <- df2genind(dat[, -1], sep = "_", pop = dat[, 1]))
  expect_equivalent(alleles(datgi)$loc1, c("1", "2", "1p1"))
})


test_that("different imports sort populations in the same way", {
    skip_on_cran()

    ## read nancycats data from different formats
  expect_output({
    x.str <- read.structure(
      system.file("files/nancycats.str", package = "adegenet"),
      onerowperind = FALSE,
      n.ind = 237,
      n.loc = 9,
      col.lab = 1,
      col.pop = 2,
      ask = FALSE
    )
  }, "Converting data")
  expect_output({
    expect_warning({
      x.gen <- read.genepop(system.file("files/nancycats.gen", package = "adegenet"))
    }, "Duplicate individual names detected.")
  }, "Converting data")
  
  expect_output({
    x.dat <- read.fstat(system.file("files/nancycats.dat", package = "adegenet"))
  }, "Converting data")
  
  expect_output({
    x.gtx <- read.genetix(system.file("files/nancycats.gtx", package = "adegenet"))
  }, "Converting data")
  
  
    ## check that the pop are identical:

    ## we use 'table(pop(...))' because individuals may be sorted differently in the files, so
    ## 'pop(...)' may be different

    identical(table(pop(x.gen)), table(pop(x.str)))
    identical(table(pop(x.gen)), table(pop(x.dat)))
    identical(table(pop(x.gen)), table(pop(x.gtx)))
})


test_that("ensure importing structure files with numbers for locus names imports correectly", {
  skip_on_cran()
  
  # Column names should have no extra characters in front of them. Your IDE
  # may be adding them, so watch out!
  tmp <- tempfile(fileext = ".stru")
  cat(
"		1_25	8_54	1358_15	1363_12	1368_57	1369_41	1372_14	1373_9	1377_42	1378_53	1379_10	1382_37	1386_27	1398_46	1400_9	1401_25	1403_13	1404_17	1409_42	1416_48	1419_11	1421_14	1423_5	1424_74	1426_55	1429_46	1432_23	1435_30	1436_7	1438_9	1443_37
A_KH1584	A	1	4	4	1	1	3	2	4	4	2	3	3	2	4	1	3	1	1	2	3	1	4	4	3	2	2	3	4	4	4	2
A_KH1584	A	1	4	4	1	1	3	2	4	4	4	3	3	4	4	1	3	1	3	2	3	3	4	4	3	4	2	3	4	4	4	2
C_KH1059	C	0	4	4	1	1	3	2	4	4	2	1	3	2	4	1	3	1	3	2	3	3	2	4	3	2	2	3	2	4	4	2
C_KH1059	C	0	4	4	1	1	3	2	4	4	4	3	3	2	4	1	3	1	3	2	3	3	4	4	3	2	2	3	4	4	4	2
M_KH1834	M	0	2	2	1	1	3	2	4	4	2	3	3	2	4	1	3	1	1	2	3	3	4	4	3	2	2	3	2	4	4	2
M_KH1834	M	0	4	4	1	3	3	2	4	4	2	3	3	2	4	1	3	1	3	2	3	3	4	4	3	2	4	3	4	4	4	2
M_KH1837	M	1	4	4	1	1	3	2	4	4	0	3	3	2	2	1	3	1	1	2	3	3	4	4	3	4	2	3	4	4	4	2
M_KH1837	M	1	4	4	1	3	3	2	4	4	0	3	3	4	4	1	3	1	3	2	3	3	4	4	3	4	2	3	4	4	4	2", 
      file = tmp)
  
  xy1 <- read.structure(
    tmp,
    NA.char = "0",
    n.ind = 4,
    n.loc = 31,
    onerowperind = FALSE,
    col.lab = 1,
    col.pop = 2,
    row.marknames = 1,
    sep = "\t",
    col.others = 0,
    quiet = TRUE
  )
  
  x1 <- tab(xy1)
  # should return all 1, incorrect is NA
  expect_true(all(x1[, grepl("1401_25", colnames(x1)), drop = FALSE] == 1))

  tmp2 <- tempfile(fileext = ".stru")
  # Column names should have no extra characters in front of them. Your IDE
  # may be adding them, so watch out!
  cat(
"		X1_25	X8_54	X1358_15	X1363_12	X1368_57	X1369_41	X1372_14	X1373_9	X1377_42	X1378_53	X1379_10	X1382_37	X1386_27	X1398_46	X1400_9	X1401_25	X1403_13	X1404_17	X1409_42	X1416_48	X1419_11	X1421_14	X1423_5	X1424_74	X1426_55	X1429_46	X1432_23	X1435_30	X1436_7	X1438_9	X1443_37
A_KH1584	A	1	4	4	1	1	3	2	4	4	2	3	3	2	4	1	3	1	1	2	3	1	4	4	3	2	2	3	4	4	4	2
A_KH1584	A	1	4	4	1	1	3	2	4	4	4	3	3	4	4	1	3	1	3	2	3	3	4	4	3	4	2	3	4	4	4	2
C_KH1059	C	0	4	4	1	1	3	2	4	4	2	1	3	2	4	1	3	1	3	2	3	3	2	4	3	2	2	3	2	4	4	2
C_KH1059	C	0	4	4	1	1	3	2	4	4	4	3	3	2	4	1	3	1	3	2	3	3	4	4	3	2	2	3	4	4	4	2
M_KH1834	M	0	2	2	1	1	3	2	4	4	2	3	3	2	4	1	3	1	1	2	3	3	4	4	3	2	2	3	2	4	4	2
M_KH1834	M	0	4	4	1	3	3	2	4	4	2	3	3	2	4	1	3	1	3	2	3	3	4	4	3	2	4	3	4	4	4	2
M_KH1837	M	1	4	4	1	1	3	2	4	4	0	3	3	2	2	1	3	1	1	2	3	3	4	4	3	4	2	3	4	4	4	2
M_KH1837	M	1	4	4	1	3	3	2	4	4	0	3	3	4	4	1	3	1	3	2	3	3	4	4	3	4	2	3	4	4	4	2", 
      file = tmp2)
  
  xy2 <- read.structure(
    tmp2,
    NA.char = "0",
    n.ind = 4,
    n.loc = 31,
    onerowperind = FALSE,
    col.lab = 1,
    col.pop = 2,
    row.marknames = 1,
    sep = "\t",
    col.others = 0,
    quiet = TRUE
  )
  
  x2 <- tab(xy2)
  
  # should return all 1
  expect_true(all(x2[, grepl("1401_25", colnames(x2)), drop = FALSE] == 1))
  
})


test_that("df2genind throws a warning when the user borks the ncode", {
  dat <- data.frame(stringsAsFactors = FALSE,
                    A = c("A5A5", "A5A4", "A5A5", "A5A5", "A5A4", "A5A4", "A5A5",
                          "A5A4", "A5A4", "A5A5"),
                    B = c("B1B1", "B1B1", "B1B1", "B1B1", "B1B1", "B1B1", "B1B1",
                          "B1B1", "B1B2", "B1B1"),
                    C = c("C2C2", "C2C2", "C2C2", "C2C2", "C2C2", "C2C2", "C2C2",
                          "C2C2", "C2C2", "C2C2"),
                    D = c("D2D5", "D5D3", "D5D5", "D2D5", "D2D1", "D5D3", "D1D1",
                          "D2D2", "D2D5", "D2D4")
  )
  dat
  
  expect_warning(df2genind(dat, ploidy = 2, ncode = 1), "observed allele dosage \\(4-4\\)")
  
})