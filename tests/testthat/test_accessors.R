context("Accessor tests")

data("microbov")
strata(microbov) <- data.frame(other(microbov))

test_that("individual accessors work as expected", {
  skip_on_cran()
  expect_that(nInd(microbov), equals(704))
  indNames(microbov)[1] <- "replacement"
  expect_that(indNames(microbov)[1], equals("replacement"))
})

test_that("population accessors work for genind objects", {
  skip_on_cran()
  expect_that(nPop(microbov), equals(15))
  expect_that(popNames(microbov), is_equivalent_to(levels(pop(microbov))))
  popNames(microbov)[1] <- "replacement"
  expect_that(popNames(microbov)[1], equals("replacement"))
  expect_that(unique(head(pop(microbov))), is_equivalent_to(factor("replacement")))
})

test_that("population accessors work for genpop objects", {
  skip_on_cran()
  micpop <- genind2genpop(microbov, quiet = TRUE)
  expect_that(nPop(micpop), equals(15))
  expect_that(popNames(micpop), is_equivalent_to(rownames(micpop@tab)))
  popNames(micpop)[1] <- "replacement"
  expect_that(popNames(micpop)[1], equals("replacement"))
  expect_that(rownames(micpop@tab)[1], equals("replacement"))
})

test_that("locus accessors work for genind objects", {
  skip_on_cran()
  expect_that(nLoc(microbov), equals(30))
  locNames(microbov)[1] <- "replacement"
  expect_that(locNames(microbov)[1], equals("replacement"))
})

test_that("locus accessors work for genpop objects", {
  skip_on_cran()
  micpop <- genind2genpop(microbov, quiet = TRUE)
  expect_that(nLoc(micpop), equals(30))
  locNames(micpop)[1] <- "replacement"
  expect_that(locNames(micpop)[1], equals("replacement"))
})

test_that("'[' method works for genind objects", {
  skip_on_cran()
  two_random_loci    <- sample(locNames(microbov), 2)
  ten_random_samples <- sample(nInd(microbov), 10)
#   cat("\nLoci:", dput(two_random_loci), "\nSamples:",
#       dput(ten_random_samples), "\n")
  pops <- levels(pop(microbov))
  pops <- pops[pops %in% pop(microbov)[ten_random_samples]]
  loci <- microbov@loc.fac[microbov@loc.fac %in% two_random_loci]
  loci <- factor(loci)

  mic10     <- microbov[ten_random_samples]
  mic2Loc   <- microbov[loc = two_random_loci]
  mic2Loc10 <- microbov[ten_random_samples, loc = two_random_loci]

  names(two_random_loci) <- two_random_loci
  two_random_loci <- two_random_loci[levels(loci)]

  expect_that(nInd(mic10), equals(10))
  expect_that(nInd(mic2Loc), equals(nInd(microbov)))
  expect_that(nInd(mic2Loc10), equals(10))

  expect_that(nLoc(mic10), equals(nLoc(microbov)))
  expect_that(nLoc(mic2Loc), equals(2))
  expect_that(nLoc(mic2Loc10), equals(2))

  expect_equivalent(popNames(mic10), pops)
  expect_equivalent(popNames(mic2Loc10), pops)
  expect_that(nPop(mic10), equals(length(pops)))
  expect_that(nPop(mic2Loc10), equals(length(pops)))

  expect_that(length(mic10@loc.fac), equals(ncol(tab(microbov))))
  expect_that(mic2Loc@loc.fac, equals(loci))
  expect_that(mic2Loc10@loc.fac, equals(loci))

  expect_that(mic10@loc.n.all, equals(microbov@loc.n.all))
  expect_that(mic2Loc@loc.n.all, equals(microbov@loc.n.all[two_random_loci]))
  expect_that(mic2Loc10@loc.n.all, equals(microbov@loc.n.all[two_random_loci]))
})

test_that("'[' method works for genind objects with drop = TRUE", {
  skip_on_cran()
  two_random_loci    <- sample(locNames(microbov), 2)
  ten_random_samples <- sample(nInd(microbov), 10)
#   cat("\nLoci:", dput(two_random_loci),
#       "\nSamples:", dput(ten_random_samples), "\n")
  pops <- levels(pop(microbov))
  pops <- pops[pops %in% pop(microbov)[ten_random_samples]]
  j    <- microbov@loc.fac %in% two_random_loci
  loci <- microbov@loc.fac[j]
  loci <- factor(loci)
  loci <- loci[colSums(tab(microbov)[, j], na.rm = TRUE) > 0]
  ten_ind_loci <- colSums(tab(microbov[ten_random_samples, ]), na.rm = TRUE) > 0


  mic10     <- microbov[ten_random_samples, , drop = TRUE]
  mic2Loc   <- microbov[loc = two_random_loci]
  mic2Loc10 <- microbov[ten_random_samples, loc = two_random_loci, drop = TRUE]

  names(two_random_loci) <- two_random_loci
  two_random_loci <- two_random_loci[levels(loci)]

  expect_that(nInd(mic10), equals(10))
  expect_that(nInd(mic2Loc), equals(nInd(microbov)))
  expect_that(nInd(mic2Loc10), equals(10))

  expect_that(nLoc(mic10), equals(nLoc(microbov)))
  expect_that(nLoc(mic2Loc), equals(2))
  expect_that(nLoc(mic2Loc10), equals(2))

  expect_equivalent(popNames(mic10), pops)
  expect_equivalent(popNames(mic2Loc10), pops)
  expect_that(nPop(mic10), equals(length(pops)))
  expect_that(nPop(mic2Loc10), equals(length(pops)))

  ten_ind_loci2 <- factor(mic10@loc.fac[mic10@loc.fac %in% levels(loci)])
  expect_that(length(mic10@loc.fac), equals(ncol(tab(microbov)[, ten_ind_loci])))
  expect_that(mic2Loc@loc.fac, equals(loci))
  expect_that(mic2Loc10@loc.fac, equals(ten_ind_loci2))

  expect_true(all(mic10@loc.n.all <= microbov@loc.n.all))
  expect_that(mic2Loc@loc.n.all, equals(microbov@loc.n.all[two_random_loci]))
  expect_that(names(mic2Loc10@loc.n.all), equals(names(microbov@loc.n.all[two_random_loci])))
  expect_true(all(mic2Loc10@loc.n.all <= microbov@loc.n.all[two_random_loci]))
})

test_that("tab will retain dimensions", {
  skip_on_cran()
  micpop <- genind2genpop(microbov[pop(microbov) %in% popNames(microbov)[1]], quiet = TRUE)
  tabdim <- dim(micpop@tab)
  expect_that(tabdim, equals(dim(tab(micpop))))
  expect_that(tabdim, equals(dim(tab(micpop, freq = TRUE))))
})

test_that("tab will return frequencies for PA data", {
  skip_on_cran()
  x <- read.table(system.file("files/AFLP.txt", package = "adegenet"))
  aflp <- df2genind(x, type = "PA", ploidy = 1, pop = c(rep(1, 4), rep(2, 3)))
  apop <- genind2genpop(aflp)
  atab <- tab(apop, freq = TRUE)
  res  <- tab(apop)/rowSums(tab(apop))
  expect_equivalent(atab, res)
})
