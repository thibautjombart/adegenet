context("Accessor tests")

data("microbov")
strata(microbov) <- data.frame(other(microbov))

test_that("individual accessors work as expected", {
  skip_on_cran()
  expect_equal(nInd(microbov), 704)
  indNames(microbov)[1] <- "replacement"
  expect_equal(indNames(microbov)[1], "replacement")
})

test_that("population accessors work for genind objects", {
  skip_on_cran()
  expect_equal(nPop(microbov), 15)
  expect_equivalent(popNames(microbov), levels(pop(microbov)))
  popNames(microbov)[1] <- "replacement"
  expect_equal(popNames(microbov)[1], "replacement")
  expect_equivalent(unique(head(pop(microbov))), factor("replacement"))
})

test_that("population accessors work for genpop objects", {
  skip_on_cran()
  micpop <- genind2genpop(microbov, quiet = TRUE)
  expect_equal(nPop(micpop), 15)
  expect_equivalent(popNames(micpop), rownames(micpop@tab))
  popNames(micpop)[1] <- "replacement"
  expect_equal(popNames(micpop)[1], "replacement")
  expect_equal(rownames(micpop@tab)[1], "replacement")
})

test_that("locus accessors work for genind objects", {
  skip_on_cran()
  expect_equal(nLoc(microbov), 30)
  locNames(microbov)[1] <- "replacement"
  expect_equal(locNames(microbov)[1], "replacement")
})

test_that("locus accessors work for genpop objects", {
  skip_on_cran()
  micpop <- genind2genpop(microbov, quiet = TRUE)
  expect_equal(nLoc(micpop), 30)
  locNames(micpop)[1] <- "replacement"
  expect_equal(locNames(micpop)[1], "replacement")
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

  expect_equal(nInd(mic10), 10)
  expect_equal(nInd(mic2Loc), nInd(microbov))
  expect_equal(nInd(mic2Loc10), 10)

  expect_equal(nLoc(mic10), nLoc(microbov))
  expect_equal(nLoc(mic2Loc), 2)
  expect_equal(nLoc(mic2Loc10), 2)

  expect_equivalent(popNames(mic10), pops)
  expect_equivalent(popNames(mic2Loc10), pops)
  expect_equal(nPop(mic10), length(pops))
  expect_equal(nPop(mic2Loc10), length(pops))

  expect_equal(length(mic10@loc.fac), ncol(tab(microbov)))
  expect_equal(mic2Loc@loc.fac, loci)
  expect_equal(mic2Loc10@loc.fac, loci)

  expect_equal(mic10@loc.n.all, microbov@loc.n.all)
  expect_equal(mic2Loc@loc.n.all, microbov@loc.n.all[two_random_loci])
  expect_equal(mic2Loc10@loc.n.all, microbov@loc.n.all[two_random_loci])
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

  expect_equal(nInd(mic10), 10)
  expect_equal(nInd(mic2Loc), nInd(microbov))
  expect_equal(nInd(mic2Loc10), 10)

  expect_equal(nLoc(mic10), nLoc(microbov))
  expect_equal(nLoc(mic2Loc), 2)
  expect_equal(nLoc(mic2Loc10), 2)

  expect_equivalent(popNames(mic10), pops)
  expect_equivalent(popNames(mic2Loc10), pops)
  expect_equal(nPop(mic10), length(pops))
  expect_equal(nPop(mic2Loc10), length(pops))

  ten_ind_loci2 <- factor(mic10@loc.fac[mic10@loc.fac %in% levels(loci)])
  expect_equal(length(mic10@loc.fac), ncol(tab(microbov)[, ten_ind_loci]))
  expect_equal(mic2Loc@loc.fac, loci)
  expect_equal(mic2Loc10@loc.fac, ten_ind_loci2)

  expect_true(all(mic10@loc.n.all <= microbov@loc.n.all))
  expect_equal(mic2Loc@loc.n.all, microbov@loc.n.all[two_random_loci])
  expect_equal(names(mic2Loc10@loc.n.all), names(microbov@loc.n.all[two_random_loci]))
  expect_true(all(mic2Loc10@loc.n.all <= microbov@loc.n.all[two_random_loci]))
})

test_that("tab will retain dimensions", {
  skip_on_cran()
  micpop <- genind2genpop(microbov[pop(microbov) %in% popNames(microbov)[1]], quiet = TRUE)
  tabdim <- dim(micpop@tab)
  expect_equal(tabdim, dim(tab(micpop)))
  expect_equal(tabdim, dim(tab(micpop, freq = TRUE)))
})

test_that("tab will return frequencies for PA data", {
  skip_on_cran()
  x <- read.table(system.file("files/AFLP.txt", package = "adegenet"))
  aflp <- df2genind(x, type = "PA", ploidy = 1, pop = c(rep(1, 4), rep(2, 3)))
  apop <- genind2genpop(aflp, quiet = TRUE)
  atab <- tab(apop, freq = TRUE)
  res  <- tab(apop)/rowSums(tab(apop))
  expect_equivalent(atab, res)
})

test_that("subsettors give one warning for individuals", {
  skip_on_cran()
  expect_warning(microbov[c("bippity", "hop", "bop")], "the following specified individuals do not exist: bippity, hop, bop")
  expect_warning(microbov[c("bippity", "hop", "bop")], "no individual selected - ignoring")
})


test_that("subsettors give one warning for loci", {
  skip_on_cran()
  expect_warning(microbov[loc = c("bippity", "hop", "bop")], "the following specified loci do not exist: bippity, hop, bop")
  expect_warning(microbov[loc = c("bippity", "hop", "bop")], "no loci selected - ignoring")
})
  
test_that("subsettors give one warning for populations", {
  skip_on_cran()
  expect_warning(microbov[pop = c("bippity", "hop", "bop")], "the following specified populations do not exist: bippity, hop, bop")
  expect_warning(microbov[pop = c("bippity", "hop", "bop")], "no populations selected - ignoring")
})

test_that("subsettors give one warning for genpop objects", {
  skip_on_cran()
  micropop <- genind2genpop(microbov, quiet = TRUE)
  expect_warning(micropop[c("bippity", "hop", "bop")], "the following specified populations do not exist: bippity, hop, bop")
  expect_warning(micropop[c("bippity", "hop", "bop")], "no population selected - ignoring")
})