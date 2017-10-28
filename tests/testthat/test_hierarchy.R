context("Strata methods")

test_that("strata methods work for genind objects.", {
  skip_on_cran()

  data(microbov, package = "adegenet")
  expect_null(strata(microbov))
  strata(microbov) <- data.frame(other(microbov))
  breeds <- c("Borgou", "Zebu", "Lagunaire", "NDama", "Somba", "Aubrac", "Bazadais",
              "BlondeAquitaine", "BretPieNoire", "Charolais", "Gascon", "Limousin",
              "MaineAnjou", "Montbeliard", "Salers")

  expect_equal(length(strata(microbov)), 3)
  expect_equal(popNames(microbov), breeds)
  expect_warning(
    expect_error({
      microbovsplit <- splitStrata(microbov, ~Pop/Subpop)
      })
    )
  nameStrata(microbov) <- ~Country/Breed/Species
  expect_equal(names(strata(microbov)), c("Country", "Breed", "Species"))
  setPop(microbov) <- ~Country/Species
  expect_equal(popNames(microbov), c("AF_BI", "AF_BT", "FR_BT"))
})

test_that("strata produce proper errors", {
  skip_on_cran()
  expect_warning(setPop(microbov, ~bippity/boppity/boo))
  strata(microbov) <- data.frame(other(microbov))
  expect_error({strata(microbov) <- data.frame(a = 1)})
  expect_error({addStrata(microbov) <- data.frame(a = 1:10)})
  expect_error(setPop(microbov, ~bippity/boppity/boo))
  expect_error({strata(microbov) <- "a stratum"})
  expect_error({setPop(microbov) <- "thepop"})
})

test_that("strata methods work for genlight objects", {
  skip_on_cran()

  michier <- data.frame(other(microbov))
  make_gl <- function(n = 10, hier = michier){
    objs <- lapply(seq(n), function(x) sample(c(0, 1, NA), 10, replace = TRUE, prob = c(0.49, 0.49, 0.01)))
    return(new("genlight", objs, strata = hier[sample(704, 10), sample(3, 2)], parallel = FALSE))
  }
  set.seed(9999)
  glTest <- lapply(1:10, function(x, y, z) make_gl(y, z), 10, michier)
  res <- do.call("rbind.genlight", c(glTest, parallel = FALSE))
  expect_is(res, "genlight")
  expect_equal(nInd(res), 100)
  expect_equal(nLoc(res), 10)
  expect_equal(length(strata(res)), 3)
  nameStrata(res) <- ~Hickory/Dickory/Doc
  expect_equal(names(strata(res)), c("Hickory", "Dickory", "Doc"))
})
