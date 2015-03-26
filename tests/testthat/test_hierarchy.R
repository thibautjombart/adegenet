context("strata methods")

test_that("strata methods work for genind objects.", {
  skip_on_cran()

  data(microbov, package = "adegenet")
  expect_null(strata(microbov))
  strata(microbov) <- data.frame(other(microbov))
  breeds <- c("Borgou", "Zebu", "Lagunaire", "NDama", "Somba", "Aubrac", "Bazadais", 
              "BlondeAquitaine", "BretPieNoire", "Charolais", "Gascon", "Limousin", 
              "MaineAnjou", "Montbeliard", "Salers")

  expect_that(length(strata(microbov)), equals(3))
  expect_that(microbov@pop.names, equals(breeds))
  expect_that({microbovsplit <- splitstrata(microbov, ~Pop/Subpop)}, throws_error())
  namestrata(microbov) <- ~Country/Breed/Species
  expect_that(names(strata(microbov)), equals(c("Country", "Breed", "Species")))
  setpop(microbov) <- ~Country/Species
  expect_that(microbov@pop.names, equals(c("AF_BI", "AF_BT", "FR_BT")))
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
  res <- do.call("rbind.genlight", glTest)
  expect_that(res, is_a("genlight"))
  expect_that(nInd(res), equals(100))
  expect_that(nLoc(res), equals(10))
  expect_that(length(strata(res)), equals(3))
  namestrata(res) <- ~Hickory/Dickory/Doc
  expect_that(names(strata(res)), equals(c("Hickory", "Dickory", "Doc")))
})