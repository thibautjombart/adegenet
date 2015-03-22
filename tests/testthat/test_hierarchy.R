context("Hierarchy methods")

test_that("Hierarchy methods work for genind objects.", {
  skip_on_cran()

  data(microbov, package = "adegenet")
  expect_null(gethierarchy(microbov))
  sethierarchy(microbov) <- data.frame(other(microbov))
  breeds <- c("Borgou", "Zebu", "Lagunaire", "NDama", "Somba", "Aubrac", "Bazadais", 
              "BlondeAquitaine", "BretPieNoire", "Charolais", "Gascon", "Limousin", 
              "MaineAnjou", "Montbeliard", "Salers")

  expect_that(length(gethierarchy(microbov)), equals(3))
  expect_that(microbov@pop.names, equals(breeds))
  expect_that({microbovsplit <- splithierarchy(microbov, ~Pop/Subpop)}, throws_error())
  namehierarchy(microbov) <- ~Country/Breed/Species
  expect_that(names(gethierarchy(microbov)), equals(c("Country", "Breed", "Species")))
  setpop(microbov) <- ~Country/Species
  expect_that(microbov@pop.names, equals(c("AF_BI", "AF_BT", "FR_BT")))
})