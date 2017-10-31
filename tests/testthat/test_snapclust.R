context("Test snapclust")


test_that("snapclust gives decent results for F1 Zebu-Salers", {
    skip_on_cran()


    set.seed(1)

    ## Simulate hybrids F1 Zebu/Salers
    data(microbov)
    zebu <- microbov[pop="Zebu"]
    salers <- microbov[pop="Salers"]
    hyb <- hybridize(zebu, salers, n=30)
    x <- repool(zebu, salers, hyb)

    ## run analysis
    res.hyb <- snapclust(x, k=2, hybrids=TRUE)

    ## check results
    expect_true(res.hyb$converged)
    expect_equal(unname(apply(table(pop(x), res.hyb$group),1,max)),
                 c(50,50,30)) # indiv from original pop all in one group

})




test_that("snapclust gives decent results for F1 & back-cross Zebu-Salers", {
    skip_on_cran()

    set.seed(1)

    ## Simulate hybrids F1 Zebu/Salers
    data(microbov)
    zebu <- microbov[pop="Zebu"]
    salers <- microbov[pop="Salers"]
    hyb <- hybridize(zebu, salers, n=30)
    x <- repool(zebu, salers, hyb)


    ## Simulate hybrids backcross (F1 / parental)
    f1.zebu <- hybridize(hyb, zebu, 20, pop = "f1.zebu")
    f1.salers <- hybridize(hyb, salers, 25, pop = "f1.salers")
    y <- repool(x, f1.zebu, f1.salers)


    ## method with back-cross
    res.back <- snapclust(y, k=2, hybrids = TRUE, hybrid.coef = c(.25,.5))
    tab <- table(pop(y), res.back$group)

    ## check results
    expect_true(res.back$converged)
    expect_true(tab[1,1] > 47)
    expect_true(tab[2,2] > 47)
    expect_true(tab[3,4] > 25)
    expect_true(tab[4,3] > 10)
    expect_true(tab[5,5] > 10)

})


test_that("snapclust.choose.k will recognize objects that inherit genind objects", {
  skip_on_cran()
  skip_if_not_installed("poppr")
  requireNamespace("poppr", quietly = TRUE)
  data(microbov)
  zebu <- poppr::as.genclone(microbov[pop = "Zebu"])
  expect_is(zebu, "genclone")
  res <- snapclust.choose.k(2, zebu)
  expect_length(res, 2)
  expect_is(res, "numeric")
  unloadNamespace("poppr")
})

test_that("snapclust.choose.k will ignore any extra genind objects supplied", {
  skip_on_cran()
  data(microbov)
  expect_warning({
    res <- snapclust.choose.k(2, microbov[pop = "Zebu"], microbov[pop = "Salers"])  
  }, "Too many genind objects provided")
  expect_is(res, "numeric")
  expect_length(res, 2)
})

test_that("snapclust.choose.k will throw an error if there are no genind objects", {
  skip_on_cran()
  expect_error(snapclust.choose.k(2, 1), "No genind provided")
})