context("Test genclust.em")


test_that("genclust.em gives decent results for F1 Zebu-Salers", {
    skip_on_cran()


    set.seed(1)

    ## Simulate hybrids F1 Zebu/Salers
    data(microbov)
    zebu <- microbov[pop="Zebu"]
    salers <- microbov[pop="Salers"]
    hyb <- hybridize(zebu, salers, n=30)
    x <- repool(zebu, salers, hyb)

    ## run analysis
    res.hyb <- genclust.em(x, k=2, hybrids=TRUE)

    ## check results
    expect_true(res.hyb$converged)
    expect_equal(unname(apply(table(pop(x), res.hyb$group),1,max)),
                 c(50,50,30)) # indiv from original pop all in one group

})




test_that("genclust.em gives decent results for F1 & back-cross Zebu-Salers", {
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
    res.back <- genclust.em(y, k=2, hybrids = TRUE, hybrid.coef = c(.25,.5))
    tab <- table(pop(y), res.back$group)

    ## check results
    expect_true(res.back$converged)
    expect_true(tab[1,1] > 47)
    expect_true(tab[2,2] > 47)
    expect_true(tab[3,4] > 25)
    expect_true(tab[4,3] > 10)
    expect_true(tab[5,5] > 10)

})

