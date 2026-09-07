test_that("optim.a.score returns values from a non-consecutive PC grid", {
    set.seed(901)
    grp <- factor(rep(LETTERS[1:3], each = 40))
    predictors <- matrix(rnorm(120 * 8), 120, 8)
    predictors[, 1:3] <- predictors[, 1:3] + as.integer(grp) * 0.35
    fit <- dapc(predictors, grp, n.pca = 8, n.da = 2)

    set.seed(1)
    result <- optim.a.score(fit, n.pca = c(2, 5, 6), smart = FALSE,
                            n.sim = 3, n.da = 2, plot = FALSE)

    expect_equal(result$best, 6)
    expect_true(result$best %in% c(2, 5, 6))
    expect_named(result$mean, c("2", "5", "6"))
})

test_that("optim.a.score applies the requested discriminant dimension", {
    set.seed(902)
    grp <- factor(rep(LETTERS[1:3], each = 35))
    predictors <- matrix(rnorm(105 * 6), 105, 6)
    predictors[, 1:2] <- predictors[, 1:2] + as.integer(grp) * 0.4
    fit <- dapc(predictors, grp, n.pca = 6, n.da = 2)

    set.seed(903)
    one.axis <- optim.a.score(fit, n.pca = c(2, 4, 6), smart = FALSE,
                              n.sim = 4, n.da = 1, plot = FALSE)
    set.seed(903)
    two.axes <- optim.a.score(fit, n.pca = c(2, 4, 6), smart = FALSE,
                              n.sim = 4, n.da = 2, plot = FALSE)

    expect_false(isTRUE(all.equal(one.axis$mean, two.axes$mean)))

    set.seed(903)
    omitted <- optim.a.score(fit, n.pca = c(2, 4, 6), smart = FALSE,
                             n.sim = 4, plot = FALSE)
    expect_equal(omitted$mean, two.axes$mean)
})

test_that("optim.a.score validates dimensions and controls", {
    fit <- dapc(iris[, 1:4], iris$Species, n.pca = 4, n.da = 2)

    expect_error(optim.a.score(fit, n.pca = c(2, 5), plot = FALSE),
                 "cannot exceed")
    expect_error(optim.a.score(fit, n.pca = c(1, 2), n.da = 3,
                               plot = FALSE),
                 "n.da must be")
    expect_error(optim.a.score(fit, n.pca = c(1, 2), n.sim = 0,
                               plot = FALSE),
                 "n.sim must be")
    expect_error(optim.a.score(fit, n.pca = c(1, 2), smart = NA,
                               plot = FALSE),
                 "smart must be")
})
