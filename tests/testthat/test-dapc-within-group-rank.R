test_that("DAPC respects the pooled within-group rank bound", {
    set.seed(302)
    x <- matrix(rnorm(30 * 200), 30, 200,
                dimnames = list(NULL, paste0("v", seq_len(200))))
    grp <- factor(rep(c("A", "B", "C"), each = 10))
    expect_warning(fit <- dapc(x, grp, n.pca = 29, n.da = 2),
                   "N - K")
    expect_equal(fit$n.pca, 27)
    expect_equal(fit$n.da, 2)
    expect_true(all(is.finite(fit$posterior)))
    expect_lt(max(fit$eig), 1e10)
    expect_lt(max(abs(fit$loadings)), 1e6)
})

test_that("DAPC reduces data-specific within-group rank deficiency", {
    set.seed(303)
    grp <- factor(rep(c("A", "B", "C"), each = 10))
    signal <- rep(c(0, 1, 2), each = 10)
    x <- cbind(nuisance = rnorm(30), signal = signal)
    expect_warning(fit <- dapc(x, grp, n.pca = 2, n.da = 2),
                   "rank deficient")
    expect_equal(fit$n.pca, 1)
    expect_equal(fit$n.da, 1)
    expect_true(all(is.finite(fit$posterior)))

    expect_error(dapc(matrix(signal, ncol = 1,
                             dimnames = list(NULL, "signal")),
                      grp, n.pca = 1, n.da = 1),
                 "constant within groups")
})

test_that("unused group levels do not alter DAPC dimensions", {
    set.seed(304)
    x <- matrix(rnorm(90), 30, 3,
                dimnames = list(NULL, paste0("v", 1:3)))
    grp <- factor(rep(c("A", "B", "C"), each = 10),
                  levels = c("A", "B", "C", "ghost"))
    fit <- dapc(x, grp, n.pca = 3, n.da = 4)
    expect_equal(levels(fit$grp), c("A", "B", "C"))
    expect_equal(fit$n.da, 2)
    expect_equal(ncol(fit$posterior), 3)
})

test_that("genlight DAPC applies the same rank safeguards", {
    set.seed(305)
    dosage <- matrix(sample(0:2, 30 * 40, replace = TRUE), 30, 40)
    gl <- new("genlight", dosage)
    pop(gl) <- factor(rep(c("A", "B", "C"), each = 10),
                      levels = c("A", "B", "C", "ghost"))
    expect_warning(fit <- dapc(gl, n.pca = 29, n.da = 4,
                               var.contrib = FALSE), "N - K")
    expect_equal(fit$n.pca, 27)
    expect_equal(fit$n.da, 2)
    expect_equal(levels(fit$grp), c("A", "B", "C"))
    expect_true(all(is.finite(fit$posterior)))
    expect_false(is.null(fit$pca.loadings))
    prediction <- predict(fit, gl[1:3, ])
    expect_equal(nrow(prediction$posterior), 3)

    lean <- dapc(gl, n.pca = 5, n.da = 2, pca.info = FALSE,
                 var.contrib = FALSE, var.loadings = FALSE)
    expect_null(lean$pca.loadings)
    expect_true(all(is.finite(lean$posterior)))
})
