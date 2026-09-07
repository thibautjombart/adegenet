test_that("validation individuals do not contribute to PCA or DAPC fits", {
    set.seed(104)
    x <- matrix(rnorm(48 * 6), 48, 6)
    grp <- factor(rep(c("A", "B", "C"), each = 16))
    keep <- unlist(lapply(split(seq_len(48), grp), head, 10))
    dat <- list(DATA = x, GRP = grp, KEEP = keep,
                CENTER = TRUE, SCALE = TRUE)
    # Instrument a private copy so this test does not modify the namespace.
    worker <- adegenet:::.boot_dapc_pred
    env <- new.env(parent = environment(worker))
    seen <- list()
    fits <- 0L
    env$dudi.pca <- function(df, ...) {
        seen[[length(seen) + 1L]] <<- df
        ade4::dudi.pca(df, ...)
    }
    env$dapc <- function(...) {
        fits <<- fits + 1L
        adegenet::dapc(...)
    }
    environment(worker) <- env
    answer <- worker(dat, n.pca = 2, n.da = 1)
    expect_length(seen, 1)
    expect_equal(seen[[1]], x[keep, , drop = FALSE])
    expect_equal(fits, 1L)
    reference <- dapc(x[keep, ], grp[keep], n.pca = 2, n.da = 1,
                      center = TRUE, scale = TRUE)
    expected <- mean(predict(reference, x[-keep, ])$assign == grp[-keep])
    expect_equal(answer, expected)
    dat$DATA[-keep, ] <- dat$DATA[-keep, ] + 1000
    worker(dat, n.pca = 2, n.da = 1)
    expect_equal(seen[[2]], seen[[1]])
})

test_that("final DAPC preserves requested preprocessing and n.da", {
    set.seed(105)
    x <- matrix(rnorm(60 * 6), 60, 6)
    x <- sweep(x, 2, seq_len(6), "*") + 10
    grp <- factor(rep(c("A", "B", "C"), each = 20))
    ans <- xvalDapc(x, grp, n.pca = 2, n.pca.max = 2, n.da = 1,
                    training.set = 0.7, center = FALSE, scale = TRUE,
                    n.rep = 1, xval.plot = FALSE)
    reference <- dapc(x, grp, n.pca = 2, n.da = 1,
                      center = FALSE, scale = TRUE)
    expect_equal(ans$DAPC$pca.cent, reference$pca.cent)
    expect_equal(ans$DAPC$pca.norm, reference$pca.norm)
    expect_equal(ans$DAPC$n.da, 1)
    expect_equal(ans$DAPC$posterior, reference$posterior)
    expect_equal(nrow(ans[[1]]), 1L)
    expect_equal(as.numeric(ans[[5]]), abs(ans[[1]]$success - 1))
})

test_that("rank-infeasible candidates are not scored under a false PC count", {
    set.seed(106)
    z <- rnorm(60)
    x <- cbind(locus1 = z, locus2 = z, locus3 = z)
    grp <- factor(rep(c("A", "B"), each = 30))
    expect_warning(ans <- xvalDapc(x, grp, n.pca = c(1, 2),
                                   n.pca.max = 2, n.rep = 2,
                                   xval.plot = FALSE), "excluded")
    expect_true(all(is.na(ans[[1]]$success[ans[[1]]$n.pca == 2])))
    expect_equal(ans$DAPC$n.pca, 1)
    expect_true(is.na(ans[[3]]["2"]))
    expect_true(is.na(ans[[5]]["2"]))
    expect_error(xvalDapc(x, grp, n.pca = 2, n.pca.max = 2,
                          n.rep = 2, xval.plot = FALSE), "No candidate")
})

test_that("PC candidates respect the training within-group LDA rank", {
    set.seed(107)
    x <- matrix(rnorm(30 * 40), 30, 40)
    grp <- factor(rep(c("A", "B", "C"), each = 10))
    # With 21 training individuals in 3 groups, pooled within-group
    # covariance has rank at most 18, even though PCA supports 20 axes.
    ans <- xvalDapc(x, grp, n.pca = c(18, 19, 20), n.pca.max = 20,
                    training.set = 0.7, n.rep = 1, xval.plot = FALSE)
    expect_equal(unique(ans[[1]]$n.pca), 18)
    expect_lte(ans$DAPC$n.pca, 18)
})

test_that("missing data are rejected explicitly and splits retain groups", {
    x <- matrix(seq_len(40), 10, 4)
    x[1, 1] <- NA
    expect_error(xvalDapc(x, rep(1:2, each = 5)), "missing values are not imputed")
    grp <- factor(c("singleton", rep("A", 3), rep("B", 4)))
    for (fraction in c(0.01, 0.99)) {
        keep <- unlist(lapply(levels(grp), adegenet:::.group_sampler,
                               grp = grp, training.set = fraction))
        expect_true(1L %in% keep)
        expect_equal(length(unique(grp[keep])), 3L)
        expect_equal(length(unique(grp[-keep])), 2L)
    }
})
