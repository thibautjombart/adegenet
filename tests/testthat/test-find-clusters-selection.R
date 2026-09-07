test_that("interactive percentage selection is actually applied", {
    set.seed(401)
    x <- as.data.frame(matrix(rnorm(40 * 5), 40, 5))
    input <- textConnection("80")
    old <- getOption("adegenet.testcon")
    options(adegenet.testcon = input)
    on.exit({
        options(adegenet.testcon = old)
        close(input)
    }, add = TRUE)
    set.seed(402)
    interactive <- find.clusters(x, pca.select = "percVar",
                                 n.clust = 2, n.start = 5)
    set.seed(402)
    explicit <- find.clusters(x, pca.select = "percVar", perc.pca = 80,
                              n.clust = 2, n.start = 5)
    expect_equal(interactive$grp, explicit$grp)
    expect_equal(interactive$size, explicit$size)
})

test_that("cluster and PCA dimensions are validated against the data", {
    set.seed(403)
    x <- as.data.frame(matrix(rnorm(8 * 3), 8, 3))
    expect_warning(ans <- find.clusters(x, n.pca = 2,
                                        max.n.clust = 100,
                                        choose.n.clust = FALSE,
                                        criterion = "min", n.start = 2),
                   "Reducing max.n.clust")
    expect_lte(length(ans$Kstat), 8)
    expect_error(find.clusters(x, n.pca = 2, n.clust = 8),
                 "n.clust must be an integer")
    expect_error(find.clusters(x, n.pca = 0, n.clust = 2),
                 "n.pca must be a positive integer")
    expect_error(find.clusters(x, n.pca = 2, n.clust = 2.5),
                 "n.clust must be an integer")
    x[1, 1] <- NA
    expect_error(find.clusters(x, n.pca = 2, n.clust = 2),
                 "complete, finite")
})

test_that("automatic cluster criteria return actual valid K values", {
    expect_equal(.choose_n_clusters(c(10, 7, 5), "min"), 3)
    expect_warning(
        expect_equal(.choose_n_clusters(c(10, 7, 5), "goesup"), 3),
        "no stopping point"
    )
    expect_warning(
        expect_equal(.choose_n_clusters(c(10, 7),
                                                    "smoothNgoesup"), 2),
        "no stopping point"
    )
    # The first position represents K = 1; it must never become K = 0.
    expect_equal(.choose_n_clusters(c(1, 10, 20), "goodfit"), 1)
    selected <- .choose_n_clusters(c(15, 8, 7, 6.8),
                                              "diffNgroup")
    expect_true(selected %in% 1:4)
})

test_that("fixed seeds reproduce k-means cluster memberships", {
    set.seed(404)
    x <- as.data.frame(matrix(rnorm(50 * 4), 50, 4))
    set.seed(405)
    first <- find.clusters(x, n.pca = 3, n.clust = 3, n.start = 10)
    set.seed(405)
    second <- find.clusters(x, n.pca = 3, n.clust = 3, n.start = 10)
    expect_equal(first$grp, second$grp)
    expect_equal(first$size, second$size)
})
