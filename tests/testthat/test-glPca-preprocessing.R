test_that("glPca handles invariant loci consistently when scaling", {
    genotypes <- matrix(
        c(0, 0, 0, 0,
          0, 1, 2, 1,
          0, NA, 2, 1),
        nrow = 4
    )
    colnames(genotypes) <- c("fixed", "variable", "incomplete")
    rownames(genotypes) <- paste0("individual", seq_len(nrow(genotypes)))
    x <- new("genlight", genotypes)
    indNames(x) <- rownames(genotypes)
    locNames(x) <- colnames(genotypes)

    expect_warning(
        compiled <- glPca(x, center = TRUE, scale = TRUE, nf = 2,
                          loadings = TRUE, useC = TRUE),
        "Null variances"
    )
    expect_warning(
        interpreted <- glPca(x, center = TRUE, scale = TRUE, nf = 2,
                             loadings = TRUE, useC = FALSE),
        "Null variances"
    )

    expect_true(all(is.finite(compiled$scores)))
    expect_true(all(is.finite(compiled$loadings)))
    expect_equal(compiled$loadings[1, ], c(Axis1 = 0, Axis2 = 0))
    expect_equal(abs(compiled$scores), abs(interpreted$scores), tolerance = 1e-10)
    expect_equal(abs(compiled$loadings), abs(interpreted$loadings),
                 tolerance = 1e-10)
})

test_that("unscaled invariant loci agree between glPca implementations", {
    genotypes <- matrix(
        c(2, 2, 2, 2,
          0, 1, 2, 1),
        nrow = 4
    )
    x <- new("genlight", genotypes)

    expect_warning(
        compiled <- glPca(x, center = FALSE, scale = TRUE, nf = 2,
                          loadings = FALSE, useC = TRUE),
        "Null variances"
    )
    expect_warning(
        interpreted <- glPca(x, center = FALSE, scale = TRUE, nf = 2,
                             loadings = FALSE, useC = FALSE),
        "Null variances"
    )

    expect_true(all(is.finite(compiled$scores)))
    expect_equal(abs(compiled$scores), abs(interpreted$scores), tolerance = 1e-10)
})
