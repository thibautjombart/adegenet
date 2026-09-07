test_that("genind prediction uses training centres for missing alleles", {
    set.seed(19)
    groups <- factor(rep(c("A", "B"), each = 30))
    dosage <- sapply(seq_len(12), function(j) {
        stats::rbinom(60, 2, ifelse(groups == "A", 0.2, 0.8))
    })
    genotypes <- as.data.frame(apply(dosage, 2, function(z) {
        c("a/a", "a/b", "b/b")[z + 1]
    }))
    g <- df2genind(genotypes, sep = "/", ploidy = 2, pop = groups)

    for (scaled in c(FALSE, TRUE)) {
        fit <- dapc(g, n.pca = 4, n.da = 1, scale = scaled)
        for (n_missing in c(11L, 0L, 1L, 6L, 12L)) {
            batch_a <- g[c(1, 1:20), , drop = FALSE]
            batch_b <- g[c(1, 31:50), , drop = FALSE]
            mask <- locFac(g) %in% locNames(g)[seq_len(n_missing)]
            batch_a@tab[1, mask] <- NA
            batch_b@tab[1, mask] <- NA
            alone <- batch_a[1, , drop = FALSE]
            pa <- predict(fit, batch_a)
            pb <- predict(fit, batch_b)
            expect_equal(unname(pa$posterior[1, ]),
                         unname(pb$posterior[1, ]))
            ps <- predict(fit, alone)
            expect_equal(unname(pa$posterior[1, ]),
                         unname(ps$posterior[1, ]))
            expect_equal(unname(pa$ind.scores[1, ]),
                         unname(pb$ind.scores[1, ]))

            # Independent reference: replace missing frequencies explicitly
            # with the centres recorded during training, before prediction.
            reference <- tab(batch_a, freq = TRUE, NA.method = "asis")
            for (j in seq_len(ncol(reference))) {
                reference[is.na(reference[, j]), j] <- fit$pca.cent[j]
            }
            expected <- predict(fit, reference)
            expect_equal(pa$posterior, expected$posterior)
            expect_equal(pa$ind.scores, expected$ind.scores)
        }
    }
})
