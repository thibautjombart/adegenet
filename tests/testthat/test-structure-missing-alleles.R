test_that("STRUCTURE missing codes match complete allele tokens in both layouts", {
    for (missing in c("0", "-9", ".")) {
        for (one.row in c(TRUE, FALSE)) {
            # L1 has complete, partially missing, and fully missing genotypes.
            # L2 keeps all individuals in the output for checking missing calls.
            a <- matrix(c("118", "120", "126", "132", missing, "128",
                          missing, missing), ncol=2, byrow=TRUE)
            lines <- character()
            for (i in 1:4) {
                if (one.row) {
                    lines <- c(lines, paste(c(paste0("ind",i), 1, a[i,], "140", "142"), collapse=" "))
                } else {
                    lines <- c(lines, paste(c(paste0("ind",i), 1, a[i,1], "140"), collapse=" "),
                               paste(c(paste0("ind",i), 1, a[i,2], "142"), collapse=" "))
                }
            }
            f <- tempfile(fileext=".stru")
            writeLines(lines, f)
            g <- read.structure(f, n.ind=4, n.loc=2, onerowperind=one.row,
                                col.lab=1, col.pop=2, row.marknames=0,
                                NA.char=missing, ask=FALSE, quiet=TRUE)
            unlink(f)
            expect_identical(indNames(g), paste0("ind",1:4))
            expect_true(all(c("118","120","126","132") %in% alleles(g)[[1]]))
            counts <- seploc(g, res.type="matrix")[[1]]
            expect_equal(unname(rowSums(counts)[1:2]), c(2,2))
            expect_true(all(is.na(counts[3:4,,drop=FALSE])))
            expect_equal(unname(rowSums(seploc(g,res.type="matrix")[[2]])), rep(2,4))
        }
    }
})

test_that("two-row STRUCTURE input preserves single-locus dimensions", {
    f <- tempfile(fileext=".stru")
    on.exit(unlink(f))
    writeLines(c("ind1 1 118", "ind1 1 120", "ind2 1 126", "ind2 1 128"), f)
    g <- read.structure(f, n.ind=2, n.loc=1, onerowperind=FALSE,
                        col.lab=1, col.pop=2, row.marknames=0,
                        NA.char="0", ask=FALSE, quiet=TRUE)
    expect_identical(indNames(g), c("ind1","ind2"))
    expect_equal(nLoc(g), 1L)
    expect_equal(unname(rowSums(tab(g))), c(2,2))
})
