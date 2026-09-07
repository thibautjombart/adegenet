test_that("split-allele export preserves incomplete genotypes and rows", {
    x <- df2genind(data.frame(L1=c("A/A","A/B","B/B","A/B"),
                             L2=c("C/D","C/C","D/D","C/D")), sep="/", ploidy=2)
    cols <- which(locFac(x) == locNames(x)[1])
    x@tab[2,cols] <- c(0L,1L)
    x@tab[3,cols] <- 0L
    x@tab[4,cols] <- NA_integer_
    y <- genind2df(x, oneColPerAll=TRUE, usepop=FALSE)
    expect_identical(rownames(y), indNames(x))
    expect_identical(y[[1]], c("A","B",NA_character_,NA_character_))
    expect_identical(y[[2]], c("A",NA_character_,NA_character_,NA_character_))
    expect_identical(y[[3]], c("C","C","D","C"))
    expect_identical(y[[4]], c("D","C","D","D"))
    x@tab[1,cols] <- c(3L,0L)
    expect_error(genind2df(x, oneColPerAll=TRUE), "exceed ploidy")
})

test_that("rupica exports every individual without recycling allele copies", {
    data(rupica, package="adegenet", envir=environment())
    expect_warning(y <- genind2df(rupica, oneColPerAll=TRUE, usepop=FALSE), NA)
    expect_equal(nrow(y), nInd(rupica))
    expect_equal(ncol(y), 2*nLoc(rupica))
    expect_identical(rownames(y), indNames(rupica))
    for (j in seq_len(nLoc(rupica))) {
        counts <- seploc(rupica,res.type="matrix")[[j]]
        for (i in seq_len(nInd(rupica))) {
            emitted <- unlist(y[i, c(2*j-1,2*j)], use.names=FALSE)
            expected <- if(anyNA(counts[i,])) character() else
                rep(alleles(rupica)[[j]], counts[i,])
            expect_identical(sort(emitted[!is.na(emitted)]), sort(unname(expected)))
        }
    }
})

test_that("presence/absence exports honour population overrides", {
    x <- df2genind(data.frame(L1=c(1,0,1), L2=c(0,1,1)), type="PA",
                  ncode=1, pop=c("old","old","old"))
    replacement <- factor(c("a","b","c"))
    expect_identical(genind2df(x,pop=replacement)$pop, replacement)
    expect_identical(genind2df(x)$pop, pop(x))
    expect_equal(genind2df(x,pop=replacement,usepop=FALSE), tab(x))
})

test_that("split export preserves literal NA and separator allele labels", {
    x <- df2genind(data.frame(L1=c("AA:BB","AA:AA")),sep=":",ploidy=2)
    x@all.names[[1]] <- c("NA","A/B")
    y <- genind2df(x,oneColPerAll=TRUE,usepop=FALSE)
    expect_identical(y[[1]], c("NA","NA"))
    expect_identical(y[[2]], c("A/B","NA"))
})
