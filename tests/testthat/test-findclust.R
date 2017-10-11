context("find.clusters tests")

test_that("find.clusters works with pre-defined clusters", {
  skip_on_cran()
  data(nancycats)
  
  # set connection
  f <- file()
  options(adegenet.testcon = f)
  
  # add data to connection
  twos <- paste(rep(2, nPop(nancycats)), collapse = "\n")
  write(twos, f)
  
  # test function
  expect_warning(capture.output(res <- find.clusters(nancycats, clust = pop(nancycats), n.pca = 100)))
  # We expect each group to receive two clusters
  expect_equal(length(levels(res$grp)), nPop(nancycats) * 2)
  # We expect all individuals accounted for
  expect_equal(length(res$grp), nInd(nancycats))
  
  # reset variable
  options(adegenet.testcon = stdin())
  # close connection
  close(f)
})

