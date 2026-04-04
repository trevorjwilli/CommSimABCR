
test_that("run_single_sim params are correct", {
  tol = 0.00001
  # Make Neutral priors object
  NT_dl_priors <- make_priors(n.spec = 30, n.site = 20)
  
  NT_dl_priors$jdist <- 1
  NT_dl_priors$jparams <- c(100, 1000)
  NT_dl_priors$seldist <- 1
  NT_dl_priors$selparams <- c(1, 1)
  NT_dl_priors$fddist <- 1
  NT_dl_priors$fdparams <- c(0, 0)
  NT_dl_priors$migdist <- 1
  NT_dl_priors$migdistparams <- c(10, 50)
  NT_dl_priors$migprobdist <- 1
  NT_dl_priors$migprobparams <- c(0.05, 0.2)
  
  NT_test <- run_single_sim(t=10, priors = NT_dl_priors,
                            eqpop = TRUE, eqmig = TRUE, output = TRUE)
  
  # Check that selection coefficient matrix is correct size
  testthat::expect_equal(dim(NT_test$input[[2]]$s), c(20, 30))
  # Check that selection coefficients are all 1
  testthat::expect_equal(NT_test$input[[2]]$s, matrix(1, nrow=20, ncol=30))
  # Check that fd vector is right size
  testthat::expect_length(NT_test$input[[2]]$fd, 30)
  # Check that fd vector is all 0s
  testthat::expect_equal(NT_test$input[[2]]$fd, rep(0, 30))
  # Check that migration matrix list is correct size
  testthat::expect_length(NT_test$input[[2]]$mig, 30)
  # Check that all migration matrices are the same
  testthat::expect_length(unique(NT_test$input[[2]]$mig), 1)
  # Check that communities all have same size
  testthat::expect_length(unique(rowSums(NT_test$metacommunity)), 1)
  
  
  # Make Species sorting priors object
  SS_ndl_priors <- make_priors(n.spec = 30, n.site = 20)
  
  SS_ndl_priors$jdist <- 1
  SS_ndl_priors$jparams <- c(100, 1000)
  SS_ndl_priors$seldist <- 2
  SS_ndl_priors$selparams <- c(0.5, 0.25)
  SS_ndl_priors$fddist <- 1
  SS_ndl_priors$fdparams <- c(-1, 0)
  SS_ndl_priors$migdist <- 1
  SS_ndl_priors$migdistparams <- c(100, 200)
  SS_ndl_priors$migprobdist <- 1
  SS_ndl_priors$migprobparams <- c(0.05, 0.2)
  
  SS_test <- run_single_sim(t=10, priors = SS_ndl_priors,
                            eqpop = FALSE, eqmig = FALSE,
                            output = TRUE)
  
  # Check that selection coefficients are not all 1
  testthat::expect_false(sum(SS_test$input[[2]]$s) == sum(matrix(1, nrow=20, ncol=30)))
  # Check that fd coefficients are 0 or less
  testthat::expect_equal(all(SS_test$input[[2]]$fd <= 0), TRUE)
  # Check that all migration matrices are the same
  testthat::expect_gt(length(unique(SS_test$input[[2]]$mig)), 1)
  # Check that communities all have same size
  testthat::expect_gt(length(unique(rowSums(SS_test$metacommunity))), 1)
})
