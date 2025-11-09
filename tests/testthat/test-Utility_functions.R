
test_that("normalize_meta outputs are correct", {
  
  input = matrix(c(10, 0, 5, 5,
                   0, 20, 0, 0, 
                   8, 2, 8, 2),
                 nrow = 3, byrow = T)
  
  output = matrix(c(0.5, 0.0, 0.25, 0.25,
                    0.0, 1.0, 0.0, 0.0, 
                    0.4, 0.1, 0.4, 0.1),
                  nrow = 3, byrow = T)
  
  # Output is as expected
  expect_equal(normalize_meta(input), output)
  # Rows sum to 1
  expect_equal(rowSums(normalize_meta(input)), rep(1,nrow(input)))
  
})


test_that("normalize_meta warnings are correct", {
  
  input = matrix(c(10, 0, 5, 5,
                   0, 0, 0, 0, 
                   8, 2, 8, 2),
                 nrow = 3, byrow = T)
  
  expect_warning(normalize_meta(input), "Removing empty sites")
  
})


test_that("normalize_meta_errors are correct", {
  
  input_1 <- matrix(LETTERS, nrow=2)
  
  input_2 <- list(a=letters, b=1:4)
  
  input_3 <- plot(1:4, 1:4*2)
  
  expect_error(normalize_meta(input_1), class='matrix_error_incompatible_type')
  expect_error(normalize_meta(input_2), class='matrix_error_incompatible_type')
  expect_error(normalize_meta(input_3), class='matrix_error_uncoercible')
  
})

test_that("ecosumstats outputs are correct", {
  tol <- 0.000001
  meta <- create_20_30_meta()
  metas <- create_20_30_metas()
  
  test_output <- c(2.28059426390603, 2.36259551061325, 0.703920694658288,
                   0.857559875810992, 0.89716049382716, 0.131219473180127,
                   0.961790528427438, 0.961049262043089, 0.0209900987734729,
                   13.45, 12.5, 7.85711209836768, 8.96666666666667, 10,
                   2.00831604418561, 0.762858754112658, 0.777404718693285,
                   0.149547933873534, 0.302100177456275, 13.45, 2.23048327137546,
                   30, 9.78249206447442, 2.82954458718889, 27.6799974702518,
                   7.02049373863978, 3.68962227805471, 25.902970101029,
                   10.1150551032724, 2.70031439531347, 27.3138289047555,
                   30.2631578947368, 0.770036327580812, 481.25313283208,
                   12, 1.07014333638304, 232, 2015, 1.27040816326531)
  
  names(test_output) <- c("mean.shannon", "median.shannon", "sd.shannon",
                          "mean.simpson", "median.simpson", "sd.simpson",
                          "mean.evenness", "median.evenness", "sd.evenness",
                          "mean.richness", "median.richness", "sd.richness",
                          "mean.sites", "median.sites", "sd.sites", "beta.bray",
                          "beta.bray.med", "beta.bray.sd", "beta.bray.ssid",
                          "alpha.0", "beta.0", "gamma.0", "alpha.1", "beta.1",
                          "gamma.1", "alpha.2", "beta.2", "gamma.2", 
                          "alpha.1.weight", "beta.1.weight", "gamma.1.weight",
                          "c.score", "c.score.skew", "c.score.var",
                          "checkerscore", "v.ratio", "coherence", "turnov",
                          "morisit")
  
  
  # Test that output is correct when passing in one matrix
  expect_warning(expect_equal(ecosumstats(meta), test_output))
  expect_warning(expect_true(is.vector(ecosumstats(meta))))
  
  # Test that dimensions are correct when passing in list of matrices
  expect_warning(expect_equal(dim(ecosumstats(metas)), c(20, 39)))
  

})