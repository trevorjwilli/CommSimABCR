
test_that("d_alpha produces expected output", {
  tol = 0.00001
  # With balanced data and no weights
  balanced_meta <- create_balanced_meta()
  
  expect_equal(d_alpha(balanced_meta, q=0), 2.333333, tolerance = tol)
  expect_equal(d_alpha(balanced_meta, q=1), 2.21388, tolerance = tol)
  expect_equal(d_alpha(balanced_meta, q=2), 2.120141, tolerance = tol)
  expect_equal(d_alpha(balanced_meta, q=3), 2.045507, tolerance = tol)
  
  # With unbalanced data and no weights
  unbalanced_meta <- create_unbalanced_meta()
  
  expect_equal(d_alpha(unbalanced_meta, q=0), 2.666667, tolerance = tol)
  expect_equal(d_alpha(unbalanced_meta, q=1), 2.521085, tolerance = tol)
  expect_equal(d_alpha(unbalanced_meta, q=2), 2.379903, tolerance = tol)
  expect_equal(d_alpha(unbalanced_meta, q=3), 2.248595, tolerance = tol)
  
  # With unbalanced data and weights
  expect_equal(d_alpha(unbalanced_meta, q=0, weight=TRUE), 2.509836, tolerance = tol)
  expect_equal(d_alpha(unbalanced_meta, q=1, weight=TRUE), 2.509836, tolerance = tol)
  expect_equal(d_alpha(unbalanced_meta, q=2, weight=TRUE), 2.509836, tolerance = tol)
  expect_equal(d_alpha(unbalanced_meta, q=3, weight=TRUE), 2.509836, tolerance = tol)
  
  # That balanced data with and without weights are the same
  expect_equal(d_alpha(balanced_meta, q=1), d_alpha(balanced_meta, q=1, weight = TRUE), tolerance = tol)
})

test_that("d_gamma produces expected output", {
  tol = 0.00001
  # With balanced data and no weights
  balanced_meta <- create_balanced_meta()
  
  expect_equal(d_gamma(balanced_meta, q=0), 4, tolerance = tol)
  expect_equal(d_gamma(balanced_meta, q=1), 3.653728, tolerance = tol)
  expect_equal(d_gamma(balanced_meta, q=2), 3.41556, tolerance = tol)
  expect_equal(d_gamma(balanced_meta, q=3), 3.260196, tolerance = tol)
  
  # With unbalanced data and no weights
  unbalanced_meta <- create_unbalanced_meta()
  
  expect_equal(d_gamma(unbalanced_meta, q=0), 4, tolerance = tol)
  expect_equal(d_gamma(unbalanced_meta, q=1), 3.916231, tolerance = tol)
  expect_equal(d_gamma(unbalanced_meta, q=2), 3.839772, tolerance = tol)
  expect_equal(d_gamma(unbalanced_meta, q=3), 3.770926, tolerance = tol)
})

test_that("d_beta produces expected output", {
  tol = 0.00001
  # With balanced data and no weights
  # With balanced data and no weights
  balanced_meta <- create_balanced_meta()
  
  expect_equal(d_beta(balanced_meta, q=0), 1.714286, tolerance = tol)
  expect_equal(d_beta(balanced_meta, q=1), 1.650373, tolerance = tol)
  expect_equal(d_beta(balanced_meta, q=2), 1.611006, tolerance = tol)
  expect_equal(d_beta(balanced_meta, q=3), 1.593832, tolerance = tol)
  
  # With unbalanced data and no weights
  unbalanced_meta <- create_unbalanced_meta()
  
  expect_equal(d_beta(unbalanced_meta, q=0), 1.5, tolerance = tol)
  expect_equal(d_beta(unbalanced_meta, q=1), 1.553391, tolerance = tol)
  expect_equal(d_beta(unbalanced_meta, q=2), 1.613416, tolerance = tol)
  expect_equal(d_beta(unbalanced_meta, q=3), 1.677014, tolerance = tol)
  
  # With unbalanced data and weights
  expect_equal(d_beta(unbalanced_meta, q=0, weight=TRUE), 1.579678, tolerance = tol)
  expect_equal(d_beta(unbalanced_meta, q=1, weight=TRUE), 1.579678, tolerance = tol)
  expect_equal(d_beta(unbalanced_meta, q=2, weight=TRUE), 1.579678, tolerance = tol)
  expect_equal(d_beta(unbalanced_meta, q=3, weight=TRUE), 1.579678, tolerance = tol)
  
})



