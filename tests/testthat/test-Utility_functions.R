
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
  
  input_1 = matrix(LETTERS, nrow=2)
  
  input_2 = list(a=letters, b=1:4)
  
  input_3 = plot(1:4, 1:4*2)
  
  expect_error(normalize_meta(input_1), class='matrix_error_incompatible_type')
  expect_error(normalize_meta(input_2), class='matrix_error_incompatible_type')
  expect_error(normalize_meta(input_3), class='matrix_error_uncoercible')
  
})