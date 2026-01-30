
test_that("unit tests for samp_com_for_birth", {
  
  set.seed(99)
  meta <- rand_meta(3, 3, 25)
  
  t_obj <- samp_com_for_birth(meta)
  
  expect_length(t_obj, 1)
  expect_equal(t_obj, 2)
  expect_error(samp_com_for_birth(c(0, 3, 2)),
               regexp = "Input must be an abundance")
  
})