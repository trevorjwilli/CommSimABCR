

test_that("unit tests for samp_com_for_birth", {
  
  meta <- create_balanced_meta()
  
  t_obj <- samp_com_for_birth(meta)
  
  expect_length(t_obj, 1)
  expect_equal(t_obj, 2)
  expect_error(samp_com_for_birth(c(0, 3, 2)),
               regexp = "Input must be an abundance")
  
})

test_that("unit tests for samp_com_birth_cpp", {
  
  meta <- create_balanced_meta()
  
  t_obj <- samp_com_for_birth_cpp(meta)
  
  expect_length(t_obj, 1)
  expect_equal(t_obj, 1)
})

test_that("unit tests for calculate coefficient", {
  tol = 0.00001
  meta <- create_balanced_meta()
  params <- create_params()
  
  cf <- calculate_species_probs(meta, 2, params)
  expected <- c(0.476443, 0.000000, 0.523557, 0.000000)
  expect_equal(cf, expected, tolerance = tol)
  
  meta <- create_20_30_meta()
  params <- create_20_30_params()
  cf <- calculate_species_probs(meta, 7, params)
  expected <- c(0.00000000,0.00000000,0.01663403,0.00000000,0.06288491,0.03367888,0.00000000,0.08218243,0.04739270,
                0.01701057,0.06683575,0.08273971,0.00000000,0.01462438,0.02825836,0.11000528,0.09154158,0.01342758,
                0.03252077,0.00000000,0.00000000,0.06208842,0.05059997,0.01722286,0.04249182,0.01307424,0.02939186,
                0.01556106,0.04097731,0.02885554)

  expect_equal(cf, expected, tolerance = tol)
})

test_that("unit tests for calculate_species_probs_cpp", {
  tol = 0.00001
  meta <- create_balanced_meta()
  params <- create_params()
  
  cf <- calculate_species_probs_cpp(meta, 1, params)
  expected <- c(0.476443, 0.000000, 0.523557, 0.000000)
  expect_equal(cf, expected, tolerance = tol)
  
  meta <- create_20_30_meta()
  params <- create_20_30_params()
  cf <- calculate_species_probs_cpp(meta, 6, params)
  expected <- c(0.00000000,0.00000000,0.01663403,0.00000000,0.06288491,0.03367888,0.00000000,0.08218243,0.04739270,
                0.01701057,0.06683575,0.08273971,0.00000000,0.01462438,0.02825836,0.11000528,0.09154158,0.01342758,
                0.03252077,0.00000000,0.00000000,0.06208842,0.05059997,0.01722286,0.04249182,0.01307424,0.02939186,
                0.01556106,0.04097731,0.02885554)
  
  expect_equal(cf, expected, tolerance = tol)
})

test_that("birth_death_process unit tests", {
  x <- create_unbalanced_meta()
  inparam <- create_params()
  
  set.seed(42)
  b_com <- samp_com_for_birth(x)
  s_probs <- calculate_species_probs(x, b_com, inparam)
  s_birth <- samp_species_for_birth(x, s_probs)
  d_com <- samp_com_for_death(x, s_birth, b_com, inparam)
  s_death <- samp_species_for_death(x, d_com)
  x <- update_meta(x, d_com, s_birth, s_death)
  
  y <- create_unbalanced_meta()
  
  set.seed(42)
  y <- birth_death_process(y, inparam)
  
  expect_equal(x, y)
  
  z <- create_unbalanced_meta()
  
  set.seed(42)
  b_com <- samp_com_for_birth_cpp(z)
  s_probs <- calculate_species_probs_cpp(z, b_com, inparam)
  s_birth <- samp_species_for_birth_cpp(z, s_probs)
  d_com <- samp_com_for_death_cpp(z, s_birth, b_com, inparam)
  s_death <- samp_species_for_death_cpp(z, d_com)
  z <- update_meta_cpp(z, d_com, s_birth, s_death)
  
  expect_equal(x, z)
  
})

test_that("birth_death_process integration test", {
  
  x <- create_unbalanced_meta()
  inparam <- create_params()
  
  birth_com <- c()
  species_birth <- c()
  death_com <- c()
  species_death <- c()
  
  for(i in 1:1000) {
    b_com <- samp_com_for_birth(x)
    birth_com[i] <- b_com
    s_probs <- calculate_species_probs(x, b_com, inparam)
    s_birth <- samp_species_for_birth(x, s_probs)
    species_birth[i] <- s_birth
    d_com <- samp_com_for_death(x, s_birth, b_com, inparam)
    death_com[i] <- d_com
    s_death <- samp_species_for_death(x, d_com)
    species_death[i] <- s_death
  }
  
  df <- data.frame(birth_com, species_birth, death_com, species_death)
  
  obs_birth_probs <- as.vector(table(df$birth_com)/length(df$birth_com))
  exp_birth_probs <- rowSums(x)/sum(x)
  
  expect_equal(obs_birth_probs, exp_birth_probs, tolerance = 0.08)
  
  c1_s_probs <- calculate_species_probs(x, 1, inparam)
  c2_s_probs <- calculate_species_probs(x, 2, inparam)
  c3_s_probs <- calculate_species_probs(x, 3, inparam)
  
  s_probs_test <- df |>
    dplyr::group_by(birth_com, species_birth) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    tidyr::complete(birth_com, species_birth, fill=list(n=0)) |>
    dplyr::group_by(birth_com) |>
    dplyr::mutate(obs_prob = n/sum(n)) |>
    dplyr::ungroup()
  
  s_probs_test$expected_prob <- c(c1_s_probs, c2_s_probs, c3_s_probs) 
  expect_equal(s_probs_test$obs_prob, s_probs_test$expected_prob, tolerance=0.08)
  
  obs_mig_df <- df |>
    dplyr::group_by(birth_com, species_birth, death_com) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    dplyr::group_by(birth_com, species_birth) |>
    dplyr::mutate(obs_probs = n/sum(n)) |>
    dplyr::ungroup()
  
  flatten_mig_mat <- function(x) {
    out <- c()
    for(j in 1:ncol(x)) {
      for(i in 1:nrow(x)) {
        out <- c(out, x[i,j])
      }
    }
    out
  }
  
  mig_list <- list()
  for(i in 1:length(inparam$mig)) {
    exp_probs <- flatten_mig_mat(inparam$mig[[i]])
    species_birth <- rep(i, length(inparam$mig[[i]]))
    birth_com <- rep(1:nrow(inparam$mig[[i]]), each=nrow(inparam$mig[[i]]))
    death_com <- rep(1:nrow(inparam$mig[[i]]), times=nrow(inparam$mig[[i]]))
    mig_list[[i]] <- data.frame(species_birth, birth_com, death_com, exp_probs)
  }
  
  mig_df <- do.call(rbind, mig_list)
  
  obs_mig_df <- obs_mig_df |>
    dplyr::left_join(mig_df)
  
  expect_equal(obs_mig_df$obs_probs, obs_mig_df$exp_probs, tolerance=0.05)
})

test_that("speciation unit tests", {
  x <- create_20_30_meta()
  inparam <- create_20_30_params()
  
  x[,4] <- 0
  x[,27] <- 0
  
  new_data <- speciate(x, inparam)
  
  expect_true(all(colSums(new_data[[1]]) > 0))
  expect_equal(rowSums(new_data[[1]]), rowSums(x))
  expect_equal(sum(new_data[[2]]$s[,4] == inparam$s[,4]), 0)
  expect_equal(sum(new_data[[2]]$s[,5] == inparam$s[,5]), nrow(x))
  expect_equal(sum(new_data[[2]]$s[,27] == inparam$s[,27]), 0)
  expect_equal(sum(new_data[[2]]$fd == inparam$fd), length(inparam$fd) - 2)
  
})

test_that("speciate_cpp unit tests", {
  x <- create_20_30_meta()
  inparam <- create_20_30_params()
  
  x[,4] <- 0
  x[,27] <- 0
  
  x2 <- create_20_30_meta()
  inparam2 <- create_20_30_params()
  
  new_data <- speciate_cpp(x, inparam, 0.1, 0.1)
  expect_true(all(colSums(new_data[[1]]) > 0))
  expect_equal(sum(new_data[[2]]$s[,4] == inparam2$s[,4]), 0)
  expect_equal(sum(new_data[[2]]$s[,5] == inparam2$s[,5]), nrow(x2))
  expect_equal(sum(new_data[[2]]$s[,27] == inparam2$s[,27]), 0)
  expect_equal(sum(new_data[[2]]$fd == inparam2$fd), length(inparam2$fd) - 2)
  
  x2[,4] <- 0
  x2[,27] <- 0
  
  expect_equal(rowSums(new_data[[1]]), rowSums(x2))
  
})
