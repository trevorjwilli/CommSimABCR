
create_balanced_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(3, 4, 20)
}

create_unbalanced_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(3, 4, J=c(15, 20, 10))
}

create_params <- function() {
  set.seed(33)
  xy <- matrix(sample(25, 6, replace = TRUE), ncol = 2)
  
  inparam <- make_params(4,3)
  inparam$s <- set_sel(inparam, 'beta', 10, 1)
  inparam$fd <- set_fd(inparam, -1, -.1)
  inparam$mig <- set_mig(inparam, xy, 25, .05)
  inparam
}

create_20_30_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(20, 30, J=sample(20:100, size=20))
}

create_20_30_metas <- function() {
  set.seed(42)
  replicate(20, CommSimABC::rand_meta(20, 30, J=sample(20:100, size=20)),
            simplify = FALSE)
}

create_20_30_params <- function() {
  set.seed(33)
  xy <- matrix(sample(100, 40, replace = TRUE), ncol = 2)
  
  inparam <- make_params(30, 20)
  inparam$s <- set_sel(inparam, 'beta', 10, 1)
  inparam$fd <- set_fd(inparam, -1, -.1)
  inparam$mig <- set_mig(inparam, xy, 50, .05)
  inparam
}