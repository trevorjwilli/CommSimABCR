
create_balanced_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(3, 4, 20)
}

create_unbalanced_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(3, 4, J=c(15, 20, 10))
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