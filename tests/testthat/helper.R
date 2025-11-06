
create_balanced_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(3, 4, 20)
}

create_unbalanced_meta <- function() {
  set.seed(42)
  CommSimABC::rand_meta(3, 4, J=c(15, 20, 10))
}