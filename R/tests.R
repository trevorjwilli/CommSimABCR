# devtools::install_github('trevorjwilli/CommSimABCR')
# install.packages('~/dev/r/CommSimABCR/')

library(CommSimABC)

meta <- rand_meta(N = 2, S = 3, 20)
xy <- matrix(sample(50, 4, replace = TRUE), ncol = 2)

inparam <- make_params(3, 2)
inparam$s <- matrix(c(1, .95, .99,
                      .95, 1, 1), byrow=T, nrow = 2)
inparam$fd <- set_fd(inparam, -.5, -.2)
inparam$mig <- set_mig(inparam, xy, 100, .05)
inparam
test <- moran_deme(meta, 100, inparam)
plot(test)

meta
