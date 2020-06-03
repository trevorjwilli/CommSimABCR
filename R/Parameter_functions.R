### Functions for initializing parameters for moran_deme()

#' Create params object.
#'
#' Initializes an object of class \emph{params}.
#'
#' @param n.spec Numeric, the number of species in the metacommunity.
#' @param n.site Numeric, the number of communities in the metacommunity.
#'
#' @details This function is used to initialize a \emph{params} object.
#'   A \emph{params} object is a list with 3 named elements: s, fd, and mig.
#'   s is a matrix of selection coefficients where columns correspond to species
#'   and rows correspond to sites (or communities). fd is a vector containing
#'   the frequency dependence parameters for each species. mig is a list of
#'   migration matrices, where cell ij is the probability for each species to
#'   migrate from site j to site i.
#'
#' @return Returns an empty object of S3 class \emph{params}.
#'
#' @examples
#' make_params(5,5)
#'
#' @export

make_params <- function(n.spec, n.site) {
  x <- list(s = matrix(nrow = n.site, ncol = n.spec), fd = vector(length = n.spec), mig = NULL)
  attr(x, 'NumSite') <- n.site
  attr(x, 'NumSpec') <- n.spec
  class(x) <- 'params'
  x
}

#' Random selection matrix.
#'
#' Function to create random selection matrix using uniform, normal, or gamma distributions.
#'
#' @param paramfile \emph{params} object initialized using \code{\link{make_params}}.
#' @param distr Character, one of 'uniform', 'normal', or 'gamma'.
#' @param input1,input2 Numeric, parameters for distribution being used. See details.
#'
#' @details This function generates a random selection matrix according to a specified
#'   distribution. The specified distribution can be one of 'uniform', 'normal', or 'gamma.'
#'   If 'uniform' input1 is the minimum value and input2 is the maximum value. If 'normal'
#'   input1 is the mean and input2 is the standard deviation. If 'gamma' input1 is the
#'   shape and input2 is the scale parameters.
#'
#' @return Returns a matrix of selection coefficients where cell ij is the selection
#' coefficient for species j in community i.
#'
#' @examples
#'
#' paramfile <- make_params(5,5)
#' paramfile$s <- set_sel(paramfile, 'normal', .5, .2)
#'
#' ## For neutral model ##
#' paramfile <- make_params(5,5)
#' paramfile$s <- set_sel(paramfile, 'uniform', 1,1)
#'
#' @export

set_sel <- function(paramfile, distr, input1, input2) {
  if(!is.params(paramfile)) {
    stop("paramfile not a params class")
  }

  n.spec <- attr(paramfile, 'NumSpec') # Calculate number of species

  n.sites <- attr(paramfile, 'NumSite') # Calculate number of communities

  if(distr == 'uniform') { # If prior for selection is uniform
    ##print("uniform")
    sel <- matrix(ncol = n.spec, nrow = n.sites)
    if(input1 < 0 | input2 < 0) {
      stop('Selection coefficient cannot be less than 0')
    }
    for(j in 1:n.sites) {
      coeffs <- runif(n.spec, min = input1, max = input2)
      sel[j,] <- coeffs
    }
    ##print(sel)

  } else if(distr == 'normal') { # If prior for selection is normal
    ##print("normal")
    sel <- matrix(ncol = n.spec, nrow = n.sites)
    if(input1 < 0) {
      stop('Mean selection coefficient cannot be less than zero')
    }
    for(j in 1:n.sites) {
      coeffs <- rnorm(n.spec, mean = input1, sd = input2)
      coeffs[which(coeffs < 0)] <- 0.00000001 # Make sure there are no 0s
      sel[j,] <- coeffs
    }
    ##print(sel)

  } else if(distr == 'gamma') { # If prior for selection is gamma distributed
    ##print("gamma")
    sel <- matrix(ncol = n.spec, nrow = n.sites)
    if(input1 < 0 | input2 < 0) {
      stop('Shape and scale parameters must be positive')
    }
    for(j in 1:n.sites) {
      coeffs <- rgamma(n.spec, shape = input1, scale = input2)
      coeffs[which(coeffs < 0)] <- 0.00000001 # Make sure there are no 0s
      sel[j,] <- coeffs
    }
    ##print(sel)
  }
  return(sel)
}

#' Random frequency dependence parameters.
#'
#' Function to create random vector of frequency dependence parameters.
#'
#' @param paramfile \emph{params} object initialized using \code{\link{make_params}}.
#' @param input1,input2 Numeric, parameters for uniform distribution.
#'
#' @details This function randomly generates frequency dependence parameters by sampling
#' from a random uniform distribution where input1 is the minimum and input2 is the
#' maximum.
#'
#' @return Returns a numeric matrix of selection coefficients where cell ij is the selection
#' coefficient for species j in community i.
#'
#' @examples
#' paramfile <- make_params(5,5)
#' paramfile$fd <- set_fd(paramfile, -1, 0)
#'
#' @export

set_fd <- function(paramfile, input1, input2) {
  if(!is.params(paramfile)) {
    stop("paramfile not a params class")
  }

  n.spec <- attr(paramfile, 'NumSpec') # Calculate number of species
  fd <- runif(n.spec, input1, input2)
  return(fd)
}

#' Create migration matrices.
#'
#' Function to create list of migration matrices.
#'
#' @param paramfile \emph{params} object initialized using \code{\link{make_params}}.
#' @param site.arrange Matrix or Distance object. Either a two column matrix with x-values
#'   \(i.e. longitude\) in the first column and y-values (i.e. latitude) in the second
#'   column or a pairwise distance object between each site.
#' @param max.dist Numeric vector containing the maximum distance species can migrate.
#' @param tot Numeric vector containing values between 0 and 1 of the probability
#'   for a species to migrate.
#'
#' @details The parameters max.dist and tot can be either single values or vectors. If
#'   a single value the code assumes that each species has the same migration matrix,
#'   if a vector, it creates separate migration matrices for each species according
#'   to the inputed parameters. Species specific parameters for max.dist and tot must
#'   be in the same order as in the metacommunity matrix.
#'
#' @return Returns a list of numeric migration matrices.
#'
#' @seealso \code{\link{make_mig}}
#'
#' @examples
#' paramfile <- make_params(5,5)
#'
#' # Make random landscape with 5 sites
#' xy <- matrix(sample(50, 10), ncol = 2)
#'
#' ## Where all species have the same migration matrix ##
#'
#' paramfile$mig <- set_mig(paramfile, xy, 10, .1)
#'
#' ## Where each species has own migration matrix
#'
#' paramfile$mig <- set_mig(paramfile, xy, c(10, 9, 15, 2, 11), c(.1, .2, .05, .1, .15))
#'
#' @export

set_mig <- function(paramfile, site.arrange, max.dist, tot) {
  if(!is.params(paramfile)) {
    stop("paramfile not a params class")
  }

  if(class(site.arrange) == "dist") {
    n.sites <- attr(site.arrange, "Size")
  } else {
    n.sites <- (length(site.arrange[,1]))
  }

  if(n.sites != attr(paramfile, 'NumSite')) {
    stop('Number of sites does not equal parameter file site number')
  }

  n.spec <- attr(paramfile, 'NumSpec')

  if(length(max.dist) == 1 & length(tot == 1)) {
    mig <- make_mig(site.arrange, max.dist, tot)
    mig.out <- replicate(n.spec, mig, simplify = F)
  }

  else if(length(max.dist) == n.spec & length(tot) == n.spec) {
    mig.out <- list()
    for(i in 1:n.spec) {
      #print(max.dist[i])
      #print(tot[i])
      mig <- make_mig(site.arrange, max.dist[i], tot[i])
      mig.out[[i]] <- mig
    }
  } else {
    stop("Number of distance and probability parameters must equal the number of species \n  or have length of 1")
  }

  return(mig.out)
}


#' @export

print.params <- function(x) {
  cat('\n')
  cat('Number of Sites:', attr(x, 'NumSite'), '\n')
  cat('Number of Species:', attr(x, 'NumSpec'), '\n\n')

  if(is.null(rownames(x$s))) {
    rownames(x$s) <- paste0("Site", 1:attr(x, 'NumSite'))
  }

  if(is.null(colnames(x$s))) {
    colnames(x$s) <- paste0("Species", LETTERS[1:attr(x, 'NumSpec')])
  }

  cat('Selection Matrix:\n\n')
  print(x$s, quote = F)
  cat('\n')

  if(is.null(names(x$fd))) {
    names(x$fd) <- paste0("Species", LETTERS[1:attr(x, 'NumSpec')])
  }

  cat('Frequency Dependence Parameters:\n\n')
  print(x$fd, quote = F)
  cat('\n')

  for(i in 1:attr(x, 'NumSpec')) {
    cat('Migration Matrix for Species', LETTERS[i], '\n\n')
    print(x$mig[[i]], quote = F)
    cat('\n')
  }
}

#' @export

is.params <- function(x) inherits(x, "params")





