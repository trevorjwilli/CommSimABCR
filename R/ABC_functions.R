### ABC functions ###

#' Initialize \emph{priors} object
#'
#' Creates an empty priors object for use in \code{\link{abc_moran_deme}}
#'
#' @param n.spec Numeric, the number of species in the metacommunity
#' @param n.site Numeric, the number of communities in the metacommunity
#'
#' @return Creates a \emph{priors} object which is a list containing 10 elements.
#'   jdist: One of 1 or 2, specifies a uniform or normal prior distribution for the community size.
#'
#'   jparams: A 2L Numeric vector. If jdist = 1 then the first value is the minimum and the
#'   second value is the maximum. If jdist = 2 then the first value is the mean and the second
#'   value is the standard deviation.
#'
#'   seldist: One of 1, 2, 3, or 4, specifies a uniform, normal, gamma, or beta prior distribution for the selection
#'   coefficients for each community.
#'
#'   selparams: A 2L Numeric vector. If seldist = 1 then the first value is the minimum and the
#'   second value is the maximum. If seldist = 2 then the first value is the mean and the second
#'   value is the standard deviation. If seldist = 3 then the first value is the shape parameter
#'   and the second value is the scale parameter. If seldist = 4 then the first value is the first
#'   shape parameter and the second value is the second shape parameter.
#'
#'   fddist: One of 1 or 2, specifies a uniform or normal prior distribution for the frequency dependence parameters.
#'
#'   fdparams: A 2L Numeric vector. If fddist = 1 then the first value is the minimum and the
#'   second value is the maximum. If fddist = 2 then the first value is the mean and the second
#'   value is the standard deviation.
#'
#'   migdist: One of 1 or 2, specifies a uniform or normal prior distribution for the migration distance.
#'   Can be a vector with length equal to the number of species in simulations to give each species its own
#'   migration distance prior distribution
#'
#'   migdistparams: A 2L Numeric vector. If migdist = 1 then the first value is the minimum and the
#'   second value is the maximum. If migdist = 2 then the first value is the mean and the second
#'   value is the standard deviation. Can also be a list of 2L numeric vectors specifying prior distribution
#'   parameters for each species individually.
#'
#'   migprobdist: One of 1 or 2, specifies a uniform or normal prior distribution for the migration probability.
#'   Can be a vector with length equal to the number of species in simulations to give each species its own
#'   migration probability prior distribution
#'
#'   migprobparams: A 2L Numeric vector. If migprobdist = 1 then the first value is the minimum and the
#'   second value is the maximum. If migprobdist = 2 then the first value is the mean and the second
#'   value is the standard deviation. Can also be a list of 2L numeric vectors specifying prior distribution
#'   parameters for each species individually.
#'
#' @examples
#' testpriors <- make_priors(5, 5)
#'
#' testpriors$jdist <- 1
#' testpriors$jparams <- c(100, 200)
#' testpriors$seldist <- 2
#' testpriors$selparams <- c(.5, .1)
#' testpriors$fddist <- 1
#' testpriors$fdparams <- c(-.5, -.2)
#' testpriors$migdist <- 1
#' testpriors$migdistparams <- c(10, 50)
#' testpriors$migprobdist <- 1
#' testpriors$migprobparams <- c(.1, .2)
#'
#' @export


make_priors <- function(n.spec, n.site) {
  x <- list(jdist = NULL,
            jparams = vector(length = 2),
            seldist = NULL,
            selparams = vector(length = 2),
            fddist = NULL,
            fdparams = vector(length = 2),
            migdist = NULL,
            migdistparams = vector(length = 2),
            migprobdist = NULL,
            migprobparams = vector(length = 2))
  attr(x, 'NumSite') <- n.site
  attr(x, 'NumSpec') <- n.spec
  class(x) <- 'priors'
  x
}

#' Create Selection Matrix from Prior
#'
#' Creates a random selection matrix according to either a Uniform, Normal,
#' or Gamma prior distribution
#'
#' @param n.spec Numeric, the number of species in the metacommunity
#' @param n.sites Numeric, the number of communities in the metacommunity
#' @param distr Either a 1, 2, 3, or 4 corresponding to Uniform, Normal,
#' Gamma, and Beta prior distributions respectively.
#' @param input1,input2 Numeric, parameters for prior distributions. See details.
#'
#' @details Creates a random selection coeffient matrix according to a prior
#' distribution. Selection coefficients are randomly selected for species in
#' each community according to the prior. If distr = 1, then input 1 is the
#' minimum and input 2 is the maximum. If distr = 2, then input 1 is the mean
#' and input 2 is the standard deviation. If distr = 3, then input 1 is the shape
#' parameter, and input 2 is the scale parameter. If seldist = 4 then the first 
#' value is the first shape parameter and the second value is the second shape parameter.
#'
#' @return A selection matrix where columns are species, rows are communities and
#' cell ij is the selection coefficient for species j in community i.
#'
#' @examples
#'
#' # For Uniform distribution
#' set_sel_priors(5, 5, distr = 1, input1 = .1, input2 = 1)
#'
#' # For Normal distribution
#' set_sel_priors(5, 5, distr = 2, input1 = .5, input2 = .2)
#'
#' # For Gamma Distribution
#' set_sel_priors(5, 5, distr = 3, input1 = 2, input2 = .2)
#' 
#' # For Beta distribution
#' set_sel_priors(5, 5, distr = 4, input1 = 8, input2 = 1)
#'
#' @export

set_sel_priors <- function(n.spec, n.sites, distr, input1, input2) {

  sel <- matrix(ncol = n.spec, nrow = n.sites)

  if(distr == 1) { # If prior for selection is uniform
    for(j in 1:n.sites) {
      coeffs <- stats::runif(n.spec, min = input1, max = input2)
      sel[j,] <- coeffs
    }

  } else if(distr == 2) { # If prior for selection is normal
    for(j in 1:n.sites) {
      coeffs <- stats::rnorm(n.spec, mean = input1, sd = input2)
      coeffs[which(coeffs < 0)] <- 0.00000001 # Make sure there are no 0s
      sel[j,] <- coeffs
    }


  } else if(distr == 3) { # If prior for selection is gamma distributed
    for(j in 1:n.sites) {
      coeffs <- stats::rgamma(n.spec, shape = input1, scale = input2)
      coeffs[which(coeffs < 0)] <- 0.00000001 # Make sure there are no 0s
      sel[j,] <- coeffs
    }
  } else if(distr == 4) { # If prior for selection is beta distributed
    for(j in 1:n.sites) {
      coeffs <- stats::rbeta(n.spec, shape1 = input1, shape2 = input2)
      sel[j,] <- coeffs
    }
  }
  return(sel)
}

#' Run Multiple Moran Community Simulations With Priors.
#'
#' This function runs multiple simulations of the Moran Community model by randomly generating parameters
#' according to user input prior distributions for use in downstream ABC analyses.
#'
#' @param nsims Numeric, The number of simulations to run.
#' @param t Numeric, the number of generations to run for each simulation.
#' @param priors Priors object specifying the priors to use for simulations. See \code{\link{make_priors}}.
#' @param x.max,y.max Numeric, maximum spatial extent in x or y direction for site placement if randomly
#' generating spatial data.
#' @param spatial Two column Numeric matrix or Distance object specifying spatial arrangement of communities.
#' @param eqpop Logical, if TRUE all community sizes are the same within a simulation, if FALSE
#' each community is randomly given a community size each simulation.
#' @param eqmig Logical, if TRUE all species will have the same migration matrix, if FALSE
#' each species will recieve it's own migration matrix each simulation according to priors file settings
#' @param outgens Integer, vector giving the generations for which metacommunities should be output
#'
#' @details This function is used to run the Moran Community model simulation multiple times in preparation
#' for Approximate Bayesian Analysis (ABC). Users specify prior distributions for community size,
#' selection coefficients, frequency dependence, and migration parameters in a \emph{priors} object. From
#' this information, parameter values are randomly drawn according to the prior distributions and simulations
#' are run from these parameter values.
#'
#' @return If all communities have the same size, then returns a list with a list of final metacommunity matrices
#' for each simulation, a list of selection matrices used for each simulation, a matrix with frequency dependence
#' parameters for each species (columns are species and rows are simulation runs), and a dataframe with the
#' parameters used for each simulation (each row corresponds to one simulation). If communities are allowed to have
#' different sizes, then the list also contains a matrix containing the community sizes for each community in each
#' simulation.
#'
#' @examples
#'
#' testpriors <- make_priors(5, 5)
#' xy <- random_points(5, 100, 100)
#'
#' testpriors$jdist <- 1
#' testpriors$jparams <- c(100, 200)
#' testpriors$seldist <- 2
#' testpriors$selparams <- c(.5, .1)
#' testpriors$fddist <- 1
#' testpriors$fdparams <- c(-.5, -.2)
#' testpriors$migdist <- 1
#' testpriors$migdistparams <- c(10, 50)
#' testpriors$migprobdist <- 1
#' testpriors$migprobparams <- c(.1, .2)
#'
#' abc_moran_deme(5, 5, testpriors, eqpop = FALSE, spatial = xy)
#'
#' @export

abc_moran_deme <- function(nsims, t, priors, x.max = NULL, y.max = NULL, spatial = NULL, eqpop = FALSE, eqmig = TRUE, outgens = NULL) {

  n.spec <- attr(priors, 'NumSpec') # Calculate number of species

  n.sites <- attr(priors, 'NumSite') # Calculate number of communities

  if(eqpop == FALSE) { # Check to see if all communities will have the same community size, if False...

    if(priors$jdist == 1) { # Check to see if the prior distribution is uniform for community size
      J <- t(replicate(nsims, round(stats::runif(n.sites, min = priors$jparams[1], max = priors$jparams[2])))) # Create a matrix of population sizes where each row is a simulation and each column is a community
    } else if(priors$jdist == 2) { # If prior distribution is normal
      J <- t(replicate(nsims, round(stats::rnorm(n.sites, mean = priors$jparams[1], sd = priors$jparams[2])))) # Create matrix
    }

    if(length(which(J <= 0)) > 0) { # Check to see if any population sizes are negative or zero
      stop("Negative values in population size. Check parameters")
    }

  } else if(eqpop == TRUE) { # If all communities are the same size...

    if(priors$jdist == 1) { # Use the uniform distribution
      J <- round(stats::runif(nsims, min = priors$jparams[1], max = priors$jparams[2])) # Create vector of community sizes, each element corresponding to a simulation
    } else if(priors$jdist == 2) { # Use the normal distribution
      J <- round(stats::rnorm(nsims, mean = priors$jparams[1], sd = priors$jparams[2])) # Create vector of community sizes
    }

    if(length(which(J <= 0)) > 0) { # Check to see if any communities have size zero or negative
      stop("Negative values in population size. Check parameters")
    }
  }

  sel <- replicate(nsims, set_sel_priors(n.spec, n.sites, priors$seldist, priors$selparams[1], priors$selparams[2]), simplify = F) # Create selection matrices for each simulation using set_sel_priors function

  if(priors$fddist == 1) { # Use Normal distribution for frequency dependence parameters
    fd <- t(replicate(nsims, stats::runif(n.spec, min = priors$fdparams[1], priors$fdparams[2]))) # Make matrix of fd parameters, each row is a simulation each column a species
  }
  else if(priors$fddist == 2) { # Use Normal distribution for frequency dependence parameters
    fd <- t(replicate(nsims, stats::rnorm(n.spec, mean = priors$fdparams[1], sd = priors$fdparams[2]))) # Make matrix of fd parameters
  }

  if(eqmig == TRUE) {
    if(priors$migdist == 1) { # Use Uniform distribution for migration distance
      max.dist <- round(stats::runif(nsims, min = priors$migdistparams[1], max = priors$migdistparams[2])) # Create vector of migration distances, each element is for one simulation
    }
    else if(priors$migdist == 2) { # Use Normal distribution for migration distance
      max.dist <- round(stats::rnorm(nsims, mean = priors$migdistparams[1], sd = priors$migdistparams[2])) # Create vector of migration distances
    }

    if(priors$migprobdist == 1) { # Use Uniform distribution for migration probability
      tot <- stats::runif(nsims, min = priors$migprobparams[1], max = priors$migprobparams[2]) # Create vector of migration probabilities, each element corresponds with a simulation
    }
    else if(priors$migprobdist == 2) { # Use Normal distribution for migration probability
      tot <- stats::rnorm(nsims, mean = priors$migprobparams[1], sd = priors$migprobparams[2]) # Create vector of migration probabilities
    }
  } else if(eqmig == FALSE) {
    if(length(priors$migdist) == 1) { # If only one input for migdist replicate for number of species
      priors$migdist <- rep(priors$migdist, attr(priors, "NumSpec"))
    }
    if(length(priors$migprobdist) == 1) { # if only one input for migprobdist, replicate for number of species
      priors$migprobdist <- rep(priors$migprobdist, attr(priors, "NumSpec"))
    }
    if(length(priors$migdistparams == 2)) {
      priors$migdistparams <- replicate(attr(priors, "NumSpec"), priors$migdistparams, simplify = F)
    }
    if(length(priors$migprobparams == 2)) {
      priors$migprobparams <- replicate(attr(priors, "NumSpec"), priors$migprobparams, simplify = F)
    }

    if(length(priors$migdist) != attr(priors, "NumSpec") | length(priors$migprobdist) != attr(priors, "NumSpec")) {
      stop("Incorrect number of migration distribution priors, check priors object")
    }

    migdistparam1 <- sapply(priors$migdistparams, "[", 1) # Make vector of all first parameter values for migration distance
    migdistparam2 <- sapply(priors$migdistparams, "[", 2) # Make vector of all second parameter values for migration distance

    migprobparam1 <- sapply(priors$migprobparams, "[", 1) # Make vector of all first parameter values for migration probability
    migprobparam2 <- sapply(priors$migprobparams, "[", 2) # Make vector of all second parameter values for migration probability

    distrmaker <- function(n, distr, param1, param2) { # Function to run multiple random number distributions according to
      if(distr == 1) {                                 # Uniform, Normal, or Gamma
        out <- stats::runif(n, min = param1, max = param2)
      } else if(distr == 2) {
        out <- stats::rnorm(n, mean = param1, sd = param2)
      }
      out
    }

    max.dist <- round(mapply(distrmaker, n = nsims, distr = priors$migdist, param1 = migdistparam1, param2 = migdistparam2)) # Create matrix of maximum distances
    tot <- mapply(distrmaker, n = nsims, distr = priors$migprobdist, param1 = migprobparam1, param2 = migprobparam2) # Create matrix of migration probabilities
    if(length(which(max.dist < 0)) > 0 | length(which(tot < 0)) > 0) {
      stop("Negative values for migration parameters, check prior settings")
    }
  }

  if(is.null(spatial)) { # Check to see if user inputed spatial information is available
    site.arrange <- random_points(n.sites, x.max, y.max) # If not randomly generate spatial data for sites (different for each simulation)
  } else { # If present
    site.arrange <- spatial # Make list of spatial data
  }

  if(is.matrix(J)) { # If different community sizes... (Needed because structure of J is different depending on this logical)
    meta <- list() # Initialize list of metacommunity inputs
    for(m in 1:nsims) { # For each simulation...
      meta[[m]] <- rand_meta(n.sites, n.spec, J = J[m,]) # Create a metacommunity
    }
  } else if(!is.matrix(J)) { # If all the same community size
    meta <- list() # Initialize list of metacommunities
    for(w in 1:nsims) { # For each simulation...
      meta[[w]] <- rand_meta(n.sites, n.spec, J = J[w]) # Create a metacommunity
    }
  }

  names(sel) <- paste0("sim", 1:nsims) # Name simulations in selection matrix list

  param.out <- list(J = J, max.dist = max.dist, mig.prob = tot, sel = sel, fd = fd) # Create list of parameter values used

  pb <- utils::txtProgressBar(min = 0, max = nsims, style = 3) # Set up progress bar

  out.meta <- list()
  out.params <- list()
  for(i in 1:nsims) { # Run simulations

    inparam <- make_params(n.spec, n.sites) # Set params object
    inparam$s <- sel[[i]] # give appropriate selection matrix
    inparam$fd <- fd[i,] # give appropriate fd vector
    if(is.matrix(max.dist) & is.matrix(tot)) {
      inparam$mig <- set_mig(inparam, site.arrange, max.dist[i,], tot[i,])
    } else {
      inparam$mig <- set_mig(inparam, site.arrange, max.dist[i], tot[i]) # Create migration matrices
    }

    run <- moran_deme(x = meta[[i]], t = t, params = inparam, output = F, outgens = outgens) # Run simulation
    out.meta[[paste0('sim',i)]] <- run$Metacommunity
    out.params[[paste0('sim',i)]] <- inparam
    utils::setTxtProgressBar(pb, i) # Update progress bar

  }
  return(list(metacommunities = out.meta, parameters = param.out, input = priors, parameterfiles = out.params, nsims = nsims, time = t))
}


#' Check if object is priors class
#'
#' Checks to see if an object is of class priors
#'
#' @param x object to test
#'
#' @return Logical
#'
#' @examples
#'
#' test <- make_priors(5,5)
#' is.priors(test)
#'
#' @export

is.priors <- function(x) inherits(x, "priors")

#' @export

print.priors <- function(x, ...) {
  cat('\n')
  cat('Number of Sites:', attr(x, 'NumSite'), '\n')
  cat('Number of Species:', attr(x, 'NumSpec'), '\n\n')

  if(is.null(x$jdist)) {
    cat('Community Size Distribution:', x$jdist, '\n')
    cat('Community Size Input 1:', x$jparams[1], '\n')
    cat('Community Size Input 2:', x$jparams[2], '\n\n')
  }
  else if(x$jdist == 1) {
    cat('Coummunity Size Distribution: Uniform\n')
    cat('Minimum Community Size:', x$jparams[1], '\n')
    cat('Maximum Community Size:', x$jparams[2], '\n\n')
  }
  else if(x$jdist == 2) {
    cat('Coummunity Size Prior Distribution: Normal\n')
    cat('Mean:', x$jparams[1], '\n')
    cat('Standard Deviation:', x$jparams[2], '\n\n')
  } else {
    cat('Community Size Distribution: ERROR\n')
    cat('Minimum Community Size: ERROR\n')
    cat('Maximum Community Size: ERROR\n\n')
  }

  if(is.null(x$seldist)) {
    cat('Selection Distribution:', x$seldist, '\n')
    cat('Selection Input 1:', x$selparams[1], '\n')
    cat('Selection Input 2:', x$selparams[2], '\n\n')
  }
  else if(x$seldist == 1) {
    cat('Selection Distribution: Uniform\n')
    cat('Minimum:', x$selparams[1], '\n')
    cat('Maximum:', x$selparams[2], '\n\n')
  }
  else if(x$seldist == 2) {
    cat('Selection Distribution: Normal\n')
    cat('Mean:', x$selparams[1], '\n')
    cat('Standard Deviation:', x$selparams[2], '\n\n')
  }
  else if(x$seldist == 3) {
    cat('Selection Distribution: Gamma\n')
    cat('Shape Parameter:', x$selparams[1], '\n')
    cat('Scale Parameter:', x$selparams[2], '\n\n')
  } else if(x$seldist == 4) {
    cat('Selection Distribution: Beta\n')
    cat('Shape 1 Parameter:', x$selparams[1], '\n')
    cat('Shape 2 Parameter:', x$selparams[2], '\n\n')
  } else {
    cat('Selection Distribution: ERROR\n')
    cat('Selection Input 1: ERROR\n')
    cat('Selection Input 2: ERROR\n\n')
  }
    
  if(is.null(x$fddist)) {
    cat('Frequency Dependence Distribution:', x$fddist, '\n')
    cat('Frequency Dependence Input 1:', x$fdparams[1], '\n')
    cat('Frequency Dependence Input 2:', x$fdparams[2], '\n\n')
  }
  else if(x$fddist == 1) {
    cat('Frequency Dependence Parameter Distribution: Uniform\n')
    cat('Minimum:', x$fdparams[1], '\n')
    cat('Maximum:', x$fdparams[2], '\n\n')
  }
  else if (x$fddist == 2) {
    cat('Frequency Dependence Parameter Distribution: Normal\n')
    cat('Mean:', x$fdparams[1], '\n')
    cat('Standard Deviation:', x$fdparams[2], '\n\n')
  } else {
    cat('Frequency Dependence Distribution: ERROR\n')
    cat('Frequency Dependence Input 1: ERROR\n')
    cat('Frequency Dependence Input 2: ERROR\n\n')
  }

  if(is.null(x$migdist)) {
    cat('Migration Distance Distribution:', x$migdist, '\n')
    cat('Migration Distance Input 1:', x$migdistparams[1], '\n')
    cat('Migration Distance Input 2:', x$migdistparams[2], '\n\n')
  }
  else if(x$migdist == 1) {
    cat('Migration Distance Distribution: Uniform\n')
    cat('Minimum:', x$migdistparams[1], '\n')
    cat('Maximum:', x$migdistparams[2], '\n\n')
  }
  else if(x$migdist == 2) {
    cat('Migration Distance Distribution: Normal\n')
    cat('Mean:', x$migdistparams[1], '\n')
    cat('Standard Deviation:', x$migdistparams[2], '\n\n')
  } else {
    cat('Migration Distance Distribution: ERROR\n')
    cat('Migration Distance Input 1: ERROR\n')
    cat('Migration Distance Input 2: ERROR\n\n')
  }

  if(is.null(x$migprobdist)) {
    cat('Migration Probability Distribution:', x$migprobdist, '\n')
    cat('Migration Probability Input 1:', x$migprobparams[1], '\n')
    cat('Migration Probability Input 2:', x$migprobparams[2], '\n\n')
  }
  else if(x$migprobdist == 1) {
    cat('Migration Probability Distribution: Uniform\n')
    cat('Minimum:', x$migprobparams[1], '\n')
    cat('Maximum:', x$migprobparams[2], '\n\n')
  }
  else if(x$migprobdist == 2) {
    cat('Migration Probability Distribution: Normal\n')
    cat('Mean:', x$migprobparams[1], '\n')
    cat('Standard Deviation:', x$migprobparams[2], '\n\n')
  } else {
    cat('Migration Probability Distribution: ERROR\n')
    cat('Migration Probability Input 1: ERROR\n')
    cat('Migration Probability Input 2: ERROR\n\n')
  }
}
