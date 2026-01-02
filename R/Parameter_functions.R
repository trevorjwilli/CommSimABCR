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
#' @param distr Character, one of 'uniform', 'normal', 'gamma', or 'beta'.
#' @param input1,input2 Numeric, parameters for distribution being used. See details.
#'
#' @details This function generates a random selection matrix according to a specified
#'   distribution. The specified distribution can be one of 'uniform', 'normal', or 'gamma.'
#'   If 'uniform' input1 is the minimum value and input2 is the maximum value. If 'normal'
#'   input1 is the mean and input2 is the standard deviation. If 'gamma' input1 is the
#'   shape and input2 is the scale parameters. if 'beta' input1 is the first shape
#'   parameter and input2 is the second shape parameter
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
  
  if(!(distr %in% c('uniform', 'normal', 'gamma', 'beta'))) {
    stop("distr must be one of: 'uniform', 'normal', 'gamma', 'beta'")
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
      coeffs <- stats::runif(n.spec, min = input1, max = input2)
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
      coeffs <- stats::rnorm(n.spec, mean = input1, sd = input2)
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
      coeffs <- stats::rgamma(n.spec, shape = input1, scale = input2)
      coeffs[which(coeffs < 0)] <- 0.00000001 # Make sure there are no 0s
      sel[j,] <- coeffs
    }
    ##print(sel)
  } else if(distr == 'beta') {
    sel <- matrix(ncol = n.spec, nrow = n.sites)
    if(input1 < 0 | input2 < 0) {
      stop('Shape parameters must be positive')
    }
    for(j in 1:n.sites) {
      coeffs <- stats::rbeta(n.spec, shape1 = input1, shape2 = input2)
      sel[j,] <- coeffs
    }
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
  fd <- stats::runif(n.spec, input1, input2)
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

  if(any(class(site.arrange) == "dist")) {
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

#' Create a Random Spatial Landscape.
#'
#' This function creates a spatially explicit environmental landscape and
#' samples random points from it for simulations of Moran Metacommunity
#' model. Environmental landscapes can be generated randomly or with spatial
#' autocorrelation using IDW interpolation.
#'
#' @param n.sites Numeric 1L, Number of communities in metacommunity.
#' @param n.meas Numeric 1L, Number of environmental measures to simulate.
#' @param x.max,y.max Numeric 1L, Maximum distance in x and y directions
#' for landscape extent.
#' @param env.maxmin List with length equal to n.meas. Each element in the list
#' is a Numeric vector with two elements, the minimum and maximum values for
#' each environmental variable.
#' @param idp Numeric, power parameter for IDW interpolation.
#' @param autocor Logical, if TRUE uses IDW interpolation to create landscape.
#' @param steps # Numeric, number to divide range by to create gradient of interpolation
#' point values (see details).
#' @param reps Numeric, Number of repetitions for each interpolation point value.
#'
#' @details This function generates a spatially explicit landscape that can be used
#' in Moran metacommunity simulations. Environmental variables are mapped to the
#' landscape in one of two ways. If autocor = TRUE, then the function generates
#' a vector of values (using the steps and reps arguments) between the maximum
#' and minimum of of each environmental variable being simulated. For example,
#' say that the user inputs that pH in the landscape ranges from 6 to 8. If
#' steps = 4 and rep = 2, the function creates a vector with the values
#' (6, 6, 6.5, 6.5, 7, 7, 7.5, 7.5, 8, 8). These values are then shuffled into a
#' random order and assigned to randomly placed points on the landscape.
#' From these points, the function uses IDW interpolation to fill in the environmental
#' variables for the rest of the landscapes. If autocor = FALSE, then environmental
#' values are randomly generated for each cell using a uniform distribution.
#' NOTE: the landscape generated is in raster format.
#'
#' @return Returns a list with the following elements:
#'
#' xy: A dataframe with the x and y coordinates of communities to be sampled
#'
#' env.vals: A dataframe with the environmental values sampled from each of the
#' communities.
#'
#' rasters: A raster brick with all of the environmental landscapes
#'
#' sp.points: A SpatialPointsDataFrame with the x and y coordinates of the
#' communities to be sampled.
#'
#' @examples
#'
#' pH <- c(6, 8)
#' eC <- c(0.01, 0.15)
#' temp <- c(70, 85)
#'
#' envlist <- list(pH, eC, temp)
#'
#' make_spatialenv(n.sites = 10, n.meas = 3, env.maxmin = envlist)
#' make_spatialenv(n.sites = 10, n.meas = 3, env.maxmin = envlist, autocor = FALSE)
#'
#' @export


make_spatialenv <- function(n.sites, n.meas, env.maxmin,
                            x.max = 100, y.max = 100, idp = 2,
                            autocor = TRUE, steps = 4, reps = 2) {

  ext <- raster::extent(0, x.max, 0, y.max) # set extent for sampling area

  rast.list <- list() # Create empty list to store rasters for each variable

  ## The following code chunk is to create a vector of values (in random order)
  ## between the min and max of each environmental variable to use as starting
  ## values for the IDW interpolation

  val.points <- list()
  for(j in 1:length(env.maxmin)) {
    val.points[[j]] <- sample(rep(seq(from = env.maxmin[[j]][1], to = env.maxmin[[j]][2],
                                      by = (env.maxmin[[j]][2]-env.maxmin[[j]][1])/steps), each = reps))
  }

  for(i in 1:n.meas) { # Loop to create landscape for each variable

    if(autocor == TRUE) { # If you set the option to give spatial autocorrelation use the following code

      ## This method creates a blank raster and fills the raster by randomly selecting
      ## 'n.points' number of points on the raster and giving those points a value
      ## according to the values given in 'val.points' (which can be randomized)
      ## It then fills in the raster through an inverse distance weighted (IDW)
      ## method of interpolation

      n.points <- (steps+1) * reps
      print(n.points)
      rast <- raster::raster(ext = ext) # Create initial raster
      x.all <- raster::xFromCol(rast) # Give all x values to be sampled
      y.all <- raster::yFromRow(rast) # Give all y values to be sampled

      x <- sample(x.all, n.points, replace = T) # Sample x values for points used for interpolation from raster
      y <- sample(y.all, n.points, replace = T) # Sample y values for points used for interpolation from raster


      val <- data.frame(v = val.points[[i]]) # Assign values for the points used for interpolation
      #print(val)
      p <- data.frame(x, y) # Make x and y values a data frame
      #print(p)
      p.sp <- sp::SpatialPointsDataFrame(p, val) # make a spatial points data frame from the x and y and values given for the points to be interpolated from

      gs <- gstat::gstat(formula=v~1, data = p.sp, set=list(idp=idp)) # Create a gstat object for the interpolation

      IDW <- raster::interpolate(rast, gs) # perform interpolation using the IDW method

      rast.list[[i]] <- IDW # Output raster for measure
    } else if(autocor == FALSE) {

      rast <- raster::raster(ext = ext, resolution = 1) # Create initial raster
      raster::values(rast) <- stats::runif(raster::ncell(rast), env.maxmin[[i]][1], env.maxmin[[i]][2]) # Randomly assign raster values according to highs and lows
      rast.list[[i]] <- rast # Output raster for measure

    }

  }

  rast.all <- raster::brick(rast.list) # make a single raster with all the environmental measures

  points <- data.frame(raster::sampleRandom(rast.all, size = n.sites, xy = T)) # Sample random points from raster
  #print(points)
  sp.points <- sp::SpatialPointsDataFrame(points[,1:2], data = points[,-c(1,2)]) # Create a spatial points data frame for output and plotting
  #print(sp.points)
  env.vals <- points[,-c(1,2)] # Extract enviornmental values for each site

  return(list(xy = points[,1:2], env.vals = env.vals, rasters = rast.all, sp.points = sp.points))
}

#'
#' @export

print.params <- function(x, ...) {
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

#' Create Selection Matrices from Environmnetal data
#'
#' This function creates random selection matrices by generating
#' resource-utlization niches for N number of species according
#' to inputted environmental data for communities.
#'
#' @param n.spec Integer, number of species in metacommunity
#' @param env Matrix or Data-Frame with environmental values measured
#' from communities.
#' @param env.maxmin List with length equal to number of environmental measures.
#' Each element in the list is a Numeric vector with two elements, the minimum
#' and maximum values for each environmental variable.
#' @param var.maxmin List with length equal to number of environmental measures.
#' Each element in the list is a Numeric vector with two elemens, the minimum
#' and maximum values that the variance can be for each environmental variable.
#' @param meanvals Matrix, user inputed value of means for each species. This is
#' an m x n matrix where columns are environmental variables and rows are species.
#' Each cell is the value of the environmental variable at which the species is
#' found at the highest frequency.
#' @param covar User inputed list of covariance matrices.
#'
#' @details This function creates a selection matrix by randomly generating
#' multivariate normal Resource Utilization niches and matching those to generated
#' environmental variables for each community using \code{\link{make_spatialenv}}.
#' First, values for each environmental variable are calculated  for each
#' species using a random uniform distribution. This value is considered the peak
#' for the species niche for that environmental variable  (i.e. the value at which
#' the species is found at the highest frequency). Each species is then generated
#' a random variance for each environmental variable according to the var.maxmin
#' argument. At this point, this function assumes that the covariance for each
#' environmental variable is 0. Alternatively, the user can input means and covariances
#' by hand. From these values, the multivariate density for each species is calculated.
#' Species coefficients are then generated by dividing the density of each species in each
#' community by the density of the species with the max density.
#'
#' @return A matrix of selection coefficients where columns are species, rows are
#' communities, and cell ij is the selection coefficient for species j in community i.
#'
#' @examples
#'
#' testmaxmin <- list(pH = c(6, 8), eC = c(.01, .15), temp = c(70, 85))
#' testenv <- make_spatialenv(n.sites = 10, n.meas = 3, env.maxmin = testmaxmin, autocor = TRUE)
#' testvar <- list(pH = c(2,5), eC = c(.05, .1), Temp = c(5, 8))
#' make_selfromenv(n.spec = 5, env = testenv$env.vals, env.maxmin = testmaxmin, var.maxmin = testvar)
#'
#' @export

make_selfromenv <- function(n.spec, env, env.maxmin, var.maxmin, meanvals = NULL, covar = NULL) {

  n.sites <- length(env[,1])
  n.meas <- length(env[1,])

  vec.env <- unlist(env.maxmin)
  val.lows <- vec.env[c(TRUE, FALSE)]
  val.highs <- vec.env[c(FALSE, TRUE)]
  #print(val.lows)
  #print(val.highs)

  if(is.null(meanvals)) {

    spec.means <- mapply(stats::runif, n = n.spec, min = val.lows, max = val.highs) # randomly sample means (i.e. environmental values species i is found at highest frequency)
    # between max and min values of environmental variables for each species

  } else if(!is.null(meanvals)) {
    spec.means <- meanvals
  }

  #print("Species means")
  #print(spec.means)

  if(is.null(covar)) {
    sig.list <- list() # Create list (1 element for each species) for covariance matrix for the environmental variables
    for(i in 1:n.spec) { # Randomly generate covariance matrices for each spieces
      vars <- vector(length = n.meas)
      for(j in 1:n.meas) {
        vars[j] <- stats::runif(1, min = var.maxmin[[j]][1], var.maxmin[[j]][2])
      }
      sig <- diag(vars)
      sig.list[[i]] <- sig
    }
  } else if(!is.null(covar)) {
    sig.list <- covar
  }
  #print(sig.list)

  dens.mat <- matrix(nrow = n.sites, ncol = n.spec) # Calculate densities for each species in each site (i.e. which species is most likely to succeed in each site)
  for(i in 1:n.sites) {
    for(j in 1:n.spec) {
      dens.mat[i,j] <- mvtnorm::dmvnorm(env[i,], spec.means[j,], sig.list[[j]]) # Calculates density using each species multivariate environmental niche
    }
  }
  #print(dens.mat)

  sel.mat <- dens.mat # convert density matrix into selection matrix by dividing each species density within a site by the density of the species with the highest density for that site
  for(i in 1:length(sel.mat[,1])) {
    sel.mat[i,] <- sel.mat[i,]/max(sel.mat[i,])
  }

  sel.mat[sel.mat < 0.01] <- 0 # Clean up selection matrix
  #print("Selection Matrix")
  #print(sel.mat)
  return(sel.mat)

}


#' Check if object is Params class
#'
#' Checks to see if an object is of class params
#'
#' @param x object to test
#'
#' @return Logical
#'
#' @examples
#'
#' test <- make_params(5,5)
#' is.params(test)
#'
#' @export

is.params <- function(x) {inherits(x, "params")}




