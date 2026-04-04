### Main Model Functions ###

#' Select community for birth
#' 
#' @param x Numeric matrix, species x site
#' 
#' @return int, index of community where birth is to occur
#' @export

samp_com_for_birth <- function(x) {
  
  if(!is.matrix(x)) {
    rlang::abort('Input must be an abundance (integer) site x species matrix',
                 class = 'input_type_error')
  }
  
  n_com <- nrow(x)
  probs <- rowSums(x)/sum(x)
  
  out <- sample(n_com, size=1, prob=probs)
  out
}

#' Calculate species probabilities
#' 
#' @param x Numeric matrix, species x site
#' @param com_index Int, row index for community to sample
#' @param params Numeric matrix, matrix of selection coefficients
#' 
#' @return int, index of species for birth
#' @export

calculate_species_probs <- function(x, com_index, params) {
  
  birth.freqs <- x[com_index,] / sum(x[com_index,])
  
  coeff <- exp(params$fd * (birth.freqs - (1/length(x[com_index,]))) + log(params$s[com_index,])) # Use the equation from Vellend 2016 to calculate density dependent selection coefficients
  probs <- (coeff*birth.freqs)/sum(coeff*birth.freqs)
  probs
}

#' Select species for birth
#' 
#' @param x Numeric matrix, species x site
#' @param probs Numeric vector, vector of species probabilities
#' 
#' @return int, index of species for birth
#' @export

samp_species_for_birth <- function(x, probs) {
  
  # Check that x is a matrix
  if(!is.matrix(x)) {
    rlang::abort('x parameter must be an abundance (integer) site x species matrix',
                 class = 'input_type_error')
  }
  
  # Check that the coefficients vector is a vector
  if(!is.vector(probs)) {
    rlang::abort('Probability vector is not a vector')
  }
  
  # Check that coefficients vector is the right length
  if(length(probs) != ncol(x)) {
    rlang::abort('Probability vector is not the correct length')
  }
  
  sample(ncol(x), size = 1, prob = probs) # Select which species reproduces
  
}

#' Select Community for death
#' 
#' @param x Numeric matrix, species x site
#' @param species_index Int, row index for community to sample
#' @param com_index Int, community of birth species
#' @param params Params object 
#' 
#' @return int, index of species for birth
#' @export

samp_com_for_death <- function(x, species_index, com_index, params) {

  sample(nrow(x), size = 1, prob = params$mig[[species_index]][,com_index])

}

#' Select species for death
#' 
#' @param x Numeric matrix, species x site
#' @param com_index Int, community index for death community
#' 
#' @return int, index of species for birth
#' @export

samp_species_for_death <- function(x, com_index) {
  
  probs = x[com_index, ]/sum(x[com_index,])
  ind.death <- sample(ncol(x), size = 1, prob = probs)
  ind.death
  
}

#' Update metacommunity matrix 
#' 
#' @param x Numeric matrix, species x site
#' @param com_index Int, community index for community where death occurs
#' @param b_species_index Int, index for species for birth
#' @param d_species_index Int, index for species for death
#' 
#' @return numeric matrix, updated matrix

update_meta <- function(x, com_index, b_species_index, d_species_index) {
  
  x[com_index, d_species_index] <- x[com_index, d_species_index] - 1
  x[com_index, b_species_index] <- x[com_index, b_species_index] + 1
  x

}

#' Conduct a single birth death process
#' 
#' @param x Numeric matrix, species x site
#' @param params params object
#' 
#' @return numeric matrix, updated matrix


birth_death_process <- function(x, params) {
  b_com <- samp_com_for_birth(x)
  s_probs <- calculate_species_probs(x, b_com, params)
  s_birth <- samp_species_for_birth(x, s_probs)
  d_com <- samp_com_for_death(x, s_birth, b_com, params)
  s_death <- samp_species_for_death(x, d_com)
  x <- update_meta(x, d_com, s_birth, s_death)
  x
}

#' Create new species
#' 
#' This function randomly creates new species in a site x species matrix if 
#' there are species (columns) with community sums equal to 0. 
#' 
#' @param x Numeric matrix, species x site matrix
#' @param params S3 class params object
#' @param prop_new float, the proportion of the community size that the new species
#' should have
#' @param diff_sd float, the standard deviation to use when evolving the selection
#' coefficients and frequency dependence parameters. These values are evolved by
#' randomly selecting a number from a normal distribution with mean 0 and standard
#' deviation per this parameter.
#' 
#' @details This function takes a site x species matrix and a params object and 
#' checks to see if there are any empty species (column) vectors. If there are
#' it randomly selects a single community (row) and adds individuals into the 
#' cell of the empty column up to the frequency specified in prop_new. To keep the 
#' community size the same, individuals are randomly subtracted from other non-zero
#' species in the community. Additionally, the selection coefficients and frequency
#' dependence parameters are modified from the extinct species parameters by 
#' adding values from a normal distribution centered at 0. 
#'
#' @return Returns a new site x species matrix with all column sums greater than
#' zero and a new params object with updated parameters for the new species
#' 
#' @export

speciate <- function(x, params, prop_new=0.1, diff_sd=0.1) {
  # Check that x is a matrix
  if(!is.matrix(x)) {
    rlang::abort('x must be an abundance (integer) site x species matrix',
                 class = 'input_type_error')
  }
  
  # Check that params is a params object
  if(!is.params(params)) {
    rlang::abort('Parameter object not of class params, check parameter input',
                 class = 'input_type_error')
  }
  
  n_coms <- nrow(x)
  n_spec <- ncol(x)
  
  if(any(colSums(x) == 0)) {
    inds <- which(colSums(x) == 0)
    
    for(i in inds) {
      # Select community for new species
      com_ind <- sample(n_coms, 1)
      com_size <- sum(x[com_ind,])
      n_to_add <- ceiling(com_size*prop_new)
      n_new <- 0
      
      while(n_new < n_to_add) {
        pos_spec <- which(x[com_ind,] > 1)
        if(i %in% pos_spec) {
          pos_spec <- pos_spec[-which(pos_spec == i)]
        }
        prob_spec <- x[com_ind,][pos_spec]/sum(x[com_ind,][pos_spec])
        if(length(pos_spec) == 1) {
          sel_spec <- pos_spec[1]
        } else {
          sel_spec <- sample(pos_spec, 1, prob=prob_spec)
        }
        x[com_ind, sel_spec] <- x[com_ind, sel_spec] - 1
        n_new <- n_new + 1
      }
      
      x[com_ind, i] <- n_new
      if(sum(x[com_ind,]) != com_size) {
        rlang::abort("New community size did not match old community size")
      }
      
      # evolve the selection coefficients
      for(j in 1:nrow(params$s)) {
        dif <- stats::rnorm(1, 0, diff_sd)
        params$s[j, i] <- params$s[j, i] + dif
        if(params$s[j, i] <= 0) {
          params$s[j, i] <- 0.00000001
        }
      }
      
      # evolve the frequency dependence parameter
      dif <- stats::rnorm(1, 0, diff_sd)
      params$fd[i] <- params$fd[i] + dif
      if(params$fd[i] > 0) {
        params$fd[i] <- 0
      }
      
    }
      
  }
  
  return(list(x, params))

}

#' Run Metacommunity Moran Simulation.
#'
#' This functions runs a single metacommunity simulation using an adapted version of
#' a population genetics Moran model with selection and migration.
#'
#' @param x Numeric matrix, starting metacommunity (species x site) matrix where columns
#' are species, rows are communities and cell ij is the count of species j in community i.
#' @param t Numeric, number of generations for simulation to run.
#' @param params S3 class params object, see \code{\link{make_params}}.
#' @param output Logical, if True outputs progress bar.
#'
#' @details Metacommunity simulations in CommSimABC are based upon Moran models with
#' selection and migration. As input, these simulations require a starting metacommunity,
#' the number of generations for the simulation to run, a selection matrix, where columns
#' are species, rows are sites, and cell ij is the selection coefficient for species j
#' in community i, frequency dependence parameters for each species, and migration matrices
#' for each species. For details on the selection matrix, frequency dependence parameters,
#' and migration matrices see the links in see also. For details on model algorithm and
#' mathematics see Williams et al. (XXXX).
#'
#' @return Returns a list containing three objects:
#' metacommunity: The simulated metacommunity
#' metafreq: A matrix where each column is a species and each row is a generation. Values
#' are the frequencies of species summed across the entire metacommunity.
#' input: a list containing the starting metacommunity and the starting params object
#'
#' @seealso
#' \code{\link{make_params}}
#' \code{\link{set_sel}}
#' \code{\link{set_fd}}
#' \code{\link{set_mig}}
#' \code{\link{make_mig}}
#'
#' @examples
#' meta <- rand_meta(N = 5, S = 5, 25)
#' xy <- matrix(sample(50, 10, replace = TRUE), ncol = 2)
#'
#' inparam <- make_params(5,5)
#' inparam$s <- set_sel(inparam, 'gamma', 2, .15)
#' inparam$fd <- set_fd(inparam, -.5, -.2)
#' inparam$mig <- set_mig(inparam, xy, 25, .05)
#'
#' moran_deme(meta, 5, inparam)
#'
#' @export

moran_deme <- function(x, t, params, output = TRUE) {
  if (requireNamespace("Rcpp", quietly = TRUE)) {

    out <- moran_deme_cpp(x=x, t=t, params=params)
    
  } else {
    
    out <- moran_deme_r(x=x, t=t, params=params, output=output)
    
  }
  
  out
  
}
  

#' @describeIn moran_deme R function that is used if Rcpp is not available
#' @export

moran_deme_r <- function(x, t, params, output = TRUE) {
  if(!is.params(params)) {
    stop('Parameter file not configured correctly')
  }
  
  params_out <- list(x, params)
  
  spp <- attr(params, 'NumSpec') # Calculate number of species

  J <- sum(x) # Calculate the total number of individuals in metacommunity
  
  out.mat <- matrix(nrow = t, ncol = spp) # Initialize the matrix of frequencies to output
  out.mat[1,] <- apply(x, 2, function(y){sum(y)/J}) # Calculate starting frequencies and add to output
  
  Gen <- 2 # Set index for generations
  if(output == T) {
    pb <- utils::txtProgressBar(min = 0, max = (t-1)*J, style = 3)
  }
  
  for(i in 1:((t-1)*J)){ # Loop through moran model until you have done t generations
    
    x <- birth_death_process(x, params)
    
    if(i %% sum(x) == 0){ # Check to see if enough iterations have occurred for a generation
      
      evol <- speciate(x, params)
      x <- evol[[1]]
      params <- evol[[2]]
      
      out.mat[Gen, ] <- apply(x, 2, function(y){sum(y)/J}) # Output current frequencies
      
      Gen <- Gen + 1 # Reset generation index
    }
    
    if(output == T) {
      utils::setTxtProgressBar(pb, i)
    }
    
  }
  
  out <- list(metacommunity = x, metafreq = out.mat, input = params_out)
  out
}


