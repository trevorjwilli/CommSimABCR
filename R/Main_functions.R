### Main Model Functions ###

#' Run Metacommunity Moran Simulation.
#'
#' This functions runs a single metacommunity simulation using an adapted version of
#' a population genetics Moran model with selection and migration.
#'
#' @param x Numeric matrix, starting metacommunity (species x site) matrix where columns
#' are species, rows are communities and cell ij is the count of species j in community i.
#' @param t Numeric, number of generations for simulation to run.
#' @param params S3 class params object, see \code{\link{make_params}}.
#' @param outgens Numeric vector of generations for which output is desired.
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
#' @return Returns an object of class simrun, list containing four objects:
#' Metacommunity: A list containing the metacommunity matrices from generations specified in
#' outgens
#' Incidence: A list containing presence-absence matrices from generations specified in
#' outgens
#' Metafreq: A matrix where each column is a species and each row is a generation. Values
#' are the frequencies of species summed across the entire metacommunity.
#' Comfreq: A list of matrices (one for each community )where each column is a species and
#' each row is a generation. Values are the frequencies of the species in that community
#' at the start of that generation.
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

moran_deme <- function(x, t, params, outgens = NULL, output = TRUE) {
  if(!is.params(params)) {
    stop('Parameter file not configured correctly')
  }
  spp <- attr(params, 'NumSpec') # Calculate number of species

  coms <- attr(params, 'NumSite') # Calculate number of communities

  comlist <- list() # make a list for output for each community (this line and next 3 lines)
  for(i in 1:coms) {
    comlist[[i]] <- matrix(nrow = t, ncol = spp)
    comlist[[i]][1,] <- x[i,]/sum(x[i,])
  }

  metas <- list() # make a list for output of metacommunity at certain time intervals
  #metas[["0"]] <- x
  incs <- list() # make a list for output of presence-absence matrices
  #incs[["0"]] <- ifelse(x>0, 1, x)

  J <- sum(x) # Calculate the total number of individuals in metacommunity

  out.mat <- matrix(nrow = t, ncol = spp) # Initialize the matrix of frequencies to output
  out.mat[1,] <- apply(x, 2, function(y){sum(y)/J}) # Calculate starting frequencies and add to output

  #print("Starting Community", quote = F)
  #print(x)

  Gen <- 2 # Set index for generations
  #print(paste("Generation", Gen), quote = F)
  if(output == T) {
    pb <- txtProgressBar(min = 0, max = (t-1)*J, style = 3)
  }

  for(i in 1:((t-1)*J)){ # Loop through moran model until you have done t generations
    #print(paste("Iteration", i), quote = F)
    #print(x)
    coms.probs <- apply(x, 1, sum) # This line and next calculate proper probabilities for calculating which community will experience a birth
    coms.probs <- coms.probs/sum(x)
    #print("Community Probabilities")
    #print(coms.probs)

    b.com <- sample(coms, size = 1, prob = coms.probs) # randomly select one community for birth process
    #print("Community Chosen")
    #print(b.com)

    com.birth <- x[b.com,] # Create vector corresponding to birth community
    #print("Birth Community")
    #print(com.birth)

    birth.freqs <- com.birth/sum(com.birth) # Calculate frequencies of each species in birth community
    #print("Species Frequencies")
    #print(birth.freqs)
    coef <- exp(params$fd * (birth.freqs - (1/length(com.birth))) + log(params$s[b.com,])) # Use the equation from Vellend 2016 to calculate density dependent selection coefficients
    #print("Selection Soefficients")
    #print(params$s[b.com,])
    #print("Frequency Dependence Parameters")
    #print(params$fd)
    #print("Computed Coefficients")
    #print(coef)
    s.probs <- (coef*birth.freqs)/sum(coef*birth.freqs) # Create vector of species birth probabilities using selection equation
    #print("Selection Probabilities")
    #print(s.probs)
    #print(paste("Sum of probs =", sum(s.probs)))
    #print(paste("Offset = ", 1 - sum(s.probs)))
    offset <- 1 - sum(s.probs) # Check to see if there is any remainder
    #print("Probabilities", quote = F)
    #print(s.probs)
    #print(paste("Offset =",offset), quote = F)
    if(offset != 0) { # If yes...
      pos <- sample(which(s.probs > 0), 1) # See which species have positive probabilities and pick one
      #print(paste("Species chosen:", pos))
      s.probs[pos] <- s.probs[pos] + offset # Add remainder to the chosen species ensures that probability weights sum to 1
      offsetcheck <- 1 - sum(s.probs)
      #print(paste("New offset:", offsetcheck))
    }

    species.birth <- sample(length(com.birth), size = 1, prob = s.probs) # Select which species reproduces
    #print("Species chosen for birth")
    #print(species.birth)
    mig.mat <- params$mig[[species.birth]] # Index the species migration matrix
    #print("Migration Matrix")
    #print(mig.mat)
    d.com <- sample(coms, size = 1, prob = mig.mat[,b.com]) # Sample to see if offspring migrates or stays in same community
    #print("Community chosen for Death")
    #print(d.com)
    com.death <- x[d.com,] #Select community where death occurs
    #print("Community for Death")
    #print(com.death)
    ind.death <- sample(length(com.death), size = 1, prob = com.death/sum(com.death))
    #print("Species chosen to die")
    #print(ind.death)
    com.death[ind.death] <- com.death[ind.death] - 1
    com.death[species.birth] <- com.death[species.birth] + 1

    x[d.com,] <- com.death
    #print(x)
    if(i %% sum(x) == 0){ # Check to see if enough iterations have occurred for a generation

      out.mat[Gen, ] <- apply(x, 2, function(y){sum(y)/J}) # Output current frequencies

      for(j in 1:coms){
        comlist[[j]][Gen,] <- x[j,]/sum(x[j,])
      }

      if(!is.null(outgens)) { # See if you want multiple metacommunity matrices in output
        if(Gen %in% outgens) { # Test to see if the current generation is to be put in output
          metas[[as.character(Gen)]] <- x # output matrix
          incs[[as.character(Gen)]] <- ifelse(x>0, 1, x)
        }
      }

      Gen <- Gen + 1 # Reset generation index
      #print("Generation")
      #print(Gen)
    }

    if(output == T) {
      setTxtProgressBar(pb, i)
    }

  }

  inc <- ifelse(x>0, 1, x)
  incs[[as.character(t)]] <- inc
  metas[[as.character(t)]] <- x
  out <- list(Metacommunity = metas, Incidence = incs, metafreq = out.mat,comfreq = comlist)
  class(out) <- c('simrun', 'list')
  out
}


#' @export

plot.simrun <- function(x, lgnd = T) {

  clrs <- ggsci::pal_igv()(51)

  # Plot whole metacommunity
  plot(x$metafreq[,1], type = "l", xlab = "Generations", ylab = "Frequency", col = clrs[1],
       ylim = c(0,1), main = "Meta-Community") # Plot frequency graph of first species
  for(j in 2:length(x$metafreq[1,])){ # Plot the rest of the spieces
    lines(x$metafreq[,j], col = clrs[j])
  }
  txt <- paste("Species", 1:length(x$metafreq[1,]))
  if(lgnd == T) {
    legend("topleft", legend = txt, col = clrs, lty = 1) # Add legend
  }

  for(i in 1:length(x$comfreq)){
    plot(x$comfreq[[i]][,1], type = "l", xlab = "Generations", ylab = "Frequency", col = clrs[1],
         ylim = c(0,1), main = paste("Community", i))
    for(j in 2:length(x$comfreq[[1]][1,])) {
      lines(x$comfreq[[i]][,j], col = clrs[j])
    }
    txt <- paste("Species", 1:length(x$metafreq[1,]))
    if(lgnd == T) {
      legend("topleft", legend = txt, col = clrs, lty = 1) # Add legend
    }
  }
}


