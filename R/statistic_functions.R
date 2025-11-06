### Statistic Functions ###

#' Calculate Community Evenness.
#'
#' Calculates Pielou's evenness of a metacommunity.
#'
#' @param x A metacommunity matrix (columns = species, rows = communities, cell values = counts).
#'
#' @return A vector of evenness values for each community.
#'
#' @examples
#' meta <- rand_meta(20, 30, 500)
#' evenness(meta)
#'
#' @export

evenness <- function(x) {
  H <- vegan::diversity(x)
  J <- H/log(vegan::specnumber(x))
  if(NaN %in% J) {
    J[is.nan(J)] <- 0
  }
  return(J)
}


#' Estimate beta-diversity using sum of squared interpoint dissimilarities
#'
#' Calculates the sum of squared interpoint dissimilarities of a metacommunity
#' matrix using the formula in Anderson et al. (2011) in Ecology Letters
#'
#' @param x metacommunity matrix
#' @param method type of dissimilarity metric used; one of the methods in vegdist()
#'
#' @details Estimates beta-diversity of a metacommunty matrix by using the sum of
#' squared interpoint dissimilarities: 1/(N(N-1))*sum(dij^2) where dij is the distance
#' between community i and j.
#'
#' @export

ssid2beta <- function(x, method = "bray") {
  N <- nrow(x)
  dist.x <- vegan::vegdist(x, method = method)
  sigma <- (1/(N*(N-1)))*sum(dist.x^2)
  sigma
}


#' Count Embedded Absences
#'
#' This function counts the number of embedded absences in a binary matrix
#' to calculate the coherence of a metacommunity
#'
#' @param x A metacommunity matrix with species as columns and sites as rows
#'
#' @details This function counts the number of embedded absences in a sitexspecies
#' matrix following Leibold and Mikkelson 2002.
#'
#' @examples
#' meta <- rand_meta(10, 15, 100)
#' embabs(meta)
#'
#' @export

embabs <- function(x) {
  countbetween <- function(zeros, ones) { # Internal function to count the number
    # of zeros between ones in a vector
    
    valrange <- function(ones) { # Internal function to make a list of ranges between
      # all 1s in a vector
      if(length(ones) > 1) {
        t <- 1
        out <- list()
        while(t < length(ones)) {
          vec <- c(ones[t], ones[t+1])
          #print(vec)
          out[[paste(t)]] <- vec
          t <- t + 1
        }
      } else {
        out <- list(c(0,0))
      }
      out
    }
    
    ranges <- valrange(ones)
    #print(ranges)
    
    multibetween <- function(x, ranges) { # Internal function to run dplyrs
      #between on a list of ranges
      out <- vector()
      for(i in 1:length(ranges)) {
        test <- dplyr::between(x, ranges[[i]][1], ranges[[i]][2])
        out <- append(out, test)
      }
      out
    }
    
    numbet <- sapply(zeros, multibetween, ranges = ranges)
    #print(numbet)
    if(is.list(numbet)) {
      numbet <- 0
      return(list(count = numbet, indices = NULL))
    } else if(is.matrix(numbet)) {
      ind <- zeros[which(colSums(numbet) > 0)]
      numbet <- sum(numbet)
      return(list(count = numbet, indices = ind))
    } else {
      numbet <- sum(numbet)
      return(list(count = numbet, indices = NULL))
    }
  }
  
  x <- metacom::OrderMatrix(x)
  #print(colSums(x))
  inds <- vector()
  count <- 0
  for(i in 1:nrow(x)) {
    #print(i)
    zeros <- which(x[i,] == 0)
    ones <- which(x[i,] == 1)
    num <- countbetween(zeros, ones)
    count <- count + num$count
    #print(count)
    if(num$count > 0) {
      for(j in 1:length(num$indices)) {
        index <- paste(i, num$indices[j])
        inds <- append(inds, index)
      }
    }
    
  }
  #cat("Rows Done\n")
  for(k in 1:ncol(x)) {
    zeros <- which(x[,k] == 0)
    ones <- which(x[,k] == 1)
    num <- countbetween(zeros, ones)
    #print(num)
    if(length(num$indices) > 0) {
      for(q in 1:length(num$indices)) {
        index <- paste(num$indices[q], k)
        #print(index)
        if(index %in% inds) {
          next
        } else {
          count <- count + 1
          inds <- append(inds, index)
        }
      }
    }
  }
  return(list(count = count, inds = inds, orderdmat = x))
}


#' Calculate Species Turnover
#'
#' Calculate the species turnover of a metacommunity matrix
#'
#' @param x A metacommunity matrix with species as columns and sites as rows
#'
#' @details This function calculates species turnover following (TODO: add reference)
#'
#' @examples
#' meta <- rand_meta(10, 15, 100)
#' turnover(meta)
#'
#' @export

turnover <- function(x) {
  x <- metacom::OrderMatrix(x, scores = 1, binary = TRUE)
  for (i in 1:ncol(x)) {
    x[min(which(x[, i] == 1)):max(which(x[, i] == 1)), i] <- 1
  }
  D <- vegan::designdist(x, method = "(A-J)*(B-J)", terms = "minimum")
  return(sum(D))
}


#' Calculate Morisita's Index
#'
#' Calculate Morisita's index for a metacommunity
#'
#' @param x A metacommunity matrix with species as columns and sites as rows
#'
#' @details This function calculates Morisita's index following (TODO: add reference)
#'
#' @examples
#' meta <- rand_meta(10, 15, 100)
#' morisitas(meta)
#'
#' @export

morisitas <- function (x) {
  x <- metacom::OrderMatrix(x, scores = 1)
  for (i in 1:ncol(x)) {
    x[min(which(x[, i] == 1)):max(which(x[,i] == 1)), i] <- 1
  }
  x <- t(x)
  M <- 0
  ComBnd <- rep(0, ncol(x))
  ComBndChi <- 0
  for (i in 1:nrow(x)) {
    ind1 <- which(x[i, ] == 1)
    for (j in 1:ncol(x)) {
      if (min(ind1) == j) {
        ComBnd[j] = ComBnd[j] + 1
      }
      if (max(ind1) == j) {
        ComBnd[j] = ComBnd[j] + 1
      }
    }
  }
  TotComBnd <- (nrow(x) * 2) - ComBnd[1] - ComBnd[ncol(x)]
  ExpComBnd <- TotComBnd/(ncol(x) - 2)
  df <- -1
  for (z in 2:(ncol(x) - 1)) {
    M <- M + ((ComBnd[z]/TotComBnd) * ((ComBnd[z] - 1)/(TotComBnd -
                                                          1)))
    ComBndChi <- ComBndChi + (((ComBnd[z] - ExpComBnd)^2)/ExpComBnd)
    df <- df + 1
  }
  M <- M * (ncol(x) - 2)
  M
}


#' Calculate True Alpha Diversity
#'
#' Function to calculate the true alpha diversity of a metacommunity
#'
#' @param x Numeric matrix. A meta-community matrix (species x site)  where columns
#' are species, rows are communities and cell ij is the count of species j in community i.
#' @param q Numeric. Order of the diversity measure.
#' @param weight Logical indicating if communities should be given different weights.
#'  Following Jost 2007, this is only valid when q=1 so if set to `TRUE`, q is 
#'  set to 1 and weights are set according to the sizes of individual communities.
#'
#' @details This function calculates the true alpha diversity using numbers equivalents
#' as detailed in equations 11a and 11b of Jost 2007 (see also, Hill 1973 and Jost 2006).
#'  The equation follows as:
#' \ifelse{html}{\out{{}^{q}D_{\alpha} = [{\frac{1}_{N}} \sum_{i=1}^{N} \sum_{j=1}^{S} p^{q}_{ij}]^{1/1-q}}}{\eqn{{}^{q}D_{\alpha} = [{\frac{1}_{N}} \sum_{i=1}^{N} \sum_{j=1}^{S} p^{q}_{ij}]^{1/1-q}}}
#' for \eqn{q  1} and #' \ifelse{html}{\out{{}^{q}D_{\alpha} = e^{{-\frac{1}_{N}} \sum_{i=1}^{N} \sum_{j=1}^{S} [p_{ij}*ln(p_{ij})]}}}{\eqn{{}^{q}D_{\alpha} = e^{{-\frac{1}_{N}} \sum_{i=1}^{N} \sum_{j=1}^{S} [p_{ij}*ln(p_{ij})]}}}
#' Where q is the order of the diversity measure, N is the number of communities,
#' S is the number of species, and p_ij is the frequency of species j in community i
#'
#' @return Numeric, the true alpha diversity
#'
#' @examples
#' x <- rand_meta(10, 6, 50)
#' d_alpha(x, q=2)
#'
#' @export

d_alpha <- function(x, q=1, weight=FALSE) {

  # If there is a single community turn into matrix
  if(is.vector(x)) {
    x <- t(as.matrix(as.integer(x)))
  }
  
  # Ensure input is a site x species abundance matrix
  if(!is.matrix(x)) {
    rlang::abort('Input must be an abundance site x species matrix',
                 class = 'input_type_error')
  }
  
  if(!(typeof(x) == 'integer')) {
    rlang::abort('Input must be an abundance site x species matrix with type integer',
                 class = 'input_type_error')
  }
  
  # Ensure q >= 0
  if(q < 0) {
    rlang::abort('q must be greater than or equal to 0', class = 'q_domain_error')
  }
  
  # Ensure weight is boolean
  if(!is.logical(weight)) {
    rlang::abort('Weight must be boolean, FALSE if unweighted, TRUE if weighted',
                 class = 'weight_type_error')
  }
  
  # Convert abundance matrix to frequency matrix
  x_norm <- normalize_meta(x)
  
  # If weights are to be used set q to 1
  # and w to proportion of total meta-community size
  if(weight) {
    rlang::warn('Weights only valid for q=1, setting q to 1')
    q <- 1
    w <- rowSums(x)/sum(x)
  } 
  else {
    # Set default weight
    w <- rep(1/nrow(x_norm), nrow(x_norm))
  }
  
  # Jost 2007 eqn. 11a
  if(q != 1) {
    if(q == 0) {
      x_norm <- apply(x_norm, 2, \(j) ifelse(j>0, 1, 0))
      if(is.vector(x_norm)) {
        x_norm <- t(as.matrix(x_norm))
      }
      numerator <- sum((w^q)*(rowSums(x_norm)))
    }
    else {
      numerator <- sum((w^q)*(rowSums(x_norm^q)))
    }
    denominator <- sum(w^q)
    metric <- (numerator/denominator)^(1/(1-q))
  }
  
  # Jost 2008 eqn. 11b
  else if(q == 1) {
    sums <- x_norm*log(x_norm)
    sums <- rowSums(sums, na.rm = TRUE)*-w
    metric <- exp(sum(sums))
  }
  
  metric
  
}

#' Calculate True Gamma Diversity
#'
#' Function to calculate the true gamma diversity of a metacommunity
#'
#' @param x Numeric matrix. A meta-community abundance matrix (species x site)  where columns
#' are species, rows are communities and cell ij is the count of species j in community i.
#' @param q Numeric. Order of the diversity measure.
#' @param weight Logical indicating if communities should be given different weights.
#'  Following Jost 2007, this is only valid when q=1 so if set to `TRUE`, q is 
#'  set to 1 and weights are set according to the sizes of individual communities.
#'
#' @details This function calculates true gamma diversity using numbers equivalents
#' as detailed in Jost 2008. 
#'
#' @return Numeric, the true alpha diversity
#'
#' @examples
#' x <- rand_meta(10, 6, 50)
#' d_gamma(x, q=1)
#'
#' @export

d_gamma <- function(x, q=1, weight=FALSE) {
  
  # Ensure input is a site x species abundance matrix
  if(!is.matrix(x)) {
    rlang::abort('Input must be an abundance (integer) site x species matrix',
                 class = 'input_type_error')
  }
  
  if(!(typeof(x) == 'integer')) {
    rlang::abort('Input must be an abundance (integer) site x species matrix',
                 class = 'input_type_error')
  }
  
  # Ensure q >= 0
  if(q < 0) {
    rlang::abort('q must be greater than or equal to 0', class = 'q_domain_error')
  }
  
  # If weights are to be used set q to 1
  # and w to proportion of total meta-community size
  if(weight) {
    rlang::warn('Weights only valid for q=1, setting q to 1')
    q <- 1
    w <- rowSums(x)/sum(x)
  } 
  else {
    # Set default weight
    w <- rep(1/nrow(x), nrow(x))
  }
  x <- normalize_meta(x)
  if(q == 1) { # Jost eqn. 17b
    p1 <- colSums(w*x)*-1
    p2 <- log(colSums(w*x))
    metric <- exp(sum(p1*p2))
  } else { # Jost 2007 Errata
    nsites <- nrow(x)
    tmp <- (colSums(x)/nsites)^q
    metric <- sum(tmp)^(1/(1-q))
  }
  
  metric
}

#' Calculate True Beta Diversity
#'
#' Function to calculate the true beta diversity of a metacommunity
#'
#' @param x Numeric matrix. A meta-community abundance matrix (species x site)  where columns
#' are species, rows are communities and cell ij is the count of species j in community i.
#' @param q Numeric. Order of the diversity measure.
#' @param weight Logical indicating if communities should be given different weights.
#'  Following Jost 2007, this is only valid when q=1 so if set to `TRUE`, q is 
#'  set to 1 and weights are set according to the sizes of individual communities.
#'
#' @details This function calculates true beta diversity using numbers equivalents
#' as detailed in Jost 2007 equation 17c. 
#'
#' @return Numeric, the true alpha diversity
#'
#' @examples
#' x <- rand_meta(10, 6, 50)
#' d_beta(x, q=1)
#'
#' @export

d_beta <- function(x, q=1, weight=FALSE) {
  alpha <- d_alpha(x, q, weight)
  gamma <- d_gamma(x, q, weight)
  gamma/alpha
}



