### Utility Functions ###

#' Create Migration Matrix.
#'
#' Function to create a migration matrix from x,y coordinates or pairwise distance matrix.
#'
#' @param site.arrange Matrix or Distance object. Either a two column matrix with x-values
#'   \(i.e. longitude\) in the first column and y-values (i.e. latitude) in the second
#'   column or a pairwise distance object between each site.
#' @param max.dist Numeric, the maximum distance the species can migrate.
#' @param tot Numeric, a value between 0 and 1 detailing the probability of migration.
#'
#' @details This function creates a matrix of migration probabilities where cell ij
#'   is the probability for a species to migrate from community j to community i. Migration
#'   probablities are calculated by dividing the total probability of migrating into
#'   proportions between all the communities that the species can reach according to the
#'   max.dist parameter. Proportions are set up so that they are inversely related to the
#'   distance travelled, i.e. the further the distance between sites the lower the probability
#'   of migrating between sites. If x,y coordinates are given, pairwise euclidean distances
#'   are calcuted and used.
#'
#' @return Returns a mumeric matrix where cell ij is the probability of the species migrating
#' from community j to community i.
#'
#' @examples
#' xy <- matrix(sample(50, 10), ncol = 2)
#' make_mig(xy, 20, .1)
#'
#' pairdist <- stats::dist(xy)
#' make_mig(pairdist, 20, .1)
#'
#' @export

make_mig <- function(site.arrange, max.dist, tot) {
  if(any(class(site.arrange) == "dist")) {
    n.sites <- attr(site.arrange, "Size")
    mig.dist <- site.arrange
  } else {
    n.sites <- (length(site.arrange[,1]))
    mig.dist <- stats::dist(site.arrange, upper = T, diag = T) # Calculates the euclidean distances between sites
  }

  mig <- mig.dist # copy mig.dist so as to not overwrite it in next line
  mig[mig > max.dist] <- 0 # Make all values greater than max.dist = 0
  mig <- as.matrix(mig) # Turn into a matrix
  #print(mig)

  inv.propor <- function(x, tot) { # Function to calculate the probabilities of migrating to communities with probabilities inversely proportional to distance from birth community
    y <- tot/(x*sum(1/x[which(x > 0)])) # This is the solution to the system of equations where c = xy_i where c is a constant
    # x is the distance from the birth community and y_i is the probability of migrating to community i
    # In this system c = tot/D_i * sum(1/D_i), where tot is the total probability of migrating from birth community
    # and D_i is the distance between the birth community and the community i to which the species can migrate to
    y <- ifelse(y == Inf, 0, y) # need this to convert the infinities created by dividing the total migration probability by 0 (i.e. to communities where migration can't occur)
    return(y)
  }

  mig.fin <- apply(mig, 2, inv.propor, tot = tot) # Apply inv.propor to all communities
  #print(colSums(mig.fin))
  #rownames(mig.fin) <- NULL
  #colnames(mig.fin) <- NULL

  diag(mig.fin) <- 1 - tot # Give probability of staying in local community

  for(i in 1:n.sites) {
    n.zero <- sum(mig.fin[,i] == 0) # See how many zeros are in the matrix for that column
    if(n.zero == n.sites - 1) { # If there is only positive numbers on diagonal...
      mig.fin[i,i] <- 1 # Change diagonal to 1 (i.e. can't migrate)
    }
    # Following code checks to see if probability vectors actually sum to 1 and if not fixes it
    checksum <- sum(mig.fin[,i]) # Sum column
    offset <- 1 - checksum # Calculate difference between column and 1
    nums <- which(mig.fin[,i] < 1-tot & mig.fin[,i] > 0) # See if column permits migration
    if(length(nums > 0)) { # If so...
      update <- sample(nums, 1) # Randomly select one of the migration probabilities
      mig.fin[update, i] <- mig.fin[update, i] + offset # and add the offset to it
    }
  }

  return(mig.fin)
}

#' Random metacommunity matrix.
#'
#' Creates a random metacommunity matrix where cell ij is the number of species j in
#' community i.
#'
#' @param N Numeric, the number of communities.
#' @param S Numeric, the number of species.
#' @param J Numeric vector, the number of individuals within each community. If
#' length = 1 gives all communities the same number of individuals.
#' @param min.spec Numeric, the minimum number of species within communities.
#'
#' @details This function creates a random metacommunity matrix in which species
#' absolute abundances \(counts\) within communities are recorded.
#'
#' @return Returns a numeric matrix where rows are communities, columns are species,
#' and cell ij is the count of species j in community i.
#'
#' @examples
#' rand_meta(N = 5, S = 5, J = c(100, 100, 200, 40, 300), min.spec = 2)
#' rand_meta(N = 5, S = 5, J = 500)
#'
#' \dontrun{
#' rand_meta(N = 5, S = 5, J = c(10, 10))
#' }
#'
#' @export

rand_meta <- function(N, S, J, min.spec = 2) {
  if(length(J) == 1) {
    J <- rep(J, N)
  }
  else if(length(J) != 1 & length(J) != N) {
    stop("J is not the proper length")
  }
  meta <- matrix(nrow = N, ncol = S)
  for(i in 1:N) {
    spec <- sample(min.spec:S, 1) # Choose how many species to include
    w.spec <- sample(S, spec) # Choose which species to include
    probs <- tabulate(w.spec, S)/length(w.spec) # Calculate probability vector for species chosen
    com <- sample(S, size = J[i], replace = T, prob = probs)
    meta[i,] <- tabulate(com, S)
  }
  return(meta)
}


#' Random Cartesian Points
#'
#' Creates a two-column matrix where the first column is x-values and the second column
#' is y-values and each row represents a point. Can be used to randomly generate a
#' spatially explicit list of sites.
#'
#' @param n.sites Numeric, the number of sites to generate
#' @param x.max,y.max Numeric, the maximum distance in x or y directions that sites can
#' be placed.
#'
#' @return Returns a two-column matrix of coordinate points.
#'
#' @examples
#' random_points(20, 100, 100)
#'
#' @export

random_points <- function(n.sites, x.max, y.max) {
  xpoints <- sample(x.max, n.sites)
  ypoints <- sample(y.max, n.sites)
  points <- cbind(xpoints, ypoints)
  return(points)
}

#' Calculate Community Evenness.
#'
#' Calculates Pielou's evenness of a metacommunity.
#'
#' @param y A metacommunity matrix (columns = species, rows = communities, cell values = counts).
#'
#' @return A vector of evenness values for each community.
#'
#' @examples
#' meta <- rand_meta(20, 30, 500)
#' evenness(meta)
#'
#' @export

evenness <- function(y) {
  H <- vegan::diversity(y)
  J <- H/log(vegan::specnumber(y))
  if(NaN %in% J) {
    J[is.nan(J)] <- 0
  }
  return(J)
}

#' Calculate Summary Statistics
#'
#' Calculates summary statistics on metacommunities to be used in ABC analyses.
#'
#' @param y Either a metacommunity matrix or a list of metacommunity matrices
#'
#' @details Calculates the mean, median, and standard deviations for shannon diversity index,
#' simpson diversity index, Pielou's evenness, richness, the number of sites occupied for species.
#' It also calculates beta diversity by calculating the mean and sum of squared
#' interpoint distances using pairwise Bray-Curtis dissimilarities (Anderson 2011),
#' and true alpha, beta, and gamma values using Hills numbers for q = 0 through 2 Jost (2006, 2007),
#' These values are calculated using equal weights for all values of q as well as weighted
#' for q = 1. It also calculates the C-score (as well as skew and variance of C-score), checkerboard score,
#' V-ratio. Lastly it calculates the number of embedded absences, Turnover, and Morisitas
#'  index following Leibold and Mikkelson (2002).
#'  
#'  ecosumstats2() does the same thing but is written for only a single matrix, which
#'  allows parallelization through the parallel package.
#'
#' @return Either a vector of summary statistics for a single metacommunity or a dataframe where columns
#' are different summary statistics and rows are metacommunities.
#'
#' @examples
#' meta <- CommSimABC::rand_meta(20, 30, 500)
#' ecosumstats(meta)
#'
#' metas <- replicate(20, rand_meta(20, 30, 500), simplify = FALSE)
#' ecosumstats(metas)
#'
#' @export

ecosumstats <- function(y) {


  if(is.matrix(y)) {

    y <- y[,colSums(y != 0) > 0] # remove empty columns
    mean.shannon <- mean(vegan::diversity(y, index = "shannon"))
    median.shannon <- stats::median(vegan::diversity(y, index = "shannon"))
    sd.shannon <- stats::sd(vegan::diversity(y, index = "shannon"))
    mean.simpson <- mean(vegan::diversity(y, index = "simpson"))
    median.simpson <- stats::median(vegan::diversity(y, index = "simpson"))
    sd.simpson <- stats::sd(vegan::diversity(y, index = "simpson"))
    mean.evenness <- mean(evenness(y))
    median.evenness <- stats::median(evenness(y))
    sd.evenness <- stats::sd(evenness(y))
    mean.richness <- mean(vegan::specnumber(y))
    median.richness <- stats::median(vegan::specnumber(y))
    sd.richness <- stats::sd(vegan::specnumber(y))
    mean.sites <- mean(colSums(ifelse(y > 0, 1, 0)))
    median.sites <- stats::median(colSums(ifelse(y > 0, 1, 0)))
    sd.sites <- stats::sd(colSums(ifelse(y > 0, 1, 0)))
    beta.bray <- mean(vegan::vegdist(y, method = "bray"))
    beta.bray.med <- stats::median(vegan::vegdist(y, method = "bray"))
    beta.bray.sd <- stats::sd(vegan::vegdist(y, method = "bray"))
    beta.bray.ssid <- ssid2beta(y, method = "bray")
    # alpha.0 <- vegetarian::d(y, lev = 'alpha', q = 0)
    # beta.0 <- vegetarian::d(y, lev = 'beta', q = 0)
    # gamma.0 <- vegetarian::d(y, lev = 'gamma', q = 0)
    # alpha.1 <- vegetarian::d(y, lev = 'alpha', q = 1)
    # beta.1 <- vegetarian::d(y, lev = 'beta', q = 1)
    # gamma.1 <- vegetarian::d(y, lev = 'gamma', q = 1)
    # alpha.2 <- vegetarian::d(y, lev = 'alpha', q = 2)
    # beta.2 <- vegetarian::d(y, lev = 'beta', q = 2)
    # gamma.2 <- vegetarian::d(y, lev = 'gamma', q = 2)
    # alpha.1.weight <- vegetarian::d(y, lev = 'alpha', wts = rowSums(y)/sum(y), q = 1)
    # beta.1.weight <- vegetarian::d(y, lev = 'beta', wts = rowSums(y)/sum(y), q = 1)
    # gamma.1.weight <- vegetarian::d(y, lev = 'gamma', wts = rowSums(y)/sum(y), q = 1)
    # c.score <- EcoSimR::c_score(ifelse(y > 0, 1, 0))
    # c.score.skew <- EcoSimR::c_score_skew(ifelse(y > 0, 1, 0))
    # c.score.var <- EcoSimR::c_score_var(ifelse(y > 0, 1, 0))
    # checkerscore <- EcoSimR::checker(ifelse(y > 0, 1, 0))
    v.ratio <- bipartite::V.ratio(y)
    coherence <- embabs(y)$count
    turnov <- turnover(y)
    morisit <- morisitas(y)

    out <- c(mean.shannon,
             median.shannon,
             sd.shannon,
             mean.simpson,
             median.simpson,
             sd.simpson,
             mean.evenness,
             median.evenness,
             sd.evenness,
             mean.richness,
             median.richness,
             sd.richness,
             mean.sites,
             median.sites,
             sd.sites,
             beta.bray,
             beta.bray.med,
             beta.bray.sd,
             beta.bray.ssid,
             # alpha.0,
             # beta.0,
             # gamma.0,
             # alpha.1,
             # beta.1,
             # gamma.1,
             # alpha.2,
             # beta.2,
             # gamma.2,
             # alpha.1.weight,
             # beta.1.weight,
             # gamma.1.weight,
             # c.score,
             # c.score.skew,
             # c.score.var,
             # checkerscore,
             v.ratio,
             coherence,
             turnov,
             morisit)
    names(out) <- c(
      "mean.shannon",
      "median.shannon",
      "sd.shannon",
      "mean.simpson",
      "median.simpson",
      "sd.simpson",
      "mean.evenness",
      "median.evenness",
      "sd.evenness",
      "mean.richness",
      "median.richness",
      "sd.richness",
      "mean.sites",
      "median.sites",
      "sd.sites",
      "beta.bray",
      "beta.bray.med",
      "beta.bray.sd",
      "beta.bray.ssid",
      # "alpha.0",
      # "beta.0",
      # "gamma.0",
      # "alpha.1",
      # "beta.1",
      # "gamma.1",
      # "alpha.2",
      # "beta.2",
      # "gamma.2",
      # "alpha.1.weight",
      # "beta.1.weight",
      # "gamma.1.weight",
      # "c.score",
      # "c.score.skew",
      # "c.score.var",
      # "checkerscore",
      "v.ratio",
      "coherence",
      "turnov",
      "morisit")

  } else {

    rmempty <- function(x) {
      x <- x[,colSums(x != 0) > 0]
      x
    }

    y <- lapply(y, rmempty) # Remove empty columns
    incs <- lapply(y, function(x) ifelse(x > 0, 1, 0))

    mean.shannon <- sapply(y, function(x) mean(vegan::diversity(x, index = "shannon")))
    median.shannon <- sapply(y, function(x) stats::median(vegan::diversity(x, index = "shannon")))
    sd.shannon <- sapply(y, function(x) stats::sd(vegan::diversity(x, index = "shannon")))
    mean.simpson <- sapply(y, function(x) mean(vegan::diversity(x, index = "simpson")))
    median.simpson <- sapply(y, function(x) stats::median(vegan::diversity(x, index = "simpson")))
    sd.simpson <- sapply(y, function(x) stats::sd(vegan::diversity(x, index = "simpson")))
    mean.evenness <- sapply(y, function(x) mean(evenness(x)))
    median.evenness <- sapply(y, function(x) stats::median(evenness(x)))
    sd.evenness <- sapply(y, function(x) stats::sd(evenness(x)))
    mean.richness <- sapply(y, function(x) mean(vegan::specnumber(x)))
    median.richness <- sapply(y, function(x) stats::median(vegan::specnumber(x)))
    sd.richness <- sapply(y, function(x) stats::sd(vegan::specnumber(x)))
    mean.sites <- sapply(y, function(x) mean(colSums(ifelse(x > 0, 1, 0))))
    median.sites <- sapply(y, function(x) stats::median(colSums(ifelse(x > 0, 1, 0))))
    sd.sites <- sapply(y, function(x) stats::sd(colSums(ifelse(x > 0, 1, 0))))
    beta.bray <- sapply(y, function(x) mean(vegan::vegdist(x, method = "bray")))
    beta.bray.med <- sapply(y, function(x) stats::median(vegan::vegdist(x, method = "bray")))
    beta.bray.sd <- sapply(y, function(x) stats::sd(vegan::vegdist(x, method = "bray")))
    beta.bray.ssid <- sapply(y, function(x) ssid2beta(x, method = "bray"))
    # alpha.0 <- sapply(y, function(x) vegetarian::d(x, lev = 'alpha', q = 0))
    # beta.0 <- sapply(y, function(x) vegetarian::d(x, lev = 'beta', q = 0))
    # gamma.0 <- sapply(y, function(x) vegetarian::d(x, lev = 'gamma', q = 0))
    # alpha.1 <- sapply(y, function(x) vegetarian::d(x, lev = 'alpha', q = 1))
    # beta.1 <- sapply(y, function(x) vegetarian::d(x, lev = 'beta', q = 1))
    # gamma.1 <- sapply(y, function(x) vegetarian::d(x, lev = 'gamma', q = 1))
    # alpha.2 <- sapply(y, function(x) vegetarian::d(x, lev = 'alpha', q = 2))
    # beta.2 <- sapply(y, function(x) vegetarian::d(x, lev = 'beta', q = 2))
    # gamma.2 <- sapply(y, function(x) vegetarian::d(x, lev = 'gamma', q = 2))
    # alpha.1.weight <- sapply(y, function(x) vegetarian::d(x, lev = 'alpha', wts = rowSums(x)/sum(x), q = 1))
    # beta.1.weight <- sapply(y, function(x) vegetarian::d(x, lev = 'beta', wts = rowSums(x)/sum(x), q = 1))
    # gamma.1.weight <- sapply(y, function(x) vegetarian::d(x, lev = 'gamma', wts = rowSums(x)/sum(x), q = 1))
    # c.score <- sapply(incs, function(x) EcoSimR::c_score(x))
    # c.score.skew <- sapply(incs, function(x) EcoSimR::c_score_skew(x))
    # c.score.var <- sapply(incs, function(x) EcoSimR::c_score_var(x))
    # checkerscore <- sapply(incs, function(x) EcoSimR::checker(x))
    v.ratio <- sapply(y, function(x) bipartite::V.ratio(x))
    coherence <- sapply(y, function(x) embabs(x)$count)
    turnov <- sapply(y, function(x) turnover(x))
    morisit <- sapply(y, function(x) morisitas(x))

    out <- data.frame(
      mean.shannon,
      median.shannon,
      sd.shannon,
      mean.simpson,
      median.simpson,
      sd.simpson,
      mean.evenness,
      median.evenness,
      sd.evenness,
      mean.richness,
      median.richness,
      sd.richness,
      mean.sites,
      median.sites,
      sd.sites,
      beta.bray,
      beta.bray.med,
      beta.bray.sd,
      beta.bray.ssid,
      # alpha.0,
      # beta.0,
      # gamma.0,
      # alpha.1,
      # beta.1,
      # gamma.1,
      # alpha.2,
      # beta.2,
      # gamma.2,
      # alpha.1.weight,
      # beta.1.weight,
      # gamma.1.weight,
      # c.score,
      # c.score.skew, 
      # c.score.var,
      # checkerscore,
      v.ratio,
      coherence,
      turnov,
      morisit)
  }
  return(out)
}

#' @describeIn ecosumstats Same as ecosumstats but only for single matrices to allow parallelization
#' @export

ecosumstats2 <- function(y) {
    
    y <- y[,colSums(y != 0) > 0] # remove empty columns
    mean.shannon <- mean(vegan::diversity(y, index = "shannon"))
    median.shannon <- stats::median(vegan::diversity(y, index = "shannon"))
    sd.shannon <- stats::sd(vegan::diversity(y, index = "shannon"))
    mean.simpson <- mean(vegan::diversity(y, index = "simpson"))
    median.simpson <- stats::median(vegan::diversity(y, index = "simpson"))
    sd.simpson <- stats::sd(vegan::diversity(y, index = "simpson"))
    mean.evenness <- mean(evenness(y))
    median.evenness <- stats::median(evenness(y))
    sd.evenness <- stats::sd(evenness(y))
    mean.richness <- mean(vegan::specnumber(y))
    median.richness <- stats::median(vegan::specnumber(y))
    sd.richness <- stats::sd(vegan::specnumber(y))
    mean.sites <- mean(colSums(ifelse(y > 0, 1, 0)))
    median.sites <- stats::median(colSums(ifelse(y > 0, 1, 0)))
    sd.sites <- stats::sd(colSums(ifelse(y > 0, 1, 0)))
    beta.bray <- mean(vegan::vegdist(y, method = "bray"))
    beta.bray.med <- stats::median(vegan::vegdist(y, method = "bray"))
    beta.bray.sd <- stats::sd(vegan::vegdist(y, method = "bray"))
    beta.bray.ssid <- ssid2beta(y, method = "bray")
    # alpha.0 <- vegetarian::d(y, lev = 'alpha', q = 0)
    # beta.0 <- vegetarian::d(y, lev = 'beta', q = 0)
    # gamma.0 <- vegetarian::d(y, lev = 'gamma', q = 0)
    # alpha.1 <- vegetarian::d(y, lev = 'alpha', q = 1)
    # beta.1 <- vegetarian::d(y, lev = 'beta', q = 1)
    # gamma.1 <- vegetarian::d(y, lev = 'gamma', q = 1)
    # alpha.2 <- vegetarian::d(y, lev = 'alpha', q = 2)
    # beta.2 <- vegetarian::d(y, lev = 'beta', q = 2)
    # gamma.2 <- vegetarian::d(y, lev = 'gamma', q = 2)
    # alpha.1.weight <- vegetarian::d(y, lev = 'alpha', wts = rowSums(y)/sum(y), q = 1)
    # beta.1.weight <- vegetarian::d(y, lev = 'beta', wts = rowSums(y)/sum(y), q = 1)
    # gamma.1.weight <- vegetarian::d(y, lev = 'gamma', wts = rowSums(y)/sum(y), q = 1)
    # c.score <- EcoSimR::c_score(ifelse(y > 0, 1, 0))
    # c.score.skew <- EcoSimR::c_score_skew(ifelse(y > 0, 1, 0))
    # c.score.var <- EcoSimR::c_score_var(ifelse(y > 0, 1, 0))
    # checkerscore <- EcoSimR::checker(ifelse(y > 0, 1, 0))
    v.ratio <- bipartite::V.ratio(y)
    coherence <- embabs(y)$count
    turnov <- turnover(y)
    morisit <- morisitas(y)
    
    out <- c(mean.shannon,
             median.shannon,
             sd.shannon,
             mean.simpson,
             median.simpson,
             sd.simpson,
             mean.evenness,
             median.evenness,
             sd.evenness,
             mean.richness,
             median.richness,
             sd.richness,
             mean.sites,
             median.sites,
             sd.sites,
             beta.bray,
             beta.bray.med,
             beta.bray.sd,
             beta.bray.ssid,
             # alpha.0,
             # beta.0,
             # gamma.0,
             # alpha.1,
             # beta.1,
             # gamma.1,
             # alpha.2,
             # beta.2,
             # gamma.2,
             # alpha.1.weight,
             # beta.1.weight,
             # gamma.1.weight,
             # c.score,
             # c.score.skew,
             # c.score.var,
             # checkerscore,
             v.ratio,
             coherence,
             turnov,
             morisit)
    names(out) <- c(
      "mean.shannon",
      "median.shannon",
      "sd.shannon",
      "mean.simpson",
      "median.simpson",
      "sd.simpson",
      "mean.evenness",
      "median.evenness",
      "sd.evenness",
      "mean.richness",
      "median.richness",
      "sd.richness",
      "mean.sites",
      "median.sites",
      "sd.sites",
      "beta.bray",
      "beta.bray.med",
      "beta.bray.sd",
      "beta.bray.ssid",
      # "alpha.0",
      # "beta.0",
      # "gamma.0",
      # "alpha.1",
      # "beta.1",
      # "gamma.1",
      # "alpha.2",
      # "beta.2",
      # "gamma.2",
      # "alpha.1.weight",
      # "beta.1.weight",
      # "gamma.1.weight",
      # "c.score",
      # "c.score.skew",
      # "c.score.var",
      # "checkerscore",
      "v.ratio",
      "coherence",
      "turnov",
      "morisit")
    return(out)
}


#' Write Lists to External File
#'
#' Takes a list and outputs the contents into multiple text files.
#'
#' @param x The list to be written.
#' @param dir Character, directory to which the list will be written.
#' @param sep Character, the separator to use for text delimitation.
#'
#' @details Takes the contents of a list and writes them as separate
#' files into the chosen directory. If a list is nested inside a list,
#' creates a new directory and writes the contents as separate files.
#' NOTE: only works with a single nesting. Multiple nestings will bring an
#' error.
#'
#' @export

write.list <- function(x, dir, sep = "NULL") {

  for(i in 1:length(x)) {
    ###print(paste("i = ", i))
    dex <- x[[i]]
    #print(dex)

    if(is.data.frame(dex)) {
      ###print("dex = dataframe")
      utils::write.table(dex, paste0(dir, names(x)[i]), sep = sep, row.names = F)


    } else if(is.list(dex)) {
      dir.out <- paste0(dir, names(x)[i],"/")
      #print(paste("Directory : ", dir.out))
      dir.create(dir.out)

      for(j in 1:length(dex)) {
        ddex <- dex[[j]]
        #print(paste("ddex = ", ddex))
        if(is.matrix(ddex)) {
          #print("ddex is a matrix")
          utils::write.table(ddex, paste0(dir.out, names(dex)[j]), sep = ",", col.names = T, row.names = F)

        } else if(is.data.frame(ddex)) {
          #print("ddex is a dataframe")
          utils::write.table(ddex, paste0(dir.out, names(ddex)[j]), sep = sep, row.names = F)

        } else {
          ###print("ddex is else")
          write(ddex, paste0(dir, names(ddex)[j]))

        }
      }

    } else if(is.matrix(dex)) {
      ###print("dex = matrix")
      utils::write.table(dex, paste0(dir, names(x)[i]), sep = sep, col.names = T, row.names = F)

    } else {
      ###print("dex = else")
      write(dex, paste0(dir, names(x)[i]))
    }
  }
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
#' This function counts the nubmer of embedded abscences in a binary matrix
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
    #zeros betwen ones in a vector

    valrange <- function(ones) { # Iternal function to make a list of ranges between
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

turnover <- function(web) {
  web <- metacom::OrderMatrix(web, scores = 1, binary = TRUE)
  for (i in 1:ncol(web)) {
    web[min(which(web[, i] == 1)):max(which(web[, i] == 1)), i] <- 1
  }
  D <- vegan::designdist(web, method = "(A-J)*(B-J)", terms = "minimum")
  return(sum(D))
}


morisitas <- function (comm) {
  comm <- metacom::OrderMatrix(comm, scores = 1)
  for (i in 1:ncol(comm)) {
    comm[min(which(comm[, i] == 1)):max(which(comm[,i] == 1)), i] <- 1
  }
  comm <- t(comm)
  M <- 0
  ComBnd <- rep(0, ncol(comm))
  ComBndChi <- 0
  for (i in 1:nrow(comm)) {
    ind1 <- which(comm[i, ] == 1)
    for (j in 1:ncol(comm)) {
      if (min(ind1) == j) {
        ComBnd[j] = ComBnd[j] + 1
      }
      if (max(ind1) == j) {
        ComBnd[j] = ComBnd[j] + 1
      }
    }
  }
  TotComBnd <- (nrow(comm) * 2) - ComBnd[1] - ComBnd[ncol(comm)]
  ExpComBnd <- TotComBnd/(ncol(comm) - 2)
  df <- -1
  for (z in 2:(ncol(comm) - 1)) {
    M <- M + ((ComBnd[z]/TotComBnd) * ((ComBnd[z] - 1)/(TotComBnd -
                                                          1)))
    ComBndChi <- ComBndChi + (((ComBnd[z] - ExpComBnd)^2)/ExpComBnd)
    df <- df + 1
  }
  M <- M * (ncol(comm) - 2)
  M
}

