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
#' pairdist <- dist(xy)
#' make_mig(pairdist, 20, .1)
#'
#' @export

make_mig <- function(site.arrange, max.dist, tot) {
  if(class(site.arrange) == "dist") {
    n.sites <- attr(site.arrange, "Size")
    mig.dist <- site.arrange
  } else {
    n.sites <- (length(site.arrange[,1]))
    mig.dist <- dist(site.arrange, upper = T, diag = T) # Calculates the euclidean distances between sites
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
#' simpson diversity index, Pielou's evenness, richness, beta diversity (whittakers and Sorensons),
#' and the number of sites occupied for species for communities in a metacommunity matrix. It also
#' calculats the mean and median C-score and the v-ratio.
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

    mean.shannon <- mean(vegan::diversity(y, index = "shannon"))
    median.shannon <- median(vegan::diversity(y, index = "shannon"))
    sd.shannon <- sd(vegan::diversity(y, index = "shannon"))
    mean.simpson <- mean(vegan::diversity(y, index = "simpson"))
    median.simpson <- median(vegan::diversity(y, index = "simpson"))
    sd.simpson <- sd(vegan::diversity(y, index = "simpson"))
    mean.evenness <- mean(evenness(y))
    median.evenness <- median(evenness(y))
    sd.evenness <- sd(evenness(y))
    mean.richness <- mean(vegan::specnumber(y))
    median.richness <- median(vegan::specnumber(y))
    sd.richness <- sd(vegan::specnumber(y))
    beta.whit <- mean(vegan::betadiver(y, method = 1))
    beta.whit.med <- median(vegan::betadiver(y, method = 1))
    beta.whit.sd <- sd(vegan::betadiver(y, method = 1))
    beta.sor <- mean(vegan::betadiver(y, method = 2))
    beta.sor.med <- median(vegan::betadiver(y, method = 2))
    beta.sor.sd <- sd(vegan::betadiver(y, method = 2))
    mean.c.score <- bipartite::C.score(y, normalise = T)
    median.c.score <- bipartite::C.score(y, normalise = T, FUN = median)
    v.ratio <- bipartite::V.ratio(y)
    mean.sites <- mean(colSums(ifelse(y > 0, 1, 0)))
    median.sites <- median(colSums(ifelse(y > 0, 1, 0)))
    sd.sites <- sd(colSums(ifelse(y > 0, 1, 0)))

    out <- c(mean.shannon, median.shannon, sd.shannon, mean.simpson, median.simpson, sd.simpson,
             mean.evenness, median.evenness, sd.evenness, mean.richness, median.richness,
             sd.richness, beta.whit, beta.whit.med, beta.whit.sd, beta.sor, beta.sor.med,
             beta.sor.sd, mean.c.score, median.c.score, v.ratio, mean.sites, median.sites, sd.sites)
    names(out) <- c("mean.shannon", "median.shannon", "sd.shannon", "mean.simpson", "median.simpson", "sd.simpson",
                    "mean.evenness", "median.evenness", "sd.evenness", "mean.richness", "median.richness",
                    "sd.richness", "beta.whit", "beta.whit.med", "beta.whit.sd", "beta.sor", "beta.sor.med",
                    "beta.sor.sd", "mean.c.score", "median.c.score", "v.ratio", "mean.sites", "median.sites", "sd.sites")

  } else {

    mean.shannon <- sapply(y, function(x) mean(vegan::diversity(x, index = "shannon")))
    median.shannon <- sapply(y, function(x) median(vegan::diversity(x, index = "shannon")))
    sd.shannon <- sapply(y, function(x) sd(vegan::diversity(x, index = "shannon")))
    mean.simpson <- sapply(y, function(x) mean(vegan::diversity(x, index = "simpson")))
    median.simpson <- sapply(y, function(x) median(vegan::diversity(x, index = "simpson")))
    sd.simpson <- sapply(y, function(x) sd(vegan::diversity(x, index = "simpson")))
    mean.evenness <- sapply(y, function(x) mean(evenness(x)))
    median.evenness <- sapply(y, function(x) median(evenness(x)))
    sd.evenness <- sapply(y, function(x) sd(evenness(x)))
    mean.richness <- sapply(y, function(x) mean(vegan::specnumber(x)))
    median.richness <- sapply(y, function(x) median(vegan::specnumber(x)))
    sd.richness <- sapply(y, function(x) sd(vegan::specnumber(x)))
    beta.whit <- sapply(y, function(x) mean(vegan::betadiver(x, method = 1)))
    beta.whit.med <- sapply(y, function(x) median(vegan::betadiver(x, method = 1)))
    beta.whit.sd <- sapply(y, function(x) sd(vegan::betadiver(x, method = 1)))
    beta.sor <- sapply(y, function(x) mean(vegan::betadiver(x, method = 2)))
    beta.sor.med <- sapply(y, function(x) median(vegan::betadiver(x, method = 2)))
    beta.sor.sd <- sapply(y, function(x) sd(vegan::betadiver(x, method = 2)))
    mean.c.score <- sapply(y, function(x) bipartite::C.score(x, normalise = F))
    median.c.score <- sapply(y, function(x) bipartite::C.score(x, normalise = F, FUN = median))
    v.ratio <- sapply(y, function(x) bipartite::V.ratio(x))
    mean.sites <- sapply(y, function(x) mean(colSums(ifelse(x > 0, 1, 0))))
    median.sites <- sapply(y, function(x) median(colSums(ifelse(x > 0, 1, 0))))
    sd.sites <- sapply(y, function(x) sd(colSums(ifelse(x > 0, 1, 0))))

    out <- data.frame(mean.shannon, median.shannon, sd.shannon, mean.simpson, median.simpson, sd.simpson,
                      mean.evenness, median.evenness, sd.evenness, mean.richness, median.richness,
                      sd.richness, beta.whit, beta.whit.med, beta.whit.sd, beta.sor, beta.sor.med,
                      beta.sor.sd, mean.c.score, median.c.score, v.ratio, mean.sites, median.sites, sd.sites)
  }
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
      write.table(dex, paste0(dir, names(x)[i]), sep = sep, row.names = F)


    } else if(is.list(dex)) {
      dir.out <- paste0(dir, names(x)[i],"/")
      #print(paste("Directory : ", dir.out))
      dir.create(dir.out)

      for(j in 1:length(dex)) {
        ddex <- dex[[j]]
        #print(paste("ddex = ", ddex))
        if(is.matrix(ddex)) {
          #print("ddex is a matrix")
          write.table(ddex, paste0(dir.out, names(dex)[j]), sep = ",", col.names = T, row.names = F)

        } else if(is.data.frame(ddex)) {
          #print("ddex is a dataframe")
          write.table(ddex, paste0(dir.out, names(ddex)[j]), sep = sep, row.names = F)

        } else {
          ###print("ddex is else")
          write(ddex, paste0(dir, names(ddex)[j]))

        }
      }

    } else if(is.matrix(dex)) {
      ###print("dex = matrix")
      write.table(dex, paste0(dir, names(x)[i]), sep = sep, col.names = T, row.names = F)

    } else {
      ###print("dex = else")
      write(dex, paste0(dir, names(x)[i]))
    }
  }
}





