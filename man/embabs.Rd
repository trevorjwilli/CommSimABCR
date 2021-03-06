% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility_functions.R
\name{embabs}
\alias{embabs}
\title{Count Embedded Absences}
\usage{
embabs(x)
}
\arguments{
\item{x}{A metacommunity matrix with species as columns and sites as rows}
}
\description{
This function counts the nubmer of embedded abscences in a binary matrix
to calculate the coherence of a metacommunity
}
\details{
This function counts the number of embedded absences in a sitexspecies
matrix following Leibold and Mikkelson 2002.
}
\examples{
meta <- rand_meta(10, 15, 100)
embabs(meta)

}
