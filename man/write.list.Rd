% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility_functions.R
\name{write.list}
\alias{write.list}
\title{Write Lists to External File}
\usage{
write.list(x, dir, sep = "NULL")
}
\arguments{
\item{x}{The list to be written.}

\item{dir}{Character, directory to which the list will be written.}

\item{sep}{Character, the separator to use for text delimitation.}
}
\description{
Takes a list and outputs the contents into multiple text files.
}
\details{
Takes the contents of a list and writes them as separate
files into the chosen directory. If a list is nested inside a list,
creates a new directory and writes the contents as separate files.
NOTE: only works with a single nesting. Multiple nestings will bring an
error.
}
