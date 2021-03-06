% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma.r
\name{rma_known_gamma}
\alias{rma_known_gamma}
\title{Simulate fitness effects with a known number of mutations and Gamma DFE}
\usage{
rma_known_gamma(shape, rate, Ve, k, p_neutral)
}
\arguments{
\item{shape}{numeric,  shape parameter for Gamma}

\item{rate}{numeric, Scale paramater for Gamma}

\item{Ve}{numeric, envrionmental variance}

\item{k}{integer, total number of mutations in each line}

\item{p_neutral, }{proportion of mutations with no fitness effect}

\item{n}{numeric, number of lines to simulate}
}
\value{
w, numeric simulate fitness of each line
}
\description{
This function simulates fitness effects under a model in which the fitness
distribution of mutations takes a Gamma distribution
}
\examples{
k <- stats::rpois(20, 9)
w<- rma_known_gamma(shape=1, rate=20, Ve=0.01, k=k, p_neutral=0.4)
mean(w)
}

