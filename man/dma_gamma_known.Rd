% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dma_gamma_known}
\alias{dma_gamma_known}
\title{Density function of gamma dfe}
\usage{
dma_gamma_known(obs, shape, rate, Ve, k, p_neutral, log = TRUE)
}
\arguments{
\item{obs, }{numeric, observed fitnesses}

\item{shape, }{numeric, shape of gamma distribution}

\item{rate, }{numeric, rate of gamma distribution}

\item{Ve, }{numeric, experimental variance}

\item{k, }{numeric, vector of knonw mutation counts}

\item{p_neutral}{numeric,proportion of all mutations that have no effect}

\item{log}{logical return log-liklihood (defaults to true)}
}
\value{
numeric (log-) liklihood of the specfified model and data
}
\description{
Density function of gamma dfe
}
\examples{
set.seed(123)
mu <- rpois(50, 10)
w <- rma_known_gamma(shape=1, rate=25, Ve=0.01, k=mu,p_neutral=0.7) 
dma_gamma_known(w, shape=1, rate=25, Ve=0.01, k=mu, p_neutral=0.75)
dma_gamma_known(w, shape=1, rate=25, Ve=0.01, k=mu, p_neutral=0.65)
}

