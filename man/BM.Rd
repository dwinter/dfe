% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bm.r
\name{BM}
\alias{BM}
\title{Bateman-mukai estimators for MA data}
\usage{
BM(obs, Ve)
}
\arguments{
\item{obs, }{a vector of ftiness values}

\item{Ve, }{an estimate of within-line ('environmental') variance}
}
\value{
U numeric BM estimator of mutation rate

Ea numeric BM estimator of effect-size
}
\description{
Bateman-mukai estimators for MA data
}
\examples{

low_var <- rma_normal(n=100, a=0.01, Va=1e-9, Ve=1e-4, Ut=10)
BM(low_var, 1e-4)
high_var <- rma_normal(n=100, a=0.01, Va=1e-3, Ve=1e-4, Ut=10)
BM(high_var, 1e-4)
}

