[![Travis-CI Build
Status](https://travis-ci.org/dwinter/dfe.png?branch=master)](https://travis-ci.org/dwinter/dfe)

#Distributions of fitness effects

This is a work-in-progress package, aiming to provide functions for users to
fit existing and new models of the distribution of fitness effects from
data arising from  mutation accumulation experiments

As of Jan 2015, everything here is bleedingly alpha and will almost certainly
change in the future. 

##Package design


###Simulating MA experiments

Functions starting `rma_*()` simulate the fitness effects arising from a mutation
accumulation study, in which the fitness-effects are distributed according to a
`*` distribution. Current options are:


```
rma_normal()
rma_gamma()
rma_FGM()
```

`rma_FGM()` simulates mutations under a paramaterization of Fisher's Geometric Model (by default fitness is determined by the squared distance form the origin, user-defined fitness functions are allowed). 

###Likelihood for observed data

Functions starting `dma_*()` calculate likelihood (densities) under the normal or
gamma distributed models:

```
dma_gamma()
dma_normal()
```

There are also ML fitting functions, `fit_ma*()`, which are... a work in progress

###Miscellaneous functions

`BM()` calculates Bateman-Mukai (method of moments) estimators, `moments_gamma()`
calculates the mean and variance of a given Gamma distribution, `moments_FGM()`
estimates the mean, variance, skewness and proportion of beneficial mutations
via simulation. 
