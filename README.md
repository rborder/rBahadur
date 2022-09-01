# rBahadur

> Efficient simulation of genotype / phenotype data under
> assortative mating by generating Bahadur order-2
> multivariate Bernoulli distributed random variates.

## Features

* Multivariate Bernoulli (MVB) distribution samplers
  * `rb_dplr`: generate Bahadur order-2 MVB variates with diagonal-plus-low-rank correlation structures
  * `rb_unstr`: generate Bahadur order-2 MVB variates with arbitrary correlation structures
* Assortative mating  modeling tools
  * Compute equilibrium parameters under univariate AM
  * Generate genotype / phenotype data given initial conditions


## Installation

We recommend installation using the `install-github` function provided by the [`remotes` library](https://github.com/r-lib/remotes):

```r
remotes::install-github("rborder/rBahadur")
```

## Usage
