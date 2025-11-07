# CommSimABCR 

## About The Project

This repository contains the development version for the R package CommSimABCR. 
CommSimABCR is a package to run stochastic metacommunity simulation models based on a modified
haploid Moran process. In short, the model is a stochastic birth-death forward simulation model 
that includes parameters for selection, migration, frequency dependence, and community sizes. 
This package also calculates standard community ecology summary statistics for use in 
downstream approximate Bayesian computation (ABC) processes, like model selection or 
parameter estimation. For more details see [link to preprint]

## Table of Contents

  * [About](##AboutTheProject)
  * [Getting Started](##GettingStarted)
    * [Installation](##Installation)
  * [Usage](##Usage)
  * [License](## License)
    

## Getting Started

### Installation

Currently, the package is only available through github, though there are future plans
to host it on cran. The source code can be downloaded from githab natively or using
the `devtools` package.

```R
devtools::install_github(trevorjwilli/CommSimABCR)
```

## Usage

Vignettes are not currently available but will be shortly. For now, please see
examples in the package documentation.

## License

Distributed under the [GPL-3.0 license](LICENSE)

## Contact

Trevor Williams - <trevorjwilli@gmail.com>

