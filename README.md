[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/adegenet.png?branch=master)](https://travis-ci.org/thibautjombart/adegenet)

# adegenet
*adegenet*: a R Package for the Multivariate Analysis of Genetic Markers

This is the development page of the *adegenet* package for the R software. Regular users probably want to look at the http://adegenet.r-forge.r-project.org/, still hosted on Rforge for historical reasons.


### To install the development version of *adegenet*
You will need the package *devtools* to be able to install the devel version of *adegenet*.
To install *devtools*:
```r
install.packages("devtools")
```

To install *adegenet* devel:
```r
library(devtools)
install_github("thibautjombart/adegenet")
library("adegenet")
```

### Tutorials
These are the lastest tutorials (devel version):
- [**basics**](https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-basics.pdf): data handling, basic population genetics, multivariate analysis, IBD, Monmonier
- [**dapc**](https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-dapc.pdf): introduction to the Discriminant Analysis of Principal Components
- [**spca**](https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-spca.pdf): introduction to the spatial Principal Component Analysis
- [**genomics**](https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-genomics.pdf): introduction to the genlight class for large SNP dataset, DAPC, large SNP data simulations
- [**strata**](https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-strata.pdf): introduction to the use of hierarchical clusters with *genind* and *genlight* objects.
