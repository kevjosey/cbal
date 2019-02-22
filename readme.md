cbal: Covariate Balancing Weights
=================================

An R package for finding covariate balancing weights for causal inference.
These functions find the solution to a convex optimization problem using
Bregman distances as the criterion function. The balancing weights are evaluated
as the generalized projection of the sampling weights into an intersection of hyperplanes
that form the basis of a linear subspace. This linear subspace constrains the moments for the 
emprical distribution of the covariates to be identical between treatment groups.

In its current form, cbal supports three Bregman distances and requires a binary treatment assignment.
In the future, this package will be extended to allow for other criterion distance functions and multivalued
treatment assignments.

## Vignettes

There is one vignette included with the package. Additional details will be added to this vignette.
For now, the vignette is empty but soon you should be able to use:

```r
vignette('cbal-pkg', package = 'cbal')
```

## Installing cbal
Options for installing cbal:

1. Install the developmental version from github.  This will require you to have
   [devtools](https://github.com/hadley/devtools) installed and, if you are
   using Windows, you'll need
   [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed as well.

```
library(devtools)
install_github("kevjosey/cbal", build_vignettes = TRUE)
```

2. Clone the repo and use `GNU make`

```bash
make install
```

3. Go to the [release page](https://github.com/kevjosey/cbal/releases) and down
   load the tar.gz file of the version you want to install.

  * Install from the command line

```bash
R CMD INSTALL cbal_<version>.tar.gz
```

  * Within R

```r
install.packages(<path_to_file>, repos = NULL, type = "source")
```