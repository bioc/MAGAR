# methQTL-package
[![Build Status](https://travis-ci.org/MPIIComputationalEpigenetics/methQTL-package.svg?branch=master)](https://travis-ci.org/MPIIComputationalEpigenetics/methQTL-package)

![](pictures/logo.png)

## Short description
R package for computing DNA methylation quantitative trait loci (methQTL) from DNA methylation and matched genotyping data. The tool employs a linear modeling strategy to determine siginifant associations between genotype variation and changes in DNA methylation of particular CpGs.

## Installation
The package can directly be installed in a R session using the ```devtools``` package
```r
if(!requireNamespace("devtools")) install.packages("devtools")
 devtools::install_github("MPIIComputationalEpigenetics/methQTL-package")
```

## Documentation
The methQTL R-package is documented [here](vignettes/methQTL.md)

## Contact
You can contact [Michael Scherer](mailto:mscherer@mpi-inf.mpg.de) for reporting bugs, feature requests or questions.
