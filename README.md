# methQTL-package
[![Build Status](https://travis-ci.org/MPIIComputationalEpigenetics/methQTL-package.svg?branch=master)](https://travis-ci.org/MPIIComputationalEpigenetics/methQTL-package)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![](pictures/logo.png)

## Short description
*methQTL* is an R package for computing DNA methylation quantitative trait loci (methQTL) from DNA methylation and matched genotyping data. The tool employs a two-stage approach to compute methQTL:

- CpG correlation block are computed to mimic DNA methylation haplotypes. The correlation blocks are computed from CpGs that are highly correlated across the samples.
- For each correlation block, a representative CpG (tag-CpG) is computed and used in a linear modeling strategy to determine significant associations between genotype variation and changes in DNA methylation of particular CpGs.

The package furthermore uses established software tools for handling DNA methylation and genotyping data, including [*RnBeads*](https://rnbeads.org), [*PLINK*](http://zzz.bwh.harvard.edu/plink/), and [*CRLMM*](https://www.bioconductor.org/packages/release/bioc/html/crlmm.html). Thus, raw DNA methylation data (IDAT, BED files) and genotyping data (IDAT, PLINK files) can be used as input.

## Installation
The package needs some Bioconductor dependencies, most notably RnBeads, which can be directly installed from the [*RnBeads* website](https://rnbeads.org). Atfer instaling RnBeads, the *methQTL* package can be installed in a R session using the ```devtools``` package
```r
source("https://rnbeads.org/data/install.R")
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("MPIIComputationalEpigenetics/methQTL-package")
```

## Documentation
The methQTL R-package is documented [here](vignettes/methQTL.md)

## Contact
You can contact [Michael Scherer](mailto:mscherer@mpi-inf.mpg.de) for reporting bugs, feature requests or questions.
