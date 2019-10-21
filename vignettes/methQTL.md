---
title: "methQTL package"
author: "Michael Scherer"
date: "October 20, 2019"
output: pdf_document
  sansfont: Georgia
---



# Introduction

This vignette describes the *methQTL* R-package (https://github.com/MPIIComputationalEpigenetics/methQTL-package) available from GitHub. The package uses DNA methylation data obtained using the Illumina BeadArrays, and genotyping data from Illumina genotyping microarrays or whole genome sequencing to compute methylation quantitative trait loci (methQTL). Using a linear modeling strategy, *methQTL* computes statistically significant interactions between single nucleotide polymorphisms (SNPs) and changes in the DNA methylation state of individual CpGs.

# Installation

The package can be directly installed from GitHub, after installing the *devtools* package.


```r
if(!requireNamespace("devtools")) install.packages("devtools")
```

```
## Loading required namespace: devtools
```

```r
if(!requireNamespace("methQTL")) devtools::install_github("MPIIComputationalEpigenetics/methQTL-package")
```

```
## Loading required namespace: methQTL
```

```
## Warning in .recacheSubclasses(def@className, def, env): undefined subclass
## "HDF5Matrix" of class "matrixOrHDF"; definition not updated
```

```
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
```

```
## Registered S3 method overwritten by 'openssl':
##   method      from
##   print.bytes Rcpp
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```
## Warning: no function found corresponding to methods exports from 'RnBeads'
## for: 'samples'
```

```r
suppressPackageStartupMessages(library(methQTL))
```

# Input data

The *methQTL* package needs two types of data as input: DNA methylation data obtained using the Illumina Infinium BeadArrays and genotyping data obtained using genotyping microarrays or whole genome sequencing.

## DNA methylation data

The package uses the widely used [*RnBeads*](http://rnbeads.org/) software package for DNA methylation data import. It requires raw IDAT files as input, since critical preprocessing is to be perfomed within the package. In addition to the raw methylation data, a sample annotation sheet, specifying the samples to be analyzed needs to be provided. The sheet contains a line for each sample and looks as follows:


```bash
SampleID,age,sex,barcode
Sample_1,14,f,209054857842_R01C01
Sample_2,42,f,209054857842_R02C01
Sample_3,45,m,209054857842_R03C01
```

For further details on the import process, we refer to the [RnBeads vignette](http://bioconductor.org/packages/release/bioc/vignettes/RnBeads/inst/doc/RnBeads.pdf). Most importantly, analysis options need to be specified for the import and preprocessing modules of *RnBeads*. *methQTL* provides a default setting, which is available in *extdata/rnbeads_options.xml*. You can use this file as a template for your own setting and then specify it to the methQTL package:


```r
opts <- rnb.xml2options(system.file("extdata/rnbeads_options.xml",package="methQTL"))
rnb.options(identifiers.column="SampleID")
xml.fi <- file.path(getwd(),"rnbeads_options.xml")
cat(rnb.options2xml(),file=xml.fi)
qtl.setOption(rnbeads.options = xml.fi)
```

## Genotyping data

We require that genotyping data has already been processed by *plink* and is available either in the form of binary *.bed*, *.bim* and *.fam* files, as *.ped* and *.map*, or as variant calling files (*.vcf*). For further processing, we use the command line tool *plink*, which comes with this package, but is only applicable on Linux systems. For Windows and MacOS users, please install the *plink* tool from [here](https://www.cog-genomics.org/plink/1.9/) and specify it using the option ```plink.path```. The sample identifier specified earlier also needs to match the sample IDs of the genotype calls.

## Perform data import

The ```do.import``` function requires the paths to the respective genotyping and DNA methylation data, as well as a sample annotation sheet as discussed earlier. You'll have to specify the paths to the corresponding *IDAT* and *plink* files. Additionally, you have to specify the sample identifier column in the sample annotation sheet that determines the samples in both the genotyping and DNA methylation data. For larger files, we recommend to activate the option to store large matrices on disk rather than in main memory ```"hdf5dump"```.


```r
idat.dir <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
plink.dir <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/"
anno.sheet <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_IL_IPC_example.tsv"
qtl.setOption(hdf5dump=TRUE)
imp.data <- do.import(data.location = c(idat.dir=idat.dir,plink.dir=plink.dir),
                      s.anno = anno.sheet,
                      s.id.col = "ind_IPC",
                      tab.sep = "\t")
```

```
## 2019-10-21 14:14:29     1.3  STATUS STARTED Import methQTL data
## 2019-10-21 14:14:29     1.3  STATUS     STARTED Processing genotyping data
## 2019-10-21 14:14:44     1.5  STATUS     COMPLETED Processing genotyping data
## 2019-10-21 14:14:44     1.5  STATUS     STARTED Processing DNA methylation data
```

```
## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.
```

```
## 2019-10-21 14:14:48     1.5  STATUS         STARTED Loading Data from IDAT Files
## 2019-10-21 14:14:49     1.5    INFO             Added column barcode to the provided sample annotation table
## 2019-10-21 14:14:49     1.5    INFO             Detected platform: MethylationEPIC
## 2019-10-21 14:15:06     2.0  STATUS         COMPLETED Loading Data from IDAT Files
## 2019-10-21 14:16:22     2.2  STATUS         STARTED Preprocessing
## 2019-10-21 14:16:22     2.2    INFO             Number of cores: 1
## 2019-10-21 14:16:22     2.2  STATUS             STARTED Filtering Procedures I
## 2019-10-21 14:16:23     2.2  STATUS                 STARTED Removal of SNP-enriched Sites
## 2019-10-21 14:16:23     2.2  STATUS                     Removed 139721 sites using SNP criterion "any"
## 2019-10-21 14:16:24     2.2  STATUS                     Saved removed sites to /local/tmp/RtmpFFIS9w/rnbeads_preprocessing/preprocessing_data/removed_sites_snp.csv
## 2019-10-21 14:16:24     2.2  STATUS                     Added a corresponding section to the report
## 2019-10-21 14:16:24     2.2  STATUS                 COMPLETED Removal of SNP-enriched Sites
## 2019-10-21 14:16:24     2.2  STATUS                 STARTED Removal of Cross-reactive Probes
## 2019-10-21 14:16:24     2.2  STATUS                     Removed 34264 sites
## 2019-10-21 14:16:24     2.2  STATUS                     Saved removed sites to /local/tmp/RtmpFFIS9w/rnbeads_preprocessing/preprocessing_data/removed_sites_cross_reactive.csv
## 2019-10-21 14:16:24     2.2  STATUS                     Added a corresponding section to the report
## 2019-10-21 14:16:24     2.2  STATUS                 COMPLETED Removal of Cross-reactive Probes
## 2019-10-21 14:16:25     2.2    INFO                 Working with a p-value threshold of 0.05
## 2019-10-21 14:16:25     2.4  STATUS                 STARTED Greedycut
## 2019-10-21 14:16:45     2.4  STATUS                     Calculated a total of 1055 iterations
## 2019-10-21 14:16:45     2.4    INFO                     Optimal number of iterations is 1055
```

```
## 2019-10-21 14:16:50     2.4  STATUS                     Created ROC plot
```

```
## 2019-10-21 14:16:55     2.4  STATUS                     Created line plots for matrix dimensions and other statistics
## 2019-10-21 14:16:55     2.4  STATUS                     Saved removed sites to /local/tmp/RtmpFFIS9w/rnbeads_preprocessing/preprocessing_data/removed_sites_greedycut.csv
## 2019-10-21 14:16:55     2.4  STATUS                 COMPLETED Greedycut
## 2019-10-21 14:16:55     2.4  STATUS                 Retained 8 samples and 691856 sites
## 2019-10-21 14:16:55     2.4  STATUS             COMPLETED Filtering Procedures I
## 2019-10-21 14:16:55     2.4  STATUS             STARTED Summary of Filtering Procedures I
## 2019-10-21 14:16:56     2.4  STATUS                 Created summary table of removed sites, samples and unreliable measurements
```

```
## 2019-10-21 14:16:57     2.4  STATUS                 Added summary table of removed and retained items
## 2019-10-21 14:16:57     2.4    INFO                 Subsampling 866895 sites for plotting density distributions
## 2019-10-21 14:16:58     2.5  STATUS                 Constructed sequences of removed and retained methylation values
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

```
## 2019-10-21 14:17:09     2.7  STATUS                 Added comparison between removed and retained beta values
## 2019-10-21 14:17:09     2.7  STATUS             COMPLETED Summary of Filtering Procedures I
## 2019-10-21 14:17:09     2.7  STATUS             STARTED Manipulating the object
## 2019-10-21 14:18:05     2.5  STATUS                 Removed 175039 sites (probes)
## 2019-10-21 14:18:05     2.5    INFO                 Retained 691856 sites and 8 samples
## 2019-10-21 14:18:05     2.5  STATUS             COMPLETED Manipulating the object
## 2019-10-21 14:18:05     2.5  STATUS             STARTED Normalization Procedure
## 2019-10-21 14:18:40     3.6  STATUS                 Performed normalization with method wm.dasen
## 2019-10-21 14:19:33     3.2  STATUS                 Performed normalization with method "wm.dasen"
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

```
## 2019-10-21 14:19:44     3.6  STATUS                 Added comparison between non-normalized and normalized beta values
```

```
## 2019-10-21 14:19:46     3.6  STATUS                 Added histogram of observed beta shifts (magnitude of correction)
```

```
## 2019-10-21 14:19:47     3.6  STATUS                 Added 2D histogram of observed beta values and shifts
## 2019-10-21 14:19:47     3.6  STATUS                 Added normalization section
## 2019-10-21 14:19:47     3.6  STATUS             COMPLETED Normalization Procedure
## 2019-10-21 14:19:47     3.6  STATUS             STARTED Filtering Procedures II
## 2019-10-21 14:19:49     3.6  STATUS                 STARTED Probe Context Removal
## 2019-10-21 14:19:49     3.6  STATUS                     Removed 1196 probe(s) having not acceptable context
## 2019-10-21 14:19:49     3.6  STATUS                     Saved removed sites to /local/tmp/RtmpFFIS9w/rnbeads_preprocessing/preprocessing_data/removed_sites_context.csv
## 2019-10-21 14:19:49     3.6  STATUS                     Added a corresponding section to the report
## 2019-10-21 14:19:49     3.6  STATUS                 COMPLETED Probe Context Removal
## 2019-10-21 14:19:49     3.6  STATUS                 STARTED Removal of Sites on Sex Chromosomes
## 2019-10-21 14:19:49     3.6  STATUS                     Removed 16598 site(s) on sex chromosomes
## 2019-10-21 14:19:49     3.6  STATUS                     Saved removed sites to /local/tmp/RtmpFFIS9w/rnbeads_preprocessing/preprocessing_data/removed_sites_sex.csv
## 2019-10-21 14:19:49     3.6  STATUS                     Added a corresponding section to the report
## 2019-10-21 14:19:49     3.6  STATUS                 COMPLETED Removal of Sites on Sex Chromosomes
## 2019-10-21 14:19:49     3.6  STATUS                 STARTED Missing Value Removal
## 2019-10-21 14:19:49     3.6  STATUS                     Using a sample quantile threshold of 0
## 2019-10-21 14:19:49     3.6  STATUS                     Removed 44 site(s) with too many missing values
## 2019-10-21 14:19:49     3.6  STATUS                     Saved removed sites to /local/tmp/RtmpFFIS9w/rnbeads_preprocessing/preprocessing_data/removed_sites_na.csv
```

```
## 2019-10-21 14:19:53     3.6  STATUS                     Added a corresponding section to the report
## 2019-10-21 14:19:53     3.6  STATUS                 COMPLETED Missing Value Removal
## 2019-10-21 14:19:53     3.6  STATUS                 Retained 8 samples and 674018 sites
## 2019-10-21 14:19:53     3.6  STATUS             COMPLETED Filtering Procedures II
## 2019-10-21 14:19:53     3.6  STATUS             STARTED Summary of Filtering Procedures II
## 2019-10-21 14:19:53     3.6  STATUS                 Created summary table of removed sites, samples and unreliable measurements
```

```
## 2019-10-21 14:19:54     3.6  STATUS                 Added summary table of removed and retained items
## 2019-10-21 14:19:54     3.6    INFO                 Subsampling 691856 sites for plotting density distributions
## 2019-10-21 14:19:55     3.7  STATUS                 Constructed sequences of removed and retained methylation values
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

```
## 2019-10-21 14:20:02     3.5  STATUS                 Added comparison between removed and retained beta values
## 2019-10-21 14:20:02     3.5  STATUS             COMPLETED Summary of Filtering Procedures II
## 2019-10-21 14:20:02     3.5  STATUS             STARTED Manipulating the object
## 2019-10-21 14:20:55     3.5  STATUS                 Removed 17838 sites (probes)
## 2019-10-21 14:20:55     3.5    INFO                 Retained 674018 sites and 8 samples
## 2019-10-21 14:20:55     3.5  STATUS             COMPLETED Manipulating the object
## 2019-10-21 14:20:55     3.5    INFO             No missing values present, imputation skipped
## 2019-10-21 14:20:55     3.5  STATUS         COMPLETED Preprocessing
## 2019-10-21 14:20:59     3.5  STATUS     COMPLETED Processing DNA methylation data
## 2019-10-21 14:20:59     3.5  STATUS COMPLETED Import methQTL data
```

# methQTL calling

Although *methQTL* conceptually splits the methQTL calling into two steps ((i) compute correlation block, (ii) call methQTL per correlation block), only a single function call is needed. The function only requires the input ```methQTLInput``` object produced in the previous step, but further options, such as covariates and the p-value cutoff can be directly specified as a function parameter.


```r
meth.qtl.res <- do.methQTL(imp.data)
```

```
## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.

## Warning in structure(x$children, class = "XMLNodeList"): Calling 'structure(NULL, *)' is deprecated, as NULL cannot have attributes.
##   Consider 'structure(list(), *)' instead.
```

```
## 2019-10-21 14:21:03     3.5  STATUS STARTED Imputation procedure knn 
## 2019-10-21 14:21:12     3.3  STATUS COMPLETED Imputation procedure knn 
## 
## 2019-10-21 14:21:14     3.3  STATUS STARTED Computing methQTLs
## 2019-10-21 14:21:14     3.3  STATUS     STARTED Computing methQTL for chromosome chr1
## 2019-10-21 14:21:14     3.3  STATUS         STARTED Compute correlation blocks
## 2019-10-21 14:21:14     3.3    INFO             Split workload, since facing 66108 CpGs (Maximum is 40000 )
## 2019-10-21 14:21:14     3.3  STATUS             STARTED Compute correlation blocks
## 2019-10-21 14:21:14     3.3  STATUS                 STARTED Compute correlation matrix
## 2019-10-21 14:21:46    19.6  STATUS                 COMPLETED Compute correlation matrix
## 2019-10-21 14:27:58    14.2  STATUS                 STARTED Compute correlation distance
```

```
## Warning in match.fun(.Generic)(a): NaNs produced
```

```
## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced
```

```
## 2019-10-21 14:35:55    12.3  STATUS                 COMPLETED Compute correlation distance
## 2019-10-21 14:40:28    13.0  STATUS                 STARTED Compute pairwise distances
## 2019-10-21 14:40:52    13.4  STATUS                 COMPLETED Compute pairwise distances
## 2019-10-21 14:41:43    10.3  STATUS                 STARTED Weight distances
## 2019-10-21 14:47:17    28.8  STATUS                 COMPLETED Weight distances
## 2019-10-21 14:47:18    18.4  STATUS                 STARTED Compute graph
## 2019-10-21 14:48:26    34.7  STATUS                 COMPLETED Compute graph
## 2019-10-21 14:48:26    34.7  STATUS                 STARTED Compute clustering
## 2019-10-21 14:48:29    34.7  STATUS                 COMPLETED Compute clustering
## 2019-10-21 14:48:29    34.7  STATUS             COMPLETED Compute correlation blocks
## 2019-10-21 14:48:29    34.7  STATUS             STARTED Compute correlation blocks
## 2019-10-21 14:48:29    34.7  STATUS                 STARTED Compute correlation matrix
## 2019-10-21 14:48:58    51.0  STATUS                 COMPLETED Compute correlation matrix
## 2019-10-21 14:54:54    64.3  STATUS                 STARTED Compute correlation distance
```

```
## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced

## Warning in match.fun(.Generic)(a): NaNs produced
```

```
## 2019-10-21 15:02:45    37.8  STATUS                 COMPLETED Compute correlation distance
## 2019-10-21 15:07:16    57.9  STATUS                 STARTED Compute pairwise distances
## 2019-10-21 15:07:39    18.6  STATUS                 COMPLETED Compute pairwise distances
## 2019-10-21 15:08:28    14.1  STATUS                 STARTED Weight distances
## 2019-10-21 15:13:50    18.5  STATUS                 COMPLETED Weight distances
## 2019-10-21 15:13:51    16.6  STATUS                 STARTED Compute graph
## 2019-10-21 15:14:59    32.9  STATUS                 COMPLETED Compute graph
## 2019-10-21 15:14:59    32.9  STATUS                 STARTED Compute clustering
## 2019-10-21 15:15:01    32.9  STATUS                 COMPLETED Compute clustering
## 2019-10-21 15:15:01    32.9  STATUS             COMPLETED Compute correlation blocks
## 2019-10-21 15:15:01    32.9  STATUS             STARTED Compute methQTL per correlation block
```

```
## Error in summary(lm.model)$coefficients["SNP", "Pr(>|t|)"]: subscript out of bounds
```

We will now present the two steps of the methQTL calling procedure in more detail.

## Compute CpG correlation blocks

Since neighboring CpGs are often highly correlated, using each CpG idenpendently as a potential methQTL candidate leads to many redundant results. We thus aimed to approximate *DNA methylation haplotypes* by determining highly correlated CpGs in close vicinity. The procedure itself is split into six steps, and is performed for each chromosome independently:

1. Compute the (Pearson) correlation matrix between all CpGs
2. Construct the distance matrix from the correlation matrix
3. Discard all interactions with a correlation lower than a given threshold (option ```cluster.cor.threshold```, default: 0.2)
4. Weight the distance according to the genomic distance between the two CpGs with a Gaussian (options ```standard.deviation.gauss```, default: 5000)
5. Discard all interactions further than the option ```absolute.distance.cutoff``` (default: 1,000,000)
6. Compute the Louvain clustering on the undirected, weighted graph induced by the distance matrix

## Call methQTL per correlation block

From the list of correlation blocks, *methQTL* computes methQTL interactions with all SNPs on the same chromosome. The process is split into three steps:

1. Compute a representative CpG per correlation block, as specified with the option ```representative.cpg.computation``` (default: *row.medians*).
2. Discard all SNPs that are further than ```absolute.distance.cutoff``` (default: 1,000,000) away from the representative CpG
3. Call methQTL by using linear models. Multiple options of methQTL calling are available and can be selected via the option ```linear.model.type``` (default: *classical.linear*).

In the latest stage, potential covariates can be specified using the option *sel.covariates*. We recommend to include at least *age* and *sex* as covariates, as they have a strong influence on the DNA methylation pattern.

# Advanced configuration

## Employ methQTL on a scientific compute cluster

*methQTL* can automatically distribute jobs across a high performance compute cluster, which has been setup using the Sun Grid Engine (SGE) technology. You can pass the option ```cluster.submit``` to ```do.methQTL``` and thus activate the compute cluster submission. Note that you'll also have to specify a path to an executable *Rscript* and potentially specify resource requirements using the option setting ```cluster.config```.


```r
qtl.setOption(cluster.config = c(h_vmem="60G",mem_free="20G"))
qtl.setOption(rscript.path = "/usr/bin/Rscript")
meth.qtl.res <- do.methQTL(meth.qtl = imp.data,
                          cluster.submit = T)
```
