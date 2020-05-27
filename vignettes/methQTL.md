---
title: "methQTL package"
author: "Michael Scherer"
date: "May 25, 2020"
header-includes:
  \usepackage{helvet}
  \hypersetup{colorlinks,
              urlcolor=blue}
output:
  pdf_document:
    number_sections: true
    latex_engine: xelatex
vignette:
  \VignetteEngine{knitr::knitr}
bibliography: biblio.bib  
---



# Introduction

This vignette describes the *methQTL* R-package (https://github.com/MPIIComputationalEpigenetics/methQTL-package) available from GitHub. The package uses DNA methylation data obtained using the Illumina BeadArrays, and genotyping data from Illumina genotyping microarrays or whole genome sequencing to compute methylation quantitative trait loci (methQTL). Using a linear modeling strategy, *methQTL* computes statistically significant interactions between single nucleotide polymorphisms (SNPs) and changes in the DNA methylation state of individual CpGs. DNA methylation values at single CpGs are first summarized into correlation blocks, and a representative of this correlation block is used for methQTL calling.

# Installation

The package can be directly installed from GitHub, after installing the *devtools* package.


```r
if(!requireNamespace("devtools")) install.packages("devtools")
```

```
## Loading required namespace: devtools
```

```r
if(!requireNamespace("methQTL")){
  devtools::install_github("MPIIComputationalEpigenetics/methQTL-package")
}  
```

```
## Loading required namespace: methQTL
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```
## Warning: no function found corresponding to methods exports from 'RnBeads' for:
## 'samples'
```

```r
suppressPackageStartupMessages(library(methQTL))
```

# Input data

The *methQTL* package requires two types of data as input: DNA methylation data obtained using the Illumina Infinium BeadArrays and genotyping data obtained using genotyping microarrays or whole genome sequencing.

## DNA methylation data

The package utilizes the widely used [*RnBeads*](http://rnbeads.org/) software package for DNA methylation data import. It requires raw IDAT files as input, since critical preprocessing is to be perfomed within the package. In addition to the raw methylation data, a sample annotation sheet specifying the samples to be analyzed needs to be provided. The sheet contains a line for each sample and looks as follows:


```bash
SampleID,age,sex,barcode
Sample_1,14,f,209054857842_R01C01
Sample_2,42,f,209054857842_R02C01
Sample_3,45,m,209054857842_R03C01
```

For further details on the import process, we refer to the [RnBeads vignette](http://bioconductor.org/packages/release/bioc/vignettes/RnBeads/inst/doc/RnBeads.pdf). Most importantly, analysis options need to be specified for the import and preprocessing modules of *RnBeads*. *methQTL* provides a default setting, which is available in *extdata/rnbeads_options.xml*. You can use this file as a template for your own setting and then specify it to the methQTL package:


```r
opts <- rnb.xml2options(system.file("extdata/rnbeads_options.xml",package="methQTL"))
rnb.options(identifiers.column="ind_IPC")
xml.fi <- file.path(getwd(),"rnbeads_options.xml")
cat(rnb.options2xml(),file=xml.fi)
qtl.setOption(rnbeads.options = xml.fi)
```

To redefine the correlation blocks, we allow including additional information such as genome-wide segmentation of the methylation landscape (see option ```use.segmentation``` and function ```qtl.run.segmentation```), and also function annotation according to the Ensembl regulatory build[@Zerbino2015].

## Genotyping data

### PLINK files

The package supports data that has already been processed by *plink* and is available either in the form of binary *.bed*, *.bim* and *.fam* files, as *.ped* and *.map*, as variant calling files (*.vcf*), or as imputed files in the dosage format (*.dos*). For further processing, we use the command line tool *plink*, which comes with this package, but is only applicable on Linux systems. For Windows and MacOS users, please install the *plink* tool from [here](https://www.cog-genomics.org/plink/1.9/) and specify it using the option ```plink.path```. The sample identifier specified earlier also needs to match the sample IDs of the genotype calls.

### IDAT files

The package also support raw IDAT files and uses the [CRLMM](https://www.bioconductor.org/packages/release/bioc/html/crlmm.html) R-package, together with PLINK to perform data import. The package requires a single sample annotation sheet in the format described in the [DNA methylation data] section. In addition to the column names specified above, a column named *GenoSentrixPosition* has to be added, which specifies the IDAT file IDs.


```bash
SampleID,age,sex,barcode,GenoSentrixPosition
Sample_1,14,f,209054857842_R01C01,9701756058_R05C01
Sample_2,42,f,209054857842_R02C01,9701756058_R07C01
Sample_3,45,m,209054857842_R03C01,9742011016_R04C01
```

### Imputation

Illumina SNP BeadArray data is typically imputed before further analysis, and the package allows for imputation through the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!). In order to be able to perform computation on the server, an account is required. After the account is created, one has to request an API token in the user settings and specify it to the package using the option ```imputation.user.token``` option. During the imputation process, the package will stall for a while an wait for entering the password send via e-mail to the user account. The imputation process has to be split according to chromosomes, which is why multiple e-mails will be send to the account, and the imputation process can take up to several days. However, after imputation, the imputed data will be available as PLINK files, such that the imputation has to be performed only once. For preprocessing the data for upload to the imputation server, the package required the [bgzip](http://www.htslib.org/doc/bgzip.html) and [tabix](http://www.htslib.org/doc/tabix.html) tools from the [htslib](http://www.htslib.org/) package. Also see further options to configure the imputation jobs (see the Michigan Imputation Server [documentation](https://imputationserver.readthedocs.io/en/latest/) for further information):


```r
qtl.setOption(
  impute.geno.data=TRUE,
  imputation.reference.panel="apps@hrc-r1.1",
  imputation.phasing.method="shapeit",
	imputation.population="eur"
)
```

## Perform data import

The ```do.import``` function requires the paths to the respective genotyping and DNA methylation data, as well as a sample annotation sheet as discussed earlier. You'll have to specify the paths to the corresponding *IDAT* and *plink* files. Additionally, you have to specify the sample identifier column in the sample annotation sheet that determines the samples in both the genotyping and DNA methylation data. For larger files, we recommend to activate the option to store large matrices on disk rather than in main memory (```hdf5dump```).


```r
idat.dir <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
plink.dir <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/"
anno.sheet <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_IL_IPC_example.tsv"
qtl.setOption(hdf5dump=TRUE)
imp.data <- do.import(data.location = c(idat.dir=idat.dir,geno.dir=plink.dir),
                      s.anno = anno.sheet,
                      s.id.col = "ind_IPC",
                      tab.sep = "\t")
```

```
## 2020-05-26 21:11:46     1.3  STATUS STARTED Import methQTL data
## 2020-05-26 21:11:46     1.3  STATUS     STARTED Processing genotyping data
## 2020-05-26 21:11:54     1.7    INFO         Loading system default for option 'plink.path'
## 2020-05-26 21:12:06     2.2  STATUS         STARTED Compute genotype PCA
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-26 21:12:06     2.2  STATUS         COMPLETED Compute genotype PCA
## 2020-05-26 21:12:08     2.2  STATUS     COMPLETED Processing genotyping data
## 2020-05-26 21:12:08     2.2  STATUS     STARTED Processing DNA methylation data
## 2020-05-26 21:12:09     2.2  STATUS         STARTED Loading Data from IDAT Files
## 2020-05-26 21:12:09     2.2    INFO             Added column barcode to the provided sample annotation table
## 2020-05-26 21:12:09     2.2    INFO             Detected platform: MethylationEPIC
## 2020-05-26 21:12:32     2.3  STATUS         COMPLETED Loading Data from IDAT Files
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:13:55     2.6  STATUS         STARTED Preprocessing
## 2020-05-26 21:13:55     2.6    INFO             Number of cores: 1
## 2020-05-26 21:13:55     2.6  STATUS             STARTED Filtering Procedures I
## 2020-05-26 21:13:56     2.6  STATUS                 STARTED Removal of SNP-enriched Sites
## 2020-05-26 21:13:56     2.6  STATUS                     Removed 139721 sites using SNP criterion "any"
## 2020-05-26 21:13:57     2.6  STATUS                     Saved removed sites to /local/tmp/Rtmpp3LH3d/rnbeads_preprocessing/preprocessing_data/removed_sites_snp.csv
## 2020-05-26 21:13:57     2.6  STATUS                     Added a corresponding section to the report
## 2020-05-26 21:13:57     2.6  STATUS                 COMPLETED Removal of SNP-enriched Sites
## 2020-05-26 21:13:57     2.6  STATUS                 STARTED Removal of Cross-reactive Probes
## 2020-05-26 21:13:57     2.6  STATUS                     Removed 34264 sites
## 2020-05-26 21:13:57     2.6  STATUS                     Saved removed sites to /local/tmp/Rtmpp3LH3d/rnbeads_preprocessing/preprocessing_data/removed_sites_cross_reactive.csv
## 2020-05-26 21:13:57     2.6  STATUS                     Added a corresponding section to the report
## 2020-05-26 21:13:57     2.6  STATUS                 COMPLETED Removal of Cross-reactive Probes
## 2020-05-26 21:13:57     2.6    INFO                 Working with a p-value threshold of 0.05
## 2020-05-26 21:13:59     2.7  STATUS                 STARTED Greedycut
## 2020-05-26 21:14:20     2.8  STATUS                     Calculated a total of 1055 iterations
## 2020-05-26 21:14:20     2.8    INFO                     Optimal number of iterations is 1055
```

```
## 2020-05-26 21:14:24     2.8  STATUS                     Created ROC plot
```

```
## 2020-05-26 21:14:28     2.8  STATUS                     Created line plots for matrix dimensions and other statistics
## 2020-05-26 21:14:28     2.8  STATUS                     Saved removed sites to /local/tmp/Rtmpp3LH3d/rnbeads_preprocessing/preprocessing_data/removed_sites_greedycut.csv
## 2020-05-26 21:14:28     2.8  STATUS                 COMPLETED Greedycut
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:14:28     2.8  STATUS                 Retained 8 samples and 691856 sites
## 2020-05-26 21:14:28     2.8  STATUS             COMPLETED Filtering Procedures I
## 2020-05-26 21:14:28     2.8  STATUS             STARTED Summary of Filtering Procedures I
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:14:28     2.8  STATUS                 Created summary table of removed sites, samples and unreliable measurements
```

```
## 2020-05-26 21:14:30     2.8  STATUS                 Added summary table of removed and retained items
## 2020-05-26 21:14:30     2.8    INFO                 Subsampling 866895 sites for plotting density distributions
## 2020-05-26 21:14:30     2.8  STATUS                 Constructed sequences of removed and retained methylation values
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

```
## 2020-05-26 21:14:44     2.8  STATUS                 Added comparison between removed and retained beta values
## 2020-05-26 21:14:44     2.8  STATUS             COMPLETED Summary of Filtering Procedures I
## 2020-05-26 21:14:44     2.8  STATUS             STARTED Manipulating the object
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:15:43     2.8  STATUS                 Removed 175039 sites (probes)
## 2020-05-26 21:15:43     2.8    INFO                 Retained 691856 sites and 8 samples
## 2020-05-26 21:15:43     2.8  STATUS             COMPLETED Manipulating the object
## 2020-05-26 21:15:43     2.8  STATUS             STARTED Normalization Procedure
```

```
## No methods found in package 'RSQLite' for request: 'dbListFields' when loading 'lumi'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:16:25     3.4  STATUS                 Performed normalization with method wm.dasen
## 2020-05-26 21:17:21     3.4  STATUS                 Performed normalization with method "wm.dasen"
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

```
## 2020-05-26 21:17:35     3.5  STATUS                 Added comparison between non-normalized and normalized beta values
```

```
## 2020-05-26 21:17:37     3.6  STATUS                 Added histogram of observed beta shifts (magnitude of correction)
```

```
## 2020-05-26 21:17:38     3.6  STATUS                 Added 2D histogram of observed beta values and shifts
## 2020-05-26 21:17:38     3.7  STATUS                 Added normalization section
## 2020-05-26 21:17:38     3.7  STATUS             COMPLETED Normalization Procedure
## 2020-05-26 21:17:38     3.7  STATUS             STARTED Filtering Procedures II
## 2020-05-26 21:17:40     3.5  STATUS                 STARTED Probe Context Removal
## 2020-05-26 21:17:40     3.5  STATUS                     Removed 1196 probe(s) having not acceptable context
## 2020-05-26 21:17:40     3.5  STATUS                     Saved removed sites to /local/tmp/Rtmpp3LH3d/rnbeads_preprocessing/preprocessing_data/removed_sites_context.csv
## 2020-05-26 21:17:40     3.5  STATUS                     Added a corresponding section to the report
## 2020-05-26 21:17:40     3.5  STATUS                 COMPLETED Probe Context Removal
## 2020-05-26 21:17:40     3.5  STATUS                 STARTED Removal of Sites on Sex Chromosomes
## 2020-05-26 21:17:40     3.5  STATUS                     Removed 16598 site(s) on sex chromosomes
## 2020-05-26 21:17:41     3.5  STATUS                     Saved removed sites to /local/tmp/Rtmpp3LH3d/rnbeads_preprocessing/preprocessing_data/removed_sites_sex.csv
## 2020-05-26 21:17:41     3.5  STATUS                     Added a corresponding section to the report
## 2020-05-26 21:17:41     3.5  STATUS                 COMPLETED Removal of Sites on Sex Chromosomes
## 2020-05-26 21:17:41     3.5  STATUS                 STARTED Missing Value Removal
## 2020-05-26 21:17:41     3.5  STATUS                     Using a sample quantile threshold of 0
## 2020-05-26 21:17:41     3.5  STATUS                     Removed 44 site(s) with too many missing values
## 2020-05-26 21:17:41     3.5  STATUS                     Saved removed sites to /local/tmp/Rtmpp3LH3d/rnbeads_preprocessing/preprocessing_data/removed_sites_na.csv
```

```
## 2020-05-26 21:17:46     3.2  STATUS                     Added a corresponding section to the report
## 2020-05-26 21:17:46     3.2  STATUS                 COMPLETED Missing Value Removal
## 2020-05-26 21:17:46     3.2  STATUS                 Retained 8 samples and 674018 sites
## 2020-05-26 21:17:46     3.2  STATUS             COMPLETED Filtering Procedures II
## 2020-05-26 21:17:46     3.2  STATUS             STARTED Summary of Filtering Procedures II
## 2020-05-26 21:17:46     3.2  STATUS                 Created summary table of removed sites, samples and unreliable measurements
```

```
## 2020-05-26 21:17:48     3.2  STATUS                 Added summary table of removed and retained items
## 2020-05-26 21:17:48     3.2    INFO                 Subsampling 691856 sites for plotting density distributions
## 2020-05-26 21:17:48     3.2  STATUS                 Constructed sequences of removed and retained methylation values
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

```
## 2020-05-26 21:17:57     3.3  STATUS                 Added comparison between removed and retained beta values
## 2020-05-26 21:17:57     3.3  STATUS             COMPLETED Summary of Filtering Procedures II
## 2020-05-26 21:17:57     3.3  STATUS             STARTED Manipulating the object
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:18:51     3.3  STATUS                 Removed 17838 sites (probes)
## 2020-05-26 21:18:51     3.3    INFO                 Retained 674018 sites and 8 samples
## 2020-05-26 21:18:51     3.3  STATUS             COMPLETED Manipulating the object
## 2020-05-26 21:18:52     3.3    INFO             No missing values present, imputation skipped
## 2020-05-26 21:18:52     3.3  STATUS         COMPLETED Preprocessing
```

```
## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
##   - in 'x': chrX, chrY
##   - in 'y': chr23, chr25
##   Make sure to always combine/compare objects based on the same reference
##   genome (use suppressWarnings() to suppress this warning).
```

```
## 2020-05-26 21:19:22     3.3  STATUS         STARTED Removing 870 CpGs overlapping with SNPs
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e067712910b.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e063350e87e.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e064007a007.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e064848f7cc.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e06601b4dbb.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e062097224a.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e0676fbd50f.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e061c9fbb88.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## opening ff /local/tmp/Rtmpp3LH3d/ff8e066dd2bf06.ff
```

```
## Found more than one class "ff_matrix" in cache; using the first, from namespace 'oligoClasses'
```

```
## Also defined by 'RnBeads'
```

```
## 2020-05-26 21:20:17     3.3  STATUS         COMPLETED Removing 870 CpGs overlapping with SNPs
## 2020-05-26 21:20:21     3.3  STATUS     COMPLETED Processing DNA methylation data
## 2020-05-26 21:20:21     3.3  STATUS COMPLETED Import methQTL data
```

For imputed data, no further processing is performed on the genotyping data and the dosage values are used as they are:


```r
idat.dir <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
geno.dir <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO_IMPUTED/"
anno.sheet <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_IL_IPC_example.tsv"
qtl.setOption(hdf5dump=TRUE)
imp.data <- do.import(data.location = c(idat.dir=idat.dir,geno.dir=geno.dir),
                      s.anno = anno.sheet,
                      s.id.col = "ind_IPC",
                      tab.sep = "\t")
```

Please note that the ```recode.allele.frequencies``` option specifies, if, according to the cohort analyzed, SNP reference and alternative allele are to be recoded according to the allele frequencies found. Alternatively, a path to a local version of dbSNP[@Sheery2001] can be provided through ```db.snp.ref```, and reference/alternative allele information will be automatically parsed from the database. This is especially crucial, if imputation is to be performed, since the Michigan Imputation Server is sensitive to reference mismatches. ```recode.allele.frequencies``` and ```db.snp.ref``` are mutually exclusive options.

# methQTL calling

Although *methQTL* conceptually splits the methQTL calling into two steps ((i) compute correlation block, (ii) call methQTL per correlation block), only a single function call is needed. The function only requires the input ```methQTLInput``` object produced in the previous step, but further options, such as covariates and the p-value cutoff can be directly specified as a function parameter, or as global parameters using ```?qtl.setOption```.


```r
meth.qtl.res <- do.methQTL(imp.data)
```

```
## 2020-05-26 21:20:21     3.3    INFO Loading default option setting
## 2020-05-26 21:20:22     3.3  STATUS STARTED Imputation procedure knn 
## 2020-05-26 21:20:28     3.4  STATUS COMPLETED Imputation procedure knn 
## 
## 2020-05-26 21:20:30     3.5  STATUS STARTED Computing methQTLs
## 2020-05-26 21:20:30     3.5  STATUS     STARTED Computing methQTL for chromosome chr1
## 2020-05-26 21:20:31     3.5  STATUS         STARTED Compute correlation blocks
## 2020-05-26 21:20:31     3.5    INFO             Split workload, since facing 66034 CpGs (Maximum is 40000 )
## 2020-05-26 21:20:31     3.5  STATUS             STARTED Compute correlation blocks
## 2020-05-26 21:20:31     3.5  STATUS                 STARTED Compute correlation matrix
## 2020-05-26 21:21:04    19.7  STATUS                 COMPLETED Compute correlation matrix
## 2020-05-26 21:31:34    13.6  STATUS                 STARTED Compute pairwise distances
## 2020-05-26 21:31:58    13.5  STATUS                 COMPLETED Compute pairwise distances
## 2020-05-26 21:32:47    10.4  STATUS                 STARTED Weight distances
## 2020-05-26 21:38:17    27.4  STATUS                 COMPLETED Weight distances
## 2020-05-26 21:38:19    16.0  STATUS                 STARTED Compute graph
## 2020-05-26 21:39:11    32.3  STATUS                 COMPLETED Compute graph
## 2020-05-26 21:39:11    32.3  STATUS                 STARTED Compute clustering
## 2020-05-26 21:39:13    32.3  STATUS                 COMPLETED Compute clustering
## 2020-05-26 21:39:13    32.3  STATUS             COMPLETED Compute correlation blocks
## 2020-05-26 21:39:13    32.3  STATUS             STARTED Compute correlation blocks
## 2020-05-26 21:39:13    32.3  STATUS                 STARTED Compute correlation matrix
## 2020-05-26 21:39:46    48.5  STATUS                 COMPLETED Compute correlation matrix
## 2020-05-26 21:50:21    38.8  STATUS                 STARTED Compute pairwise distances
## 2020-05-26 21:50:43    20.6  STATUS                 COMPLETED Compute pairwise distances
## 2020-05-26 21:51:27    10.4  STATUS                 STARTED Weight distances
## 2020-05-26 21:56:50    37.2  STATUS                 COMPLETED Weight distances
## 2020-05-26 21:56:51    18.5  STATUS                 STARTED Compute graph
## 2020-05-26 21:57:44    34.8  STATUS                 COMPLETED Compute graph
## 2020-05-26 21:57:44    34.8  STATUS                 STARTED Compute clustering
## 2020-05-26 21:57:45    34.8  STATUS                 COMPLETED Compute clustering
## 2020-05-26 21:57:45    34.8  STATUS             COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-26 21:57:46    34.8  STATUS             STARTED Compute methQTL per correlation block
## 2020-05-26 21:57:46    34.8  STATUS                 STARTED Setting up Multicore
## 2020-05-26 21:57:46    34.8    INFO                     Using 1 cores
## 2020-05-26 21:57:46    34.8  STATUS                 COMPLETED Setting up Multicore
## 2020-05-26 21:57:46    34.8    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:20    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:20    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:20    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:20    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:21    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:21    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:21    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:22    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:22    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:22    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 22:34:22    10.4    INFO                 No SNP closer than 500000
## 2020-05-26 23:07:49    10.4  STATUS             COMPLETED Compute methQTL per correlation block
## 2020-05-26 23:07:49    10.4  STATUS         COMPLETED Compute correlation blocks
## 2020-05-26 23:07:49    10.4  STATUS         STARTED Computing methQTL for chromosome chr2
## 2020-05-26 23:07:49    10.4  STATUS             STARTED Compute correlation blocks
## 2020-05-26 23:07:49    10.4    INFO                 Split workload, since facing 52088 CpGs (Maximum is 40000 )
## 2020-05-26 23:07:49    10.4  STATUS                 STARTED Compute correlation blocks
## 2020-05-26 23:07:49    10.4  STATUS                     STARTED Compute correlation matrix
## 2020-05-26 23:08:07    20.5  STATUS                     COMPLETED Compute correlation matrix
## 2020-05-26 23:14:38    25.5  STATUS                     STARTED Compute pairwise distances
## 2020-05-26 23:14:51    10.4  STATUS                     COMPLETED Compute pairwise distances
## 2020-05-26 23:15:22    10.4  STATUS                     STARTED Weight distances
## 2020-05-26 23:18:49    19.7  STATUS                     COMPLETED Weight distances
## 2020-05-26 23:18:51    13.5  STATUS                     STARTED Compute graph
## 2020-05-26 23:19:20    36.8  STATUS                     COMPLETED Compute graph
## 2020-05-26 23:19:20    36.8  STATUS                     STARTED Compute clustering
## 2020-05-26 23:19:21    36.8  STATUS                     COMPLETED Compute clustering
## 2020-05-26 23:19:21    36.8  STATUS                 COMPLETED Compute correlation blocks
## 2020-05-26 23:19:21    36.8  STATUS                 STARTED Compute correlation blocks
## 2020-05-26 23:19:21    36.8  STATUS                     STARTED Compute correlation matrix
## 2020-05-26 23:19:39    41.9  STATUS                     COMPLETED Compute correlation matrix
## 2020-05-26 23:26:16    30.6  STATUS                     STARTED Compute pairwise distances
## 2020-05-26 23:26:30    14.9  STATUS                     COMPLETED Compute pairwise distances
## 2020-05-26 23:26:59    10.4  STATUS                     STARTED Weight distances
## 2020-05-26 23:30:20    18.3  STATUS                     COMPLETED Weight distances
## 2020-05-26 23:30:21    15.4  STATUS                     STARTED Compute graph
## 2020-05-26 23:30:52    37.6  STATUS                     COMPLETED Compute graph
## 2020-05-26 23:30:52    37.6  STATUS                     STARTED Compute clustering
## 2020-05-26 23:30:53    37.6  STATUS                     COMPLETED Compute clustering
## 2020-05-26 23:30:53    37.6  STATUS                 COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-26 23:30:53    37.6  STATUS                 STARTED Compute methQTL per correlation block
## 2020-05-26 23:30:53    37.6  STATUS                     STARTED Setting up Multicore
## 2020-05-26 23:30:53    37.6    INFO                         Using 1 cores
## 2020-05-26 23:30:53    37.6  STATUS                     COMPLETED Setting up Multicore
## 2020-05-27 00:38:58    10.4  STATUS                 COMPLETED Compute methQTL per correlation block
## 2020-05-27 00:38:58    10.4  STATUS             COMPLETED Compute correlation blocks
## 2020-05-27 00:38:58    10.4  STATUS             STARTED Computing methQTL for chromosome chr3
## 2020-05-27 00:38:58    10.4  STATUS                 STARTED Compute correlation blocks
## 2020-05-27 00:38:58    10.4  STATUS                     STARTED Compute correlation matrix
## 2020-05-27 00:39:46    34.0  STATUS                     COMPLETED Compute correlation matrix
## 2020-05-27 00:54:50    22.2  STATUS                     STARTED Compute pairwise distances
## 2020-05-27 00:55:19    18.8  STATUS                     COMPLETED Compute pairwise distances
## 2020-05-27 00:56:24    14.1  STATUS                     STARTED Weight distances
## 2020-05-27 01:05:42    36.2  STATUS                     COMPLETED Weight distances
## 2020-05-27 01:05:43    25.9  STATUS                     STARTED Compute graph
## 2020-05-27 01:16:33    49.6  STATUS                     COMPLETED Compute graph
## 2020-05-27 01:16:33    49.6  STATUS                     STARTED Compute clustering
## 2020-05-27 01:16:34    49.6  STATUS                     COMPLETED Compute clustering
## 2020-05-27 01:16:34    49.6  STATUS                 COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-27 01:16:35    49.6  STATUS                 STARTED Compute methQTL per correlation block
## 2020-05-27 01:16:35    49.6  STATUS                     STARTED Setting up Multicore
## 2020-05-27 01:16:35    49.6    INFO                         Using 1 cores
## 2020-05-27 01:16:35    49.6  STATUS                     COMPLETED Setting up Multicore
## 2020-05-27 02:09:24    10.4  STATUS                 COMPLETED Compute methQTL per correlation block
## 2020-05-27 02:09:24    10.4  STATUS             COMPLETED Computing methQTL for chromosome chr3
## 2020-05-27 02:09:24    10.4  STATUS             STARTED Computing methQTL for chromosome chr4
## 2020-05-27 02:09:24    10.4  STATUS                 STARTED Compute correlation blocks
## 2020-05-27 02:09:24    10.4  STATUS                     STARTED Compute correlation matrix
## 2020-05-27 02:09:48    22.9  STATUS                     COMPLETED Compute correlation matrix
## 2020-05-27 02:17:50    11.4  STATUS                     STARTED Compute pairwise distances
## 2020-05-27 02:18:07    17.9  STATUS                     COMPLETED Compute pairwise distances
## 2020-05-27 02:18:48    10.4  STATUS                     STARTED Weight distances
## 2020-05-27 02:22:54    25.3  STATUS                     COMPLETED Weight distances
## 2020-05-27 02:22:55    16.6  STATUS                     STARTED Compute graph
## 2020-05-27 02:23:32    45.7  STATUS                     COMPLETED Compute graph
## 2020-05-27 02:23:32    45.7  STATUS                     STARTED Compute clustering
## 2020-05-27 02:23:33    45.7  STATUS                     COMPLETED Compute clustering
## 2020-05-27 02:23:33    45.7  STATUS                 COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-27 02:23:33    45.7  STATUS                 STARTED Compute methQTL per correlation block
## 2020-05-27 02:23:33    45.7  STATUS                     STARTED Setting up Multicore
## 2020-05-27 02:23:33    45.7    INFO                         Using 1 cores
## 2020-05-27 02:23:33    45.7  STATUS                     COMPLETED Setting up Multicore
## 2020-05-27 03:05:30    10.4  STATUS                 COMPLETED Compute methQTL per correlation block
## 2020-05-27 03:05:30    10.4  STATUS             COMPLETED Computing methQTL for chromosome chr4
## 2020-05-27 03:05:30    10.4  STATUS             STARTED Computing methQTL for chromosome chr5
## 2020-05-27 03:05:30    10.4  STATUS                 STARTED Compute correlation blocks
## 2020-05-27 03:05:30    10.4  STATUS                     STARTED Compute correlation matrix
## 2020-05-27 03:06:07    29.5  STATUS                     COMPLETED Compute correlation matrix
## 2020-05-27 03:18:31    39.2  STATUS                     STARTED Compute pairwise distances
## 2020-05-27 03:18:54    15.3  STATUS                     COMPLETED Compute pairwise distances
## 2020-05-27 03:19:48    11.8  STATUS                     STARTED Weight distances
## 2020-05-27 03:26:09    30.7  STATUS                     COMPLETED Weight distances
## 2020-05-27 03:26:10    21.4  STATUS                     STARTED Compute graph
## 2020-05-27 03:29:52    40.6  STATUS                     COMPLETED Compute graph
## 2020-05-27 03:29:52    40.6  STATUS                     STARTED Compute clustering
## 2020-05-27 03:29:53    40.6  STATUS                     COMPLETED Compute clustering
## 2020-05-27 03:29:53    40.6  STATUS                 COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-27 03:29:54    40.6  STATUS                 STARTED Compute methQTL per correlation block
## 2020-05-27 03:29:54    40.6  STATUS                     STARTED Setting up Multicore
## 2020-05-27 03:29:54    40.6    INFO                         Using 1 cores
## 2020-05-27 03:29:54    40.6  STATUS                     COMPLETED Setting up Multicore
## 2020-05-27 03:47:02    10.4    INFO                     No SNP closer than 500000
## 2020-05-27 03:47:03    10.4    INFO                     No SNP closer than 500000
## 2020-05-27 04:16:39    10.4  STATUS                 COMPLETED Compute methQTL per correlation block
## 2020-05-27 04:16:39    10.4  STATUS             COMPLETED Computing methQTL for chromosome chr5
## 2020-05-27 04:16:39    10.4  STATUS             STARTED Computing methQTL for chromosome chr6
## 2020-05-27 04:16:39    10.4  STATUS                 STARTED Compute correlation blocks
## 2020-05-27 04:16:39    10.4    INFO                     Split workload, since facing 41991 CpGs (Maximum is 40000 )
## 2020-05-27 04:16:39    10.4  STATUS                     STARTED Compute correlation blocks
## 2020-05-27 04:16:39    10.4  STATUS                         STARTED Compute correlation matrix
## 2020-05-27 04:16:51    17.0  STATUS                         COMPLETED Compute correlation matrix
## 2020-05-27 04:20:52    25.9  STATUS                         STARTED Compute pairwise distances
## 2020-05-27 04:21:02    10.5  STATUS                         COMPLETED Compute pairwise distances
## 2020-05-27 04:21:21    10.4  STATUS                         STARTED Weight distances
## 2020-05-27 04:23:29    18.2  STATUS                         COMPLETED Weight distances
## 2020-05-27 04:23:30    11.9  STATUS                         STARTED Compute graph
## 2020-05-27 04:23:50    25.1  STATUS                         COMPLETED Compute graph
## 2020-05-27 04:23:50    25.1  STATUS                         STARTED Compute clustering
## 2020-05-27 04:23:52    25.1  STATUS                         COMPLETED Compute clustering
## 2020-05-27 04:23:52    25.1  STATUS                     COMPLETED Compute correlation blocks
## 2020-05-27 04:23:52    25.1  STATUS                     STARTED Compute correlation blocks
## 2020-05-27 04:23:52    25.1  STATUS                         STARTED Compute correlation matrix
## 2020-05-27 04:24:07    16.9  STATUS                         COMPLETED Compute correlation matrix
## 2020-05-27 04:28:11    25.0  STATUS                         STARTED Compute pairwise distances
## 2020-05-27 04:28:20    10.5  STATUS                         COMPLETED Compute pairwise distances
## 2020-05-27 04:28:39    10.4  STATUS                         STARTED Weight distances
## 2020-05-27 04:30:47    15.0  STATUS                         COMPLETED Weight distances
## 2020-05-27 04:30:49    13.5  STATUS                         STARTED Compute graph
## 2020-05-27 04:31:09    26.6  STATUS                         COMPLETED Compute graph
## 2020-05-27 04:31:09    26.6  STATUS                         STARTED Compute clustering
## 2020-05-27 04:31:10    26.6  STATUS                         COMPLETED Compute clustering
## 2020-05-27 04:31:10    26.6  STATUS                     COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-27 04:31:11    26.6  STATUS                     STARTED Compute methQTL per correlation block
## 2020-05-27 04:31:11    26.6  STATUS                         STARTED Setting up Multicore
## 2020-05-27 04:31:11    26.6    INFO                             Using 1 cores
## 2020-05-27 04:31:11    26.6  STATUS                         COMPLETED Setting up Multicore
## 2020-05-27 05:25:07    10.4  STATUS                     COMPLETED Compute methQTL per correlation block
## 2020-05-27 05:25:07    10.4  STATUS                 COMPLETED Compute correlation blocks
## 2020-05-27 05:25:07    10.4  STATUS                 STARTED Computing methQTL for chromosome chr7
## 2020-05-27 05:25:08    10.4  STATUS                     STARTED Compute correlation blocks
## 2020-05-27 05:25:08    10.4  STATUS                         STARTED Compute correlation matrix
## 2020-05-27 05:25:49    30.4  STATUS                         COMPLETED Compute correlation matrix
## 2020-05-27 05:38:48    20.4  STATUS                         STARTED Compute pairwise distances
## 2020-05-27 05:39:14    16.3  STATUS                         COMPLETED Compute pairwise distances
## 2020-05-27 05:40:10    12.3  STATUS                         STARTED Weight distances
## 2020-05-27 05:46:54    30.6  STATUS                         COMPLETED Weight distances
## 2020-05-27 05:46:56    22.3  STATUS                         STARTED Compute graph
## 2020-05-27 05:53:39    42.4  STATUS                         COMPLETED Compute graph
## 2020-05-27 05:53:39    42.4  STATUS                         STARTED Compute clustering
## 2020-05-27 05:53:41    42.4  STATUS                         COMPLETED Compute clustering
## 2020-05-27 05:53:41    42.4  STATUS                     COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-27 05:53:42    42.4  STATUS                     STARTED Compute methQTL per correlation block
## 2020-05-27 05:53:42    42.4  STATUS                         STARTED Setting up Multicore
## 2020-05-27 05:53:42    42.4    INFO                             Using 1 cores
## 2020-05-27 05:53:42    42.4  STATUS                         COMPLETED Setting up Multicore
## 2020-05-27 06:35:51    10.4  STATUS                     COMPLETED Compute methQTL per correlation block
## 2020-05-27 06:35:51    10.4  STATUS                 COMPLETED Computing methQTL for chromosome chr7
## 2020-05-27 06:35:51    10.4  STATUS                 STARTED Computing methQTL for chromosome chr8
## 2020-05-27 06:35:51    10.4  STATUS                     STARTED Compute correlation blocks
## 2020-05-27 06:35:51    10.4  STATUS                         STARTED Compute correlation matrix
## 2020-05-27 06:36:16    24.2  STATUS                         COMPLETED Compute correlation matrix
## 2020-05-27 06:45:13    10.8  STATUS                         STARTED Compute pairwise distances
## 2020-05-27 06:45:32    19.5  STATUS                         COMPLETED Compute pairwise distances
## 2020-05-27 06:46:09    10.4  STATUS                         STARTED Weight distances
## 2020-05-27 06:50:47    22.8  STATUS                         COMPLETED Weight distances
## 2020-05-27 06:50:49    17.1  STATUS                         STARTED Compute graph
## 2020-05-27 06:51:34    30.9  STATUS                         COMPLETED Compute graph
## 2020-05-27 06:51:34    30.9  STATUS                         STARTED Compute clustering
## 2020-05-27 06:51:35    30.9  STATUS                         COMPLETED Compute clustering
## 2020-05-27 06:51:35    30.9  STATUS                     COMPLETED Compute correlation blocks
```

```
## Saving 7 x 7 in image
```

```
## 2020-05-27 06:51:36    30.9  STATUS                     STARTED Compute methQTL per correlation block
## 2020-05-27 06:51:36    30.9  STATUS                         STARTED Setting up Multicore
## 2020-05-27 06:51:36    30.9    INFO                             Using 1 cores
## 2020-05-27 06:51:36    30.9  STATUS                         COMPLETED Setting up Multicore
## 2020-05-27 07:32:13    10.4  STATUS                     COMPLETED Compute methQTL per correlation block
## 2020-05-27 07:32:13    10.4  STATUS                 COMPLETED Computing methQTL for chromosome chr8
## 2020-05-27 07:32:13    10.4  STATUS                 STARTED Computing methQTL for chromosome chr9
## 2020-05-27 07:32:13    10.4  STATUS                     STARTED Compute correlation blocks
## 2020-05-27 07:32:13    10.4  STATUS                         STARTED Compute correlation matrix
## 2020-05-27 07:32:26    16.9  STATUS                         COMPLETED Compute correlation matrix
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
## Error in H5Dwrite(h5dataset, obj, h5spaceMem = h5spaceMem, h5spaceFile = h5spaceFile) : 
##   HDF5. Dataset. Write failed.
```

```
## Error in H5Fcreate(file): HDF5. File accessibilty. Unable to open file.
```

We will now present the two steps of the methQTL calling procedure in more detail.

## Compute CpG correlation blocks

Since neighboring CpGs are often highly correlated, using each CpG independently as a potential methQTL candidate leads to many redundant results. We thus aimed to approximate *DNA methylation haplotypes* by determining highly correlated CpGs in close vicinity. The procedure itself is split into six steps, and is performed for each chromosome independently:

1. Compute the (Pearson) correlation matrix between all CpGs
2. Construct the distance matrix from the correlation matrix
3. Discard all interactions with a correlation lower than a given threshold (option: ```cluster.cor.threshold```, default: 0.25)
4. Weight the distance according to the genomic distance between the two CpGs with a Gaussian (option: ```standard.deviation.gauss```, default: 100)
5. Discard all interactions ranging longer than the option ```absolute.distance.cutoff``` (default: 500,000)
6. Compute the Louvain clustering on the undirected, weighted graph induced by the distance matrix

This will return a clustering according to the correlation structure between neighboring CpGs that we will later use for methQTL calling.

## Call methQTL per correlation block

From the list of correlation blocks, *methQTL* computes methQTL interactions with all SNPs on the same chromosome. The process is split into three steps:

1. Compute a representative CpG per correlation block, as specified with the option ```representative.cpg.computation``` (default: *row.medians*).
2. Discard all SNPs that are further than ```absolute.distance.cutoff``` (default: 1,000,000) away from the representative CpG
3. Call methQTL by using linear models. Multiple options of methQTL calling are available and can be selected via the option ```linear.model.type``` (default: *classical.linear*). Alternatively, *fastQTL* can be set as an option for
```meth.qtl.type```. This will tell the package to use the fastQTL software[@Ongen2016].

The ```meth.qtl.type``` tells, how a methQTL interaction is defined and provides three options, in addition to the already mentioned *fastQTL*:

1. *oneVSall*: A CpG can only be influenced by one SNP. We choose the one with the lowest p-value.
2. *twoVSall*: A CpG can both positively and negatively be influenced by two independent SNPs. The package will output those fulfilling the p-value cutoff.
3. *allVSall*: For each CpG, all SNPs showing a p-value lower than the p-value cutoff will be returned.

In the latest stage, potential covariates can be specified using the option *sel.covariates*. We recommend to include at least *age* and *sex* as covariates, as they have a strong influence on the DNA methylation pattern.

# Downstream analysis and interpretation

## How to use *methQTLResult*

The above procedure will create an object of class ```methQTLResult```, which contains the methQTL that are called in the previous step. To get a table of all the methQTL, you need to extract the information from the object. In most of the function below, there is the option ```type```, which takes on the values:
* 'SNP': To characterize the SNPs that influence any DNA methylation state
* 'CpG': To characterize the representative CpGs per correlation block that are influences by any genotype
* 'cor.block': To characterize all CpGs, which are part of a correlation block, whose representative CpG is influenced by any genotype

Furthermore, you can obtain genomic annotations for both the CpGs and the SNPs involved in the methQTL interactions:


```r
result.table <- getResult(meth.qtl.res)
```

```
## Error in getResult(meth.qtl.res): object 'meth.qtl.res' not found
```

```r
head(result.table)
```

```
## Error in head(result.table): object 'result.table' not found
```

```r
anno.meth <- getAnno(meth.qtl.res,"meth")
```

```
## Error in getAnno(meth.qtl.res, "meth"): object 'meth.qtl.res' not found
```

```r
head(anno.meth)
```

```
## Error in head(anno.meth): object 'anno.meth' not found
```

```r
anno.geno <- getAnno(meth.qtl.res,"geno")
```

```
## Error in getAnno(meth.qtl.res, "geno"): object 'meth.qtl.res' not found
```

```r
head(anno.geno)
```

```
## Error in head(anno.geno): object 'anno.geno' not found
```

For more detailed information about the output, also see the function ```getResults.GWASMap```.

## Plots

To visualize methQTL, the package provides some plotting functions. All functions return an object of type ```ggplot```, which can be subsequently stored or viewed. Either all methQTL can be simultaneously visualized in a single plot, or a specific methQTL can be visualized:


```r
result.table <- result.table[order(result.table$P.value,decreasing=F),]
```

```
## Error in eval(expr, envir, enclos): object 'result.table' not found
```

```r
qtl.plot.SNP.CpG.interaction(meth.qtl.res,result.table$CpG[1],result.table$SNP[1])
```

```
## Error in getAnno(meth.qtl, "meth"): object 'meth.qtl.res' not found
```

```r
distance.scatterplot(meth.qtl.res)
```

```
## Error in distance.scatterplot(meth.qtl.res): could not find function "distance.scatterplot"
```

## Interpretation functions

The package provides a bunch of interpretation functions to characterize the detected methQTLs. This includes LOLA enrichment analysis[@LOLA] (```qtl.lola.enrichment```), genomic annotation enrichment based on putative regulatory elements defined by the Ensembl Regulatory Build[@Zerbino2015] (```qtl.annotation.enrichment```), enrichment analysis of different base substitutions in SNPs (```qtl.base.substitution.enrichment```), or TFBS motif enrichment using [TFBSTools](http://bioconductor.org/packages/release/bioc/html/TFBSTools.html). Enrichment is compared for the methQTLs that are available in the provided ```methQTLResult``` (for a single input), or to the overlapping QTLs for a list of ```methQTLResult``. The background of the enrichment is defined as all the SNPs/CpGs that have been used as input to the methQTL calling.


```r
res <- qtl.base.substitution.enricment(meth.qtl.res)
```

```
## Error in qtl.base.substitution.enricment(meth.qtl.res): could not find function "qtl.base.substitution.enricment"
```

```r
qtl.plot.base.substitution(meth.qtl.res,merge=TRUE)
```

```
## Error in get.overlap.universe(meth.qtl.res, type = "SNP"): object 'meth.qtl.res' not found
```

## Lists of methQTL results

Most of the functions discussed above either support a single ```methQTLResult``` as input, or a list of such objects. In case a list is specified, the functions with typically overlap the methQTLs found and compare those with all SNPs/CpGs that have been used for methQTL calling. Additionally, there are functions that particularly work on a list of ```methQTLResult``` objects and that perform overlapping, or determine the methQTLs specific to a dataset.


```r
meth.qtl.list <- list(First=meth.qtl.res.1,Second=meth.qtl.res.2,Third=meth.qtl.res.3)
qtl.venn.plot(meth.qtl.list)
qtl.upset.plot(meth.qtl.list,type = "cor.block")
spec.first <- get.specific.qtl(meth.qtl.list$First,meth.qtl.list[-1])
```

# Advanced configuration

## Employ methQTL on a scientific compute cluster

*methQTL* can automatically distribute jobs across a high performance compute cluster, which has been setup using the Sun Grid Engine (SGE) technology. You can pass the option ```cluster.submit``` to ```do.methQTL``` and thus activate the cluster submission. Note that you'll also have to specify a path to an executable *Rscript* and potentially specify resource requirements using the option setting ```cluster.config```.


```r
qtl.setOption(cluster.config = c(h_vmem="60G",mem_free="20G"))
qtl.setOption(rscript.path = "/usr/bin/Rscript")
meth.qtl.res <- do.methQTL(meth.qtl = imp.data,
                          cluster.submit = T)
```
