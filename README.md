---
# IlluminaMouseMethylationmanifest

This R package was create for analysis of Illumina Mouse Methylation Array datasets using minfi. The annotation was downloaded from Illumina (accessed 12 Aug 2022) and follows the development of IlluminaMethylationEPIC packages. 
It can be used in combination with the [IlluminaMouseMethylationanno.12.v1.mm10](https://github.com/chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10) package, which provides additional information on the proble location.

As several CpGs have multiple probes, Illumina IDs concatenated from cgID and additinal information are used rather than cg numbers on their own.

## Installation

```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("chiaraherzog/IlluminaMouseMethylationmanifest")
```

## Requirements

R >= 4.1.2, uses minfi (>= 1.40.0)

## Use

The manifest can be applied after an RGset has been created with minfi:

```
library(IlluminaMouseMethylationmanifest)

RGset@annotation <-  c(array = "IlluminaMouseMethylation", annotation = "12.v1.mm10")
```

