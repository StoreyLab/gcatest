GCAT (genotype conditional association test)
gcatest
===

`gcatest` implements the genotype conditional association test (GCAT).

Pre-print available: http://biorxiv.org/content/early/2015/03/04/012682

Published manuscript available:  http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3244.html

Dependencies
===

`gcatest` requires the package `lfa` which is available at https://github.com/StoreyLab/lfa.  Apple OS X users with installation problems should read the FAQ on `lfa` page.

Installation
===

To install, open R and type:
```R
install.packages("devtools")
library("devtools")
install_github("Storeylab/lfa")
install_github("Storeylab/gcatest")
```

Example
===

`gcat` includes a simple example:

```R
library(gcat)
LF = lfa(sim_geno, 3)
gcat_p = gcat(sim_geno, LF, sim_trait)
```

The example is available in PLINK format at:

* http://genomics.princeton.edu/storeylab/data/gcat/demo/sim_geno.bed
* http://genomics.princeton.edu/storeylab/data/gcat/demo/sim_geno.bim
* http://genomics.princeton.edu/storeylab/data/gcat/demo/sim_geno.fam

The package `lfa` has the function `read.bed`. Example:

```R
library(gcat)
sim_geno = read.bed(bed.prefix="sim_geno")
sim_trait = read.table("sim_geno.fam")[,6]
LF = lfa(sim_geno, 3)
gcat_p = gcat(sim_geno, LF, sim_trait)
```

Checking genotype model fit
===

The main assumption that needs to verified on real data before using GCAT is that the probabilistic model of population structure fits the genotype data well.  Note that this verification does not involve the trait model, which is an important and positive aspect of GCAT.  The function `model.gof` returns a p-value for each SNP based on simulating a null distribution for the population structure model. The lower the p-value is for a given SNP, the worse the model fits that particular SNP.  Statistically significant p-values tell us which SNPs the model fails for, and those SNPs should be filtered out if necessary before using the GCAT test.  We can also adjust the value of `d` (which is the number of logistic factors included in the population structure model) to try to maximize the number of SNPs that are included in the GCAT analysis. In the example simulated data set, the last five SNPs are simulated to violate the model.

```R
set.seed(1234)
library(gcat)
sim_geno = read.bed(bed.prefix="sim_geno")
LF = lfa(sim_geno, 3)
gof = model.gof(sim_geno, LF, B=2)
filtered = gof < (1 / nrow(sim_geno))
sim_geno = sim_geno[!filtered,]

LF = lfa(sim_geno, 3)
gcat_p = gcat(sim_geno, LF, sim_trait)
```
