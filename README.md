GCAT (genotype conditional association test)
gcat
===

genotype conditional association test

Pre-print available: http://biorxiv.org/content/early/2015/03/04/012682

Published manuscript available:  http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3244.html

Dependencies
===

`gcat` requires the package `lfa` which is available at https://github.com/StoreyLab/lfa.  Apple OS X users with installation problems should read the FAQ on `lfa` page.

Installation
===

To install, open R and type:
```R
install.packages("devtools")
library("devtools")
install_github("Storeylab/lfa")
install_github("Storeylab/gcat")
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
