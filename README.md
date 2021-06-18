# GCAT (genotype conditional association test)

`gcatest` implements the genotype conditional association test (GCAT).

## Installation

To install latest version on Bioconductor, open R and type:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("gcatest")
```

You can also install development version from GitHub this way:
```R
install.packages("devtools")
library("devtools")
install_github("Storeylab/gcatest")
```

## Example

`gcatest` includes a simple example:

```R
library(gcatest)
library(lfa)
LF <- lfa(sim_geno, 3)
gcat_p <- gcat(sim_geno, LF, sim_trait)
```

The example is also available in PLINK format at:

* http://genomics.princeton.edu/storeylab/data/gcat/demo/sim_geno.bed
* http://genomics.princeton.edu/storeylab/data/gcat/demo/sim_geno.bim
* http://genomics.princeton.edu/storeylab/data/gcat/demo/sim_geno.fam

The package `genio` has the function `read_plink` to read this data.
Example:

```R
library(gcatest)
library(lfa)
library(genio)
data <- read_plink("sim_geno")
sim_geno <- data$X
sim_trait <- data$fam$pheno
LF <- lfa(sim_geno, 3)
gcat_p <- gcat(sim_geno, LF, sim_trait)
```

## Checking genotype model fit

The main assumption that needs to verified on real data before using GCAT is that the probabilistic model of population structure fits the genotype data well.  Note that this verification does not involve the trait model, which is an important and positive aspect of GCAT.  The function `model.gof` returns a p-value for each SNP based on simulating a null distribution for the population structure model. The lower the p-value is for a given SNP, the worse the model fits that particular SNP.  Statistically significant p-values tell us which SNPs the model fails for, and those SNPs should be filtered out if necessary before using the GCAT test.  We can also adjust the value of `d` (which is the number of logistic factors included in the population structure model) to try to maximize the number of SNPs that are included in the GCAT analysis. In the example simulated data set, the last five SNPs are simulated to violate the model.

```R
library(gcatest)
library(lfa)
library(genio)
data <- read_plink("sim_geno")
sim_geno <- data$X
sim_trait <- data$fam$pheno
LF <- lfa(sim_geno, 3)
gof <- sHWE(sim_geno, LF, B=2)
filtered <- gof < (1 / nrow(sim_geno))
sim_geno <- sim_geno[!filtered,]

LF <- lfa(sim_geno, 3)
gcat_p <- gcat(sim_geno, LF, sim_trait)
```

## Citations

Song, Minsun, Wei Hao, and John D. Storey. "Testing for Genetic Associations in Arbitrarily Structured Populations." Nature Genetics 47, no. 5 (May 2015): 550-54. [doi:10.1038/ng.3244](https://doi.org/10.1038/ng.3244).

Hao, Wei, Minsun Song, and John D. Storey. "Probabilistic Models of Genetic Variation in Structured Populations Applied to Global Human Studies." Bioinformatics 32, no. 5 (March 1, 2016): 713–21. [doi:10.1093/bioinformatics/btv641](https://doi.org/10.1093/bioinformatics/btv641). [arXiv](http://arxiv.org/abs/1312.2041).

