GCAT (genotype conditional association test)
gcat
===

genotype conditional association test

Pre-print available soon.

Dependencies
===

`gcat` requires the package `lfa` which is available at https://github.com/StoreyLab/lfa

Example
===

`gcat` includes a simple example:

```R
library(gcat)
LF = lfa(sim_geno, 3)
gcat_p = gcat(sim_geno, LF, sim_trait)
```
