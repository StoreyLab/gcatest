# 2020-11-13 - gcatest 2.0.0.9000

Major overhaul from last version (1.3.2, last updated 2016-10-06).
Visible differences are support for BEDMatrix and fewer cases in which association p-values are `NA`.
Internally there was major code restructuring, and added unit tests for all functions.

- User-facing changes: Functions `gcat`/`gcatest`/`gcat.stat`
  - added support for BEDMatrix objects for the genotype matrix `X`.
    - This consumes lower memory when the number of loci `m` is very large, so it enables analysis of larger datasets.
  - Fixed some cases where the test statistic (the delta deviance) and ultimately the p-values were `NA` or `NaN` and are no longer missing.
    - One common case is when fitted probabilities were zero or one, which used to lead to `NaN` deviances when their correct contribution was instead zero (because the limit of `p*log(p)` as `p` goes to zero is zero, not `0 * (-Inf) = NaN`).
	- Other `NA` and `NaN` cases are avoided in the `lfa` function `af_snp` (fixed in lfa 2.0.0.9000, 2020-09-18) used to estimate the individual-specific allele frequencies used here to compute the delta deviance.
	  However, in rare cases the logistic regression in `af_snp` fails to converge or there are other problems, resulting in `NA` values propagated to GCATest's test statistic and p-values.
	- Otherwise, the new delta deviance code (function `delta_deviance_snp`) is more numerically-stable than before.

- Internal changes
  - Separated R functions into one source file each.
  - Added more input checks to all functions. 
  - Added `.gitignore` files from another project.
  - Added unit tests for all functions using `testthat`.
  - Removed internal `assoc` C code
    - Previously only used for genotype data without missingness (so practically not on real datasets)
    - Was entirely redundant with `lfa::af_snp`, which is now called in all cases instead.
	- Had bugs concerning handling of p == 0 or 1 cases that are better handled in `assoc_snp` R code
  - Minor scattered changes solely to pass latest `R CMD check` requirements.

# 2021-02-16 - gcatest 2.0.1.9000

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).

# 2021-03-01 - gcatest 2.0.2.9000

- Added internal tests for deviance calculations against `stats::glm`.
- Deviance code (internal `delta_deviance_snp`) now returns `NA` instead of stopping when an "impossible" case is encountered (when the genotype `x` is non-zero but the fitted probabilities under either null or alternative model are zero, or the alternative allele dosage (`x-2`) has the same problem).
  These cases are clearly model fitting failures, and can arise for common ill-defined problems, particularly under binary `adjustment` variables passed to `gcat` together with rare variants; these individual cases are not handled any better by `stats:glm`, so it seemed most sensible to return `NA` at such loci and not stop.
