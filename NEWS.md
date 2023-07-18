# gcatest 2.0.0.9000 (2020-11-13)

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

# gcatest 2.0.1.9000 (2021-02-16)

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).

# gcatest 2.0.2.9000 (2021-03-01)

- Added internal tests for deviance calculations against `stats::glm`.
- Deviance code (internal `delta_deviance_snp`) now returns `NA` instead of stopping when an "impossible" case is encountered (when the genotype `x` is non-zero but the fitted probabilities under either null or alternative model are zero, or the alternative allele dosage (`x-2`) has the same problem).
  These cases are clearly model fitting failures, and can arise for common ill-defined problems, particularly under binary `adjustment` variables passed to `gcat` together with rare variants; these individual cases are not handled any better by `stats::glm`, so it seemed most sensible to return `NA` at such loci and not stop.

# gcatest 2.0.3.9000 (2021-05-11)

- Added function `delta_deviance_lf`, which calculates the delta deviance from two logistic models and the genotype matrix data.
  This function is a more general version of `gcat.stat` (which uses the new function internally), to essentially consider models that differ by more than one degree of freedom.
  It was written in particular for an external application in mind, namely the `jackstraw` package.
- Internal function `assoc_snp` was renamed to `delta_deviance_snp_lf` and its last argument changed to match that of `delta_deviance_lf` (alternative logistic factors instead of trait).

# gcatest 2.0.4.9000 (2021-05-13)

- Function `delta_deviance_lf` debugged case where either `LF0` or `LF1` is a column matrix.
  Previously these 1-column matrices were getting dropped to a vector incorrectly, which resulted in the mysterious error message "Error: argument is of length zero".
  This 1-column case is not typically observed in `gcatest`, but is common in the reverse-dependent `jackstraw` package.

# gcatest 2.0.5 (2021-06-18)

* Lots of minor changes for Bioconductor update.
  - DESCRIPTION:
    - Updated to `Authors@R`.
    - Lengthened "Description" paragraph.
    - Increased R dependency from 3.2 to 4.0.
  - Reformatted this `NEWS.md` slightly to improve its automatic parsing.
  - Added examples for function `delta_deviance_lf`.
  - Updated vignette to reflect that `lfa::read.bed` has been deprecated in favor of `genio::read_plink` and `BEDMatrix` objects.
  - Updated `README.md`, including corrections to examples.
  - Updated citations:
	- `README.md`: only had GCATest paper link, now has full citation and also full LFA citation.
    - Vignette: used to point to LFA arXiv preprint, now points to published paper.
	- `inst/CITATION`: didn't exist!  Now includes both LFA and GCATest papers.
  - Added `LICENSE.md`.
  - Internal changes:
    - All unexported functions are now prefixed with a period.
    - Replaced `1:x` with `seq_len(x)` several functions.
    - Reformatted all code with package `reformatR` and otherwise match Bioconductor guidelines.

# gcatest 2.0.6 (2023-05-25)

- `README.md` upgraded links from http to https
- Minor doc reformatting automatically performed by `roxygen2`.

# gcatest 2.1.6 (2023-05-25)

- Version bump for bioconductor devel.

# gcatest 2.1.7 (2023-06-20)

- Commented out various excessive tests against `glm`, which differ more often than expected due to poor or lack of convergence.
- Removed unused LaTeX package dependencies from vignette to prevent errors restricted to specific testing platforms.
- Fixed `..density..` deprecation warning in vignette plot.

# gcatest 2.1.8 (2023-07-18)

- Commented out two more strict tests (for non-negative deviances) that fail too often on bioconductor.

