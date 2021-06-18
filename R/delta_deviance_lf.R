#' Calculate delta deviance of logistic null/alternative models
#'
#' This function fits, at each locus of a given genotype matrix, two logistic
#' models, and under the assumption that the models are nested, calculates the
#' delta deviance between the two.
#' This general function is intended for testing models in a broad setting; for
#' the specific problem of genetic association, the interface in
#' [gcat()] and 
#' [gcat.stat()] are more user-friendly.
#'
#' @inheritParams lfa::lfa
#' @param LF0 Logistic factors for null model.
#' @param LF1 Logistic factors for alternative model.
#'
#' @return The vector of delta deviance values, one per locus of `X`.
#'
#' @examples
#' library(lfa)
#' 
#' # make example data smaller so example is fast
#' # goes from 1000 to 100 individuals
#' indexes <- sample.int( ncol(sim_geno), 100 )
#' sim_geno <- sim_geno[ , indexes ]
#' sim_trait <- sim_trait[ indexes ]
#'
#' # now run LFA and get delta deviances for trait assoc
#' # (recapitulating `gcat.stat` in this case)
#' LF <- lfa(sim_geno, 3)
#' LF0 <- LF # structure is null
#' LF1 <- cbind(LF, sim_trait) # trait is alt
#' devdiff_assoc <- delta_deviance_lf(sim_geno, LF0, LF1)
#'
#' # can instead do delta deviances for structure only
#' LF0 <- cbind(rep.int(1, ncol(sim_geno))) # intercept only is null
#' LF1 <- LF # structure is alt, no trait
#' devdiff_struc <- delta_deviance_lf(sim_geno, LF0, LF1)
#'
#' @export
delta_deviance_lf <- function(X, LF0, LF1) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(LF0))
        stop("`LF0` matrix is required!")
    if (missing(LF1))
        stop("`LF1` matrix is required!")
    # nothing can be null (`jackstraw`'s default LF0 is NULL i.e. intercept
    # only, disagreeing with `gcatest`, best to be explicit)
    if (is.null(X))
        stop("`X` cannot be null!")
    if (is.null(LF0))
        stop("`LF0` cannot be null!")
    if (is.null(LF1))
        stop("`LF1` cannot be null!")
    if (!is.matrix(X)) # check class; BEDMatrix returns TRUE
        stop("`X` must be a matrix!")
    if (methods::is(X, "BEDMatrix")) { # get dimensions
        m <- ncol(X)
        n <- nrow(X)
    } else n <- ncol(X) # m not used in this case
    if (nrow(LF0) != n) # check dimensions
        stop("Number of individuals in `X` and `LF0` disagrees!")
    if (nrow(LF1) != n)
        stop("Number of individuals in `X` and `LF1` disagrees!")
    if (anyNA(LF0)) # check LFs for missing values
        stop("`LF0` must not have missing values!")
    if (anyNA(LF1))
        stop("`LF1` must not have missing values!")
    # start actual processing
    if (methods::is(X, "BEDMatrix")) {
        devdiff <- vector("numeric", m)
        for (i in seq_len(m))  # explicit loop for BEDMatrix
            devdiff[i] <- .delta_deviance_snp_lf(X[, i], LF0, LF1)
    } else devdiff <- apply(X, 1, .delta_deviance_snp_lf, LF0, LF1)
    return(devdiff)
}

