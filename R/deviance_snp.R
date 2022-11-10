# direct deviance calculation, for tests only
.deviance_snp <- function(xi, pi) {
    if (missing(xi))
        stop("Genotype vector `xi` is required!")
    if (missing(pi))
        stop("Individual-specific allele frequencies `pi` is required!")
    # for testing, allow NAs in `xi`, but not in `pi`
    if (anyNA(xi)) {
        indexes <- !is.na(xi)
        xi <- xi[indexes]
        pi <- pi[indexes]
    }
    if (anyNA(pi))
        stop("Individual-specific allele frequencies `pi` must not have NAs")
    # for a check below let's make sure xi are integers
    if (storage.mode(xi) != "integer")
        xi <- as.integer(xi)
    # handle the sum the first way. Note xi==0 are ignored in this sum
    indexes <- xi > 0
    xi1 <- xi[indexes]
    pi1 <- pi[indexes]
    # only danger case is if any pi1 == 0, cause infinite logs
    if (any(pi1 == 0))
        stop("Observed impossible case `x > 1` and `p = 0`!")
    # add to running sum!
    deviance <- 2 * sum(xi1 * log(xi1/(2 * pi1)))
    # handle the other way first reflect
    xi2 <- 2L - xi  # keep integer
    pi2 <- 1 - pi
    # now remove zeroes in this direction (same as xi == 2L)
    indexes <- xi2 > 0
    xi2 <- xi2[indexes]
    pi2 <- pi2[indexes]
    # only danger case is if any pi2 == 0, cause infinite logs
    if (any(pi2 == 0))
        stop("Observed impossible case `x < 2` and `p = 1`!")
    # add to running sum!
    deviance <- deviance + 2 * sum(xi2 * log(xi2/(2 * pi2)))
    # now this deviance must be non-NA
    if (is.na(deviance))
        stop("deviance was NA somehow!")
    return(deviance)
}
