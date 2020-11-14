# direct deviance calculation, for tests only
deviance_snp <- function( xi, pi ) {
    if ( missing( xi ) )
        stop( 'Genotype vector `xi` is required!' )
    if ( missing( pi ) )
        stop( 'Individual-specific allele frequency vector `pi` is required!' )

    # at this point data should have been filtered to remove NAs, let's make sure of that
    if ( anyNA( xi ) )
        stop( 'Genotype vector `xi` must not have NA values' )
    if ( anyNA( pi ) )
        stop( 'Individual-specific allele frequency vector `pi` must not have NA values' )
    
    # for a check below let's make sure xi are integers
    if ( storage.mode( xi ) != 'integer' )
        xi <- as.integer( xi )
    
    # set a small non-zero value to use to fix cases (exact value doesn't matter at all)
    eps <- .Machine$double.eps
    
    # these pi values are never NA since lfa 2.x
    # however, estimates can be at the "boundaries", particularly at 1, due to how `lfa::af_snp` solves one of the numerical issues
    # deviances in those cases, as long as the observations were also zero (x or 2-x), should to go zero instead of NA
    # let's check to make sure this is what happened
    # first look at zeroes
    indexes <- pi == 0
    if ( any( indexes ) ) {
        if ( all( xi[ indexes ] == 0L ) ) {
            # this is what we expected, just set these pi values to some epsilons and it all works out
            pi[ indexes ] <- eps
        } else {
            # else there was an unexpected error!
            stop( 'Observed impossible case `p = 0` and `x != 0`!' )
        }
    }

    # repeat with ones
    indexes <- pi == 1
    if ( any( indexes ) ) {
        if ( all( xi[ indexes ] == 2L ) ) {
            # this is what we expected, just set these pi values to some epsilons and it all works out
            pi[ indexes ] <- 1 - eps
        } else {
            # else there was an unexpected error!
            stop( 'Observed impossible case `p = 1` and `x != 2`!' )
        }
    }

    # now this deviance is non-NA
    deviance <- sum( xi * log( pi ) + ( 2 - xi ) * log( 1 - pi ) )
    if ( is.na(deviance) )
        stop( 'deviance was NA somehow!' )
    return( deviance )
}
