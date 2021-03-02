# direct deviance calculation, for tests only
deviance_snp <- function( xi, pi ) {
    if ( missing( xi ) )
        stop( 'Genotype vector `xi` is required!' )
    if ( missing( pi ) )
        stop( 'Individual-specific allele frequency vector `pi` is required!' )

    # for testing, it's easier to allow NAs here too, let's filter
    # assume xi can be NAs, and remaining pi shouldn't be NA! (don't remove those if present)
    if ( anyNA( xi ) ) {
        # indexes to keep
        indexes <- !is.na( xi )
        # remove cases that are NA in xi (from both input vectors)
        xi <- xi[ indexes ]
        pi <- pi[ indexes ]
    }

    # at this point data should have been filtered to remove NAs, let's make sure of that
    if ( anyNA( xi ) )
        stop( 'Genotype vector `xi` must not have NA values' )
    if ( anyNA( pi ) )
        stop( 'Individual-specific allele frequency vector `pi` must not have NA values' )
    
    # for a check below let's make sure xi are integers
    if ( storage.mode( xi ) != 'integer' )
        xi <- as.integer( xi )

    # these pi values are never NA since lfa 2.x
    # however, estimates can be at the "boundaries", particularly at 1, due to how `lfa::af_snp` solves one of the numerical issues
    # deviances in those cases, as long as the observations were also zero (x or 2-x), should to go zero instead of NA
    # we'll check to make sure this is what happened
    
    # handle the sum the first way
    # note: xi==0 are ignored in this sum
    indexes <- xi > 0
    # just subset data
    xi1 <- xi[ indexes ]
    pi1 <- pi[ indexes ]
    # only danger case is if any pi1 == 0, cause infinite logs
    if ( any( pi1 == 0 ) )
        stop( 'Observed impossible case `x > 1` and `p = 0`!' )
    # add to running sum!
    deviance <- 2 * sum( xi1 * log( xi1 / ( 2 * pi1 ) ) )

    # handle the other way
    # first reflect
    xi2 <- 2L - xi # keep integer
    pi2 <- 1 - pi
    # now remove zeroes in this direction (same as xi == 2L)
    indexes <- xi2 > 0
    # just subset data
    xi2 <- xi2[ indexes ]
    pi2 <- pi2[ indexes ]
    # only danger case is if any pi2 == 0, cause infinite logs
    if ( any( pi2 == 0 ) )
        stop( 'Observed impossible case `x < 2` and `p = 1`!' )
    # add to running sum!
    deviance <- deviance + 2 * sum( xi2 * log( xi2 / ( 2 * pi2 ) ) )

    # now this deviance must be non-NA
    if ( is.na(deviance) )
        stop( 'deviance was NA somehow!' )
    return( deviance )
}
