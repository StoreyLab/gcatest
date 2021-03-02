# combined calculation of both deviances to maximize numerical accuracy (this is more stable than calculating them separately and then calculating the difference)
# handles some cases where model fit is poor by assuming that pi1 should fit the data better than pi0
delta_deviance_snp <- function( xi, pi0, pi1 ) {
    if ( missing( xi ) )
        stop( 'Genotype vector `xi` is required!' )
    if ( missing( pi0 ) )
        stop( 'Individual-specific allele frequency vector `pi0` is required!' )
    if ( missing( pi1 ) )
        stop( 'Individual-specific allele frequency vector `pi1` is required!' )

    # at this point genotypes should have been filtered to remove NAs, let's make sure of that
    if ( anyNA( xi ) )
        stop( 'Genotype vector `xi` must not have NA values' )

    # pi0 and pi1 may have NAs, those get handled below!
    ## if ( anyNA( pi0 ) )
    ##     stop( 'Individual-specific allele frequency vector `pi0` must not have NA values' )
    ## if ( anyNA( pi1 ) )
    ##     stop( 'Individual-specific allele frequency vector `pi1` must not have NA values' )
    
    # for a check below let's make sure xi are integers
    if ( storage.mode( xi ) != 'integer' )
        xi <- as.integer( xi )

    # want to calculate this:
    # dev0 <- sum( xi * log( pi0 ) + ( 2 - xi ) * log( 1 - pi0 ) )
    # dev1 <- sum( xi * log( pi1 ) + ( 2 - xi ) * log( 1 - pi1 ) )
    # devdiff <- - 2 * ( dev0 - dev1 )

    # rewrite to group more comparable terms, key to maximize numerical accuracy by avoiding catastrophic cancellations
    ## devdiff <- 2 * sum(
    ##     xi * log( pi1/pi0 ) + ( 2 - xi ) * log( ( 1 - pi1 ) / ( 1 - pi0 ) )
    ## )

    # try to understand NAs in pi0, pi1, and how to deal with them
    if ( anyNA( pi0 ) ) {
        indexes <- is.na( pi0 )
        if ( all( indexes ) ) {
            # LFA/glm.fit can just fail sometimes, have to accept that
            # won't be considered an error
            return( NA ) 
        } else {
            # this is harder to explain, haven't encountered yet but I'm ready if it ever happens:
            # report number and percentage of cases
            nNA <- length( indexes )
            pNA <- round( length( indexes ) / length( pi0 ) * 100, 3 )
            stop( 'pi0 had NAs: ', nNA, ', ', pNA , '%' )
        }
    }
    if ( anyNA( pi1 ) ) {
        indexes <- is.na( pi1 )
        if ( all( indexes ) ) {
            # LFA/glm.fit can just fail sometimes, have to accept that
            # won't be considered an error
            return( NA )
        } else {
            # report number and percentage of cases
            nNA <- length( indexes )
            pNA <- round( length( indexes ) / length( pi1 ) * 100, 3 )
            stop( 'pi1 had NAs: ', nNA, ', ', pNA , '%' )
        }
    }

    # let's handle other problematic edge cases I observed in toy simulations
    # sometimes LFA returns all zeroes (or all 1s?) for the allele frequencies
    # will return valid numbers (including Inf or -Inf), which result in non-NA p-values!
    # this way none of these are "errors"
    # look at zero cases
    p0_all_0 <- all( pi0 == 0 )
    p1_all_0 <- all( pi1 == 0 )
    # both equal treat as equal fit, zero deviance
    # otherwise one model was infinitely better than the other, will give p-values in correct extreme
    if ( p0_all_0 )
        return( if ( p1_all_0 ) 0 else Inf )
    # if we're here, this implies !p0_all_0
    if ( p1_all_0 )
        return( -Inf )
    # look at one cases
    p0_all_1 <- all( pi0 == 1 )
    p1_all_1 <- all( pi1 == 1 )
    if ( p0_all_1 )
        return( if ( p1_all_1 ) 0 else Inf )
    # if we're here, this implies !p0_all_1
    if ( p1_all_1 )
        return( -Inf )
    # don't expect to see p0_all_0 && p1_all_1, or other way around, won't try to handle those here
    # may end up getting calculated as NaN, meh
    
    # let's calculate sum in parts
    # add factor of 2 in the end

    # In all cases we exclude terms of the form 0 * log( 0 ), which in the limit equal 0 (lim_{p -> 0} p*log(p) = 0)
    # however, we do get pi estimates can be at the "boundaries", particularly at 1, due to how `lfa::af_snp` solves one of the numerical issues, so these must be excluded explicitly, otherwise R sets the whole thing to NA
    
    # first this sum of terms: xi * log( pi1/pi0 )
    # only xi != 0 contribute to sum, but first let's also check that the included values don't have any pi1 == 0 or pi0 == 0, if those two don't happen we're good!
    # NOTE: when those cases do occur, a bad fit is certainly to blame (probably an ill-defined problem, bad design matrix or other issues).  GLM sucks at these too.  Let's just not die and return NA
    indexes <- xi != 0L
    if ( any( pi0[indexes] == 0 ) )
        return ( NA )
    if ( any( pi1[indexes] == 0 ) )
        return ( NA )
    # what is left we are good to sum
    devdiff <- sum( xi[ indexes ] * log( pi1[ indexes ] / pi0[ indexes ] ) )
    
    # now let's do the second sum: ( 2 - xi ) * log( ( 1 - pi1 ) / ( 1 - pi0 )
    # the check is the same but "reflected"
    # only xi != 2 contribute to sum, but first let's also check that the included values don't have any pi1 == 1 or pi0 == 1, if those two don't happen we're good!
    indexes <- xi != 2L
    if ( any( pi0[indexes] == 1 ) )
        return ( NA )
    if ( any( pi1[indexes] == 1 ) )
        return ( NA )
    # what is left we are good to sum
    devdiff <- devdiff + sum( ( 2L - xi[ indexes ]) * log( ( 1 - pi1[ indexes ]) / ( 1 - pi0[ indexes ] ) ) )
    
    if ( is.na(devdiff) )
        stop( 'devdiff was NA somehow!' )
    
    # include the final factor of 2
    devdiff <- 2 * devdiff
    
    ## # theoretically, should have `devdiff >= 0`
    ## # in practice, limited numerical precision may produce small negative values
    ## # check that it's not too bad, otherwise force non-negativity now
    ## if ( devdiff < 0 ) {
    ##     if ( devdiff < - sqrt( .Machine$double.eps ) ) {
    ##         stop( 'Got a negative delta deviance: ', devdiff )
    ##     } else
    ##         devdiff <- 0
    ## }
    
    return( devdiff )
}
