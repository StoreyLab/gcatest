delta_deviance_snp_lf <- function( xi, LF0, LF1 ){
    if ( missing( xi ) )
        stop( 'Genotype vector `xi` is required!' )
    if ( missing( LF0 ) )
        stop( "`LF0` matrix is required!" )
    if ( missing( LF1 ) )
        stop( "`LF1` matrix is required!" )

    # nothing can be null either
    # (catches problems in `jackstraw`, where the default LF0 is NULL, but its defaults and `gcatest`'s are different.  It's better not to assume any particular defaults.)
    if ( is.null( xi ) )
        stop( '`xi` cannot be null!' )
    if ( is.null( LF0 ) )
        stop( '`LF0` cannot be null!' )
    if ( is.null( LF1 ) )
        stop( '`LF1` cannot be null!' )
    
    # check dimensions
    n <- length( xi )
    if ( nrow(LF0) != n )
        stop( 'Number of individuals in `xi` must equal number of individuals (rows) in `LF0`' )
    if ( nrow(LF1) != n )
        stop( 'Number of individuals in `xi` must equal number of individuals (rows) in `LF1`' )
    
    # check LFs for missing values
    if ( anyNA( LF0 ) )
        stop( '`LF0` must not have missing values!' )
    if ( anyNA( LF1 ) )
        stop( '`LF1` must not have missing values!' )
    
    # remove individuals whose genotypes are NA
    # though `lfa::af_snp` allows NAs, the resulting imputed allele frequencies do not contribute to the deviance (because the xi's are factors of the sum)!  So best to exclude them from the beginning.
    indexes_keep <- !is.na(xi)
    # TROUBLESHOOTING: though we could just return NA in these cases, they're not expected in any reasonable data
    if ( !any( indexes_keep ) )
        stop( 'All individuals at one locus were missing (unusual)!' )
    # save time if nothing gets removed
    if ( any( !indexes_keep ) ) {
        xi <- xi[ indexes_keep ]
        LF0 <- LF0[ indexes_keep, , drop = FALSE ]
        LF1 <- LF1[ indexes_keep, , drop = FALSE ]
    }
    
    # perform two logistic regressions, under the null and alternative, resulting in estimated allele frequencies under each model
    p0 <- lfa::af_snp( xi, LF0 )
    p1 <- lfa::af_snp( xi, LF1 )
    
    # now compute delta deviance
    # uses a special numerically stable algorithm
    devdiff <- delta_deviance_snp( xi, p0, p1 )
    
    return( devdiff )
}
