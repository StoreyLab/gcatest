delta_deviance_snp_lf <- function( xi, LF0, LF1 ){
    if ( missing( xi ) )
        stop( 'Genotype vector `xi` is required!' )
    if ( missing( LF0 ) )
        stop( "`LF0` matrix is required!" )
    if ( missing( LF1 ) )
        stop( "`LF1` matrix is required!" )

    # check dimensions
    n <- length( xi )
    if ( nrow(LF0) != n )
        stop( 'Number of individuals in `xi` must equal number of individuals (rows) in `LF0`' )
    if ( nrow(LF1) != n )
        stop( 'Number of individuals in `xi` must equal number of individuals (rows) in `LF1`' )
    
    # check trait and LF for missing values
    if ( anyNA( LF0 ) )
        stop( '`LF0` must not have missing values!' )
    if ( anyNA( LF1 ) )
        stop( '`LF1` must not have missing values!' )
    
    # remove individuals whose genotypes are NA
    # though `lfa::af_snp` allows NAs, the resulting imputed allele frequencies do not contribute to the deviance (because the xi's are factors of the sum)!  So best to exclude them from the beginning.
    indexes_keep <- !is.na(xi)
    xi <- xi[indexes_keep]
    LF0 <- LF0[indexes_keep,]
    LF1 <- LF1[indexes_keep,]
    
    # perform two logistic regressions, under the null and alternative, resulting in estimated allele frequencies under each model
    p0 <- lfa::af_snp( xi, LF0 )
    p1 <- lfa::af_snp( xi, LF1 )
    
    # now compute delta deviance
    # uses a special numerically stable algorithm
    devdiff <- delta_deviance_snp( xi, p0, p1 )
    
    return( devdiff )
}
