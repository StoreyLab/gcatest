assoc_snp <- function(xi, LF, trait){
    if ( missing( xi ) )
        stop( 'Genotype vector `xi` is required!' )
    if ( missing( LF ) )
        stop( "`LF` matrix is required!" )
    if ( missing( trait ) )
        stop( "`trait` is required!" )

    # check dimensions
    n <- length( xi )
    if ( nrow(LF) != n )
        stop( 'Number of individuals in `xi` must equal number of individuals (rows) in `LF`' )
    if ( length(trait) != n )
        stop( 'Number of individuals in `xi` must equal number of individuals in `trait`' )
    
    # check trait and LF for missing values
    if ( anyNA( trait ) )
        stop( '`trait` must not have missing values!' )
    if ( anyNA( LF ) )
        stop( '`LF` must not have missing values!' )
    
    # remove individuals whose genotypes are NA
    # though `lfa::af_snp` allows NAs, the resulting imputed allele frequencies do not contribute to the deviance (because the xi's are factors of the sum)!  So best to exclude them from the beginning.
    indexes_keep <- !is.na(xi)
    xi <- xi[indexes_keep]
    LF <- LF[indexes_keep,]
    trait <- trait[indexes_keep]

    # perform two logistic regressions, under the null and alternative (add trait to LF/covariates), resulting in estimated allele frequencies under each model
    p0 <- lfa::af_snp( xi, LF )
    p1 <- lfa::af_snp( xi, cbind(LF, trait) )

    # now compute delta deviance
    # uses a special numerically stable algorithm
    devdiff <- delta_deviance_snp( xi, p0, p1 )
    
    return( devdiff )
}
