#' @describeIn gcat returns the association statistics instead of the 
#' p-value.
#' @examples
#' gcat_stat <- gcat.stat(sim_geno, LF, sim_trait)
#' @export
gcat.stat <- function( X, LF, trait, adjustment = NULL ){
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( LF ) )
        stop( "`LF` matrix is required!" )
    if ( missing( trait ) )
        stop( "`trait` is required!" )
    
    # check class
    is_BEDMatrix <- FALSE
    if ( "BEDMatrix" %in% class(X) ) {
        is_BEDMatrix <- TRUE
    } else if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )

    # get dimensions
    if ( is_BEDMatrix ) {
        m <- ncol(X)
        n <- nrow(X)
    } else {
        n <- ncol(X)
        # m not used in this case
    }

    # check dimensions
    if ( nrow(LF) != n )
        stop( 'Number of individuals in `X` must equal number of individuals (rows) in `LF`' )
    if ( length(trait) != n )
        stop( 'Number of individuals in `X` must equal number of individuals in `trait`' )
    
    # check trait and LF for missing values
    if ( anyNA( trait ) )
        stop( '`trait` must not have missing values!' )
    if ( anyNA( LF ) )
        stop( '`LF` must not have missing values!' )
    
    # check adjustment var
    if ( !is.null( adjustment ) ){
        if ( !is.matrix( adjustment ) )
            stop( '`adjustment` must be matrix!' )
        if ( nrow(adjustment) != n )
            stop( 'Number of individuals in `X` must equal number of individuals (rows) in `adjustment`' )
        # if all good, add adjustments to LFs
        LF <- cbind( LF, adjustment )
    }
    
    if ( is_BEDMatrix ) {
        # explicit loop for BEDMatrix
        devdiff <- vector('numeric', m)
        for ( i in 1 : m ) {
            # get locus i genotype vector
            xi <- X[ , i ]
            # calculate and store result
            devdiff[ i ] <- assoc_snp( xi, LF, trait )
        }
    } else
        devdiff <- apply( X, 1, assoc_snp, LF, trait )
    return( devdiff )
}

