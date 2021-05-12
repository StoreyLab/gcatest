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
    
    # check class, get dimensions
    if ( "BEDMatrix" %in% class(X) ) {
        n <- nrow(X)
    } else if ( is.matrix( X ) ) {
        n <- ncol(X)
    } else 
        stop( '`X` must be a matrix!' )

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

    # bind matrices
    LF1 <- cbind(LF, trait)

    # this performs the calculations in a broader setting
    # (shared with `jackstraw` package)
    devdiff <- delta_deviance_lf( X, LF, LF1 )
    
    return( devdiff )
}

