#' Calculate delta deviance of logistic fits of a null and alternative model
#'
#' This function fits, at each locus of a given genotype matrix, two logistic models, and under the assumption that the models are nested, calculates the delta deviance between the two.
#' This general function is intended for testing models in a broad setting; for the specific problem of genetic association, the interface in [gcat()] and [gcat.stat()] are more user-friendly.
#'
#' @inheritParams lfa::lfa
#' @param LF0 Logistic factors for null model.
#' @param LF1 Logistic factors for alternative model.
#'
#' @return The vector of delta deviance values, one per locus of `X`.
#'
#' @examples
#' \dontrun{
#' devdiff <- delta_deviance_lf( X, LF0, LF1 )
#' }
#'
#' @export
delta_deviance_lf <- function( X, LF0, LF1 ){
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( LF0 ) )
        stop( "`LF0` matrix is required!" )
    if ( missing( LF1 ) )
        stop( "`LF1` matrix is required!" )
    
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
    if ( nrow(LF0) != n )
        stop( 'Number of individuals in `X` must equal number of individuals (rows) in `LF0`' )
    if ( nrow(LF1) != n )
        stop( 'Number of individuals in `X` must equal number of individuals (rows) in `LF1`' )
    
    # check LFs for missing values
    if ( anyNA( LF0 ) )
        stop( '`LF0` must not have missing values!' )
    if ( anyNA( LF1 ) )
        stop( '`LF1` must not have missing values!' )

    # start processing
    if ( is_BEDMatrix ) {
        # explicit loop for BEDMatrix
        devdiff <- vector('numeric', m)
        for ( i in 1 : m ) {
            # get locus i genotype vector
            xi <- X[ , i ]
            # calculate and store result
            devdiff[ i ] <- delta_deviance_snp_lf( xi, LF0, LF1 )
        }
    } else
        devdiff <- apply( X, 1, delta_deviance_snp_lf, LF0, LF1 )
    return( devdiff )
}

