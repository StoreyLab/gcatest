#' Genotype Conditional Association TEST
#'
#' Performs the GCAT association test between SNPs and trait, returning
#' p-values.
#' 
#' @inheritParams lfa::lfa
#' @param LF matrix of logistic factors from [lfa::lfa()]
#' @param trait vector 
#' @param adjustment matrix of adjustment variables
#' @references Song, M, Hao, W, Storey, JD (2015). Testing for genetic
#' associations in arbitrarily structured populations. Nat. Genet., 47, 5:550-4.
#' @examples
#' library(lfa)
#' 
#' # make example data smaller so example is fast
#' # goes from 1000 to 100 individuals
#' indexes <- sample.int( ncol(sim_geno), 100 )
#' sim_geno <- sim_geno[ , indexes ]
#' sim_trait <- sim_trait[ indexes ]
#'
#' # now run LFA and GCATest
#' LF <- lfa(sim_geno, 3)
#' gcat_p <- gcat(sim_geno, LF, sim_trait)
#' @return vector of p-values 
#' @export
gcat <- function(X, LF, trait, adjustment = NULL) {
    devdiff <- gcat.stat(X, LF, trait, adjustment)

    return(stats::pchisq(devdiff, 1, lower.tail = FALSE))
}

#' @describeIn gcat Alias of gcat
#' @export
gcatest <- gcat
