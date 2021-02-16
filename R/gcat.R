#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#' @title Genotype Conditional Association TEST
#' @description Performs the GCAT test for association between SNPs and trait, and returns the p-values.
#' @inheritParams lfa::lfa
#' @param LF matrix of logistic factors outputed from function [lfa::lfa()]
#' @param trait vector 
#' @param adjustment matrix of adjustment variables
#' @references Song, M, Hao, W, Storey, JD (2015). Testing for genetic associations in arbitrarily structured populations. Nat. Genet., 47, 5:550-4.
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
#' @import lfa
gcat <- function(X, LF, trait, adjustment = NULL){
    devdiff <- gcat.stat(X, LF, trait, adjustment)
    
    return( stats::pchisq( devdiff, 1, lower.tail = FALSE ) )
}

#' @describeIn gcat Alias of gcat
#' @export
gcatest <- gcat
