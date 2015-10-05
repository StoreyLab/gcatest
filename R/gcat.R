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
#' @description Performs the GCAT test for association between SNPs
#' and trait, and returns the p-values.
#' @inheritParams lfa::lfa
#' @param LF matrix of logistic factors outputed from function \code{lfa}
#' @param trait vector 
#' @param adjustment matrix of adjustment variables
#' @references Song, M, Hao, W, Storey, JD (2015). Testing for genetic associations in arbitrarily structured populations. Nat. Genet., 47, 5:550-4.
#' @examples
#' library(lfa)
#' LF = lfa(sim_geno, 3)
#' gcat_p = gcat(sim_geno, LF, sim_trait)
#' @return vector of p-values 
#' @export
#' @import lfa
#' @useDynLib gcatest
gcat <- function(X, LF, trait, adjustment=NULL){
    devdiff = gcat.stat(X, LF, trait, adjustment)
    
    return(pchisq(devdiff, 1, lower.tail=FALSE))
}

#' @describeIn gcat returns the association statistics instead of the 
#' p-value.
gcat.stat <- function(X, LF, trait, adjustment=NULL){
    if(length(trait) != ncol(X))
        stop("trait vector and genotype matrix columns must be same")
    
    #check trait for missing values
    if(sum(complete.cases(trait)) != length(trait))
        stop("NAs in trait")
    
    #check adjustment var
    if(!is.null(adjustment)){
        if(class(adjustment) != "matrix")
            stop("expecting adjustment to be matrix")
        if(nrow(adjustment) != nrow(LF))
            stop("adjustment and LF have differing number of rows")
        LF = cbind(LF, adjustment)
    }
    
    #no missing values
    if(sum(is.na(X)) == 0){
        LF = rbind(LF, LF)
        LF = cbind(LF, trait) #don't forget LFs include intercept
        devdiff = .Call("assoc", LF, X, 10, 1E-6)
        return(devdiff)
    } else{ #missing values
        devdiff = apply(X, 1, assoc_snp_na, LF, trait)
    }
    
    
}

assoc_snp_na = function(snp, LF, trait){
    ind = !is.na(snp) #logical index vector
    snp_no_na = snp[ind]
    LF_no_na = LF[ind,]
    #if(is.matrix(trait)){
    #    trait_no_na = trait[ind,]
    #} else{
        trait_no_na = trait[ind]
    #}
    
    b0 = lfa:::lreg(snp_no_na, LF_no_na)
    b1 = lfa:::lreg(snp_no_na, cbind(LF_no_na, trait_no_na))
    
    #est0 = .Call("mv", LF_no_na, b0)
    est0 = LF_no_na %*% b0
    #est1 = .Call("mv", cbind(LF_no_na, trait_no_na), b1)
    est1 = cbind(LF_no_na, trait_no_na) %*% b1
    
    p0 = exp(est0)/(1+exp(est0))
    p1 = exp(est1)/(1+exp(est1))
    
    devalt  = sum(snp_no_na*log(p1) + (2-snp_no_na)*log(1-p1))
    devnull = sum(snp_no_na*log(p0) + (2-snp_no_na)*log(1-p0))
    
    -2*(devnull-devalt)
}
