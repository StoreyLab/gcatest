#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#' @title Genotype conditional association test
#' @description Performs the GCAT test for association between SNPs
#' and trait, and returns the p-values.
#' @inheritParams lfa::lfa
#' @param trait vector 
#' @param adjustment matrix of adjustment variables
#' @export
#' @useDynLib gcat
gcat <- function(X, LF, trait, adjustment=NULL){
    devdiff = gcat.stat(X, LF, trait, adjustment)
    
    return(pchisq(devdiff, 1, lower.tail=FALSE))
}

#' @describeIn gcat returns the association statistics instead of the 
#' p-value.
gcat.stat <- function(X, LF, trait, adjustment=NULL){
    if(length(trait) != ncol(X))
        stop("trait vector and genotype matrix columns must be same")
    
    #check LFs
    
    #check adjustment var
    if(!is.null(adjustment)){
        if(class(adjustment) != "matrix")
            stop("expecting adjustment to be matrix")
        if(nrow(adjustment) != nrow(LF))
            stop("adjustment and LF have differing number of rows")
        LF = cbind(LF, adjustment)
    }
    
    LF = rbind(LF, LF)
    LF = cbind(LF, trait) #don't forget LFs include intercept
    devdiff = .Call("assoc", LF, X, 10, 1E-6)
    return(devdiff)
}
