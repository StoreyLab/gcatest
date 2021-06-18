# internal testing function, to compare our .delta_deviance_snp to. Based on
# code originally from jackstraw::dev.R
.deviance_snp_glm <- function(xi, LF) {
    # deviance as given by glm have to specify response using two columns
    # (cases vs failures), doesn't work to just pass `xi`
    x2 <- cbind(xi, 2 - xi)
    dev <- stats::glm(x2 ~ LF, family = stats::binomial())$deviance
    return(dev)
}

# another glm wrapper, this one for the most common delta deviance test present
# in our package (where trait is additional covariate)
.delta_deviance_snp_glm <- function(xi, LF, trait) {
    dev_glm1 <- .deviance_snp_glm(xi, cbind(LF, trait))
    dev_glm0 <- .deviance_snp_glm(xi, LF)
    devdiff_glm <- dev_glm0 - dev_glm1
    return(devdiff_glm)
}
