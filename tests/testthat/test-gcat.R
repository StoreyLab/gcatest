library(lfa)

# generate random data for tests

# data dimensions
n_ind <- 10
m_loci <- 300
# total data size
n_data <- n_ind * m_loci
# add missingness
miss <- 0.1

# completely unstructured genotypes
# create ancestral allele frequencies
p_anc <- runif( m_loci )
# create genotypes
X <- rbinom( n_data, 2, p_anc )
# add missing values
X[ sample( n_data, n_data * miss ) ] <- NA
# turn into matrix
X <- matrix( X, nrow = m_loci, ncol = n_ind )

# to have a reasonable dataset always, remove fixed loci and all NA loci
# first remove loci that are entirely NA (with just 10 indiviuals, very possible)
loci_keep <- rowSums( !is.na(X) ) > 0
X <- X[ loci_keep, ]
# now identify fixed loci
p_anc_hat <- rowMeans( X, na.rm = TRUE )
loci_keep <- (0 < p_anc_hat) & (p_anc_hat < 1)
X <- X[ loci_keep, ]
# update number of loci and data size
m_loci <- nrow( X )
n_data <- n_ind * m_loci

# get LFs from this data
d <- 3
LFs <- lfa( X, d )

# make adjustments for that special case
adjustment <- cbind( rbinom( n_ind, 1, 0.5) )

# simulate a completely random trait
trait <- rnorm( n_ind )

# I love my `glm` comparisons, but note these are often not "silent"
# # Warning message:
# # glm.fit: fitted probabilities numerically 0 or 1 occurred
# we must use `suppressWarnings`, otherwise these warnings show up as testing errors

# this function is not used in practice, but is kept for accuracy checks
test_that( "deviance_snp works", {
    # first draw completely random data
    # OLD: do not use X because we can't have missingness here.
    # NEW: missingness is handled by all functions!!!
    pi <- runif( n_ind )
    xi <- rbinom( n_ind, 2, pi)
    expect_silent(
        dev <- deviance_snp( xi, pi )
    )
    expect_equal( length( dev ), 1 )
    expect_true( is.numeric( dev ) )
    expect_true( !is.na( dev ) )
    expect_true( dev >= 0 )

    # repeat test where there are exact zeroes or ones
    pi[1] <- 0
    pi[n_ind] <- 1
    xi <- rbinom( n_ind, 2, pi)
    expect_silent(
        dev <- deviance_snp( xi, pi )
    )
    expect_equal( length( dev ), 1 )
    expect_true( is.numeric( dev ) )
    expect_true( !is.na( dev ) )
    expect_true( dev >= 0 )

    # direct comparison to `glm`
    suppressWarnings(
        dev_glm <- deviance_snp_glm(xi, LFs)
    )
    # replicate fit with our code, needs two stages (pi's must be fit from LFs)
    pi_fit <- lfa::af_snp(
                       xi,
                       LFs
                   )
    dev <- deviance_snp( xi, pi_fit )
    expect_equal( dev, dev_glm )

    # now test all loci, for complete assurance
    dev <- vector( 'numeric', m_loci )
    dev_glm <- vector( 'numeric', m_loci )
    for ( i in 1 : m_loci ) {
        # extract row
        xi <- X[ i, ]
        # direct comparison to `glm`
        suppressWarnings(
            dev_glm[i] <- deviance_snp_glm( xi, LFs )
        )
        # replicate fit with our code, needs two stages (pi's must be fit from LFs)
        pi_fit <- lfa::af_snp( xi, LFs )
        dev[i] <- deviance_snp( xi, pi_fit )
    }
    expect_equal( dev, dev_glm )
})

test_that( "delta_deviance_snp works", {
    # first draw completely random data
    # do not use X because we can't have missingness here
    # true allele frequencies
    pi <- runif( n_ind )
    # observed genotypes
    xi <- rbinom( n_ind, 2, pi)
    # estimate allele frequencies them in the same way as GCAT does it, through LFA, to ensure delta deviance fits as expected (that alternative fits better than null in an absolute sense, not necessarily significantly so though)
    p0 <- lfa::af_snp( xi, LFs )
    p1 <- lfa::af_snp( xi, cbind(LFs, trait) )

    expect_silent(
        devdiff <- delta_deviance_snp( xi, p0, p1 )
    )
    expect_equal( length( devdiff ), 1 )
    expect_true( is.numeric( devdiff ) )
    expect_true( !is.na( devdiff ) )
    expect_true( devdiff >= 0 ) # can fail because models are not always fit well by LFA
    # also compare to less numerically-stable but otherwise identical calculation
    devdiff2 <- deviance_snp( xi, p0 ) - deviance_snp( xi, p1 )
    expect_equal( devdiff, devdiff2 )

    # direct comparison to `glm`
    suppressWarnings(
        devdiff_glm <- delta_deviance_snp_glm( xi, LFs, trait )
    )
    expect_equal( devdiff, devdiff_glm )
    
    # repeat test where there are exact zeroes or ones in the true data
    # (doesn't ensure that for LFA, but biases the estimates certainly)
    pi[1] <- 0
    pi[n_ind] <- 1
    xi <- rbinom( n_ind, 2, pi)
    p0 <- lfa::af_snp( xi, LFs )
    p1 <- lfa::af_snp( xi, cbind(LFs, trait) )
    # manually insert trouble cases, one in each vector, in the same place as before since that way we know the genotypes are not impossible
    p0[1] <- 0
    p1[n_ind] <- 1
    expect_silent(
        devdiff <- delta_deviance_snp( xi, p0, p1 )
    )
    expect_equal( length( devdiff ), 1 )
    expect_true( is.numeric( devdiff ) )
    expect_true( !is.na( devdiff ) )
    # expect_true( devdiff >= 0 ) # can fail because models are not always fit well by LFA (seems especially bad in these artificially bad cases)
    # also compare to less numerically-stable but otherwise identical calculation
    devdiff2 <- deviance_snp( xi, p0 ) - deviance_snp( xi, p1 )
    expect_equal( devdiff, devdiff2 )

    # NOTE: because these trouble cases were not obtained from direct fit from LFs (they were hacked into the resulting pi estimates), then no direct `glm` comparison is possible here
    
    # test impossible cases, make sure there are the expected NA return values
    xi <- 0:2
    pi <- xi / 2 # perfect allele frequencies
    # impossible cases
    pi_bad1 <- pi
    pi_bad1[1] <- 1 # can't have this if data was xi[1] == 0
    pi_bad2 <- pi
    pi_bad2[3] <- 0 # can't have this if data was xi[1] == 2
    # test in all combinations
    expect_true( is.na( delta_deviance_snp( xi, pi_bad1, pi ) ) )
    expect_true( is.na( delta_deviance_snp( xi, pi_bad2, pi ) ) )
    expect_true( is.na( delta_deviance_snp( xi, pi, pi_bad1 ) ) )
    expect_true( is.na( delta_deviance_snp( xi, pi, pi_bad2 ) ) )
})

test_that("delta_deviance_snp_lf works", {
    # test a single SNP
    i <- 1
    LF1 <- cbind(LFs, trait)
    expect_silent(
        devdiff <- delta_deviance_snp_lf( X[ i, ] , LFs, LF1 )
    )
    expect_equal( length( devdiff ), 1 )
    expect_true( is.numeric( devdiff ) )
    # delta deviances can be NA if LFA/glm.fit fail to converge
    expect_true( !is.na( devdiff ) )
    expect_true( devdiff >= 0 )

    # direct comparison to `glm`
    suppressWarnings(
        devdiff_glm <- delta_deviance_snp_glm( X[ i, ], LFs, trait )
    )
    expect_equal( devdiff, devdiff_glm )
})

test_that("delta_deviance_lf works", {
    LF1 <- cbind(LFs, trait)
    expect_silent(
        devdiff <- delta_deviance_lf( X , LFs, LF1 )
    )
    expect_equal( length( devdiff ), m_loci )
    expect_true( is.numeric( devdiff ) )
    # delta deviances can be NA if LFA/glm.fit fail to converge
    expect_true( !any( is.na( devdiff ) ) )
    # theoretically this is true, but in practice it depends on LFA fitting these models well, so testing this is not appropriate here (fails sometimes)
    ## expect_true( all( devdiff >= 0 ) )
    
    # loop through `glm` version
    devdiff_glm <- vector( 'numeric', m_loci )
    for ( i in 1 : m_loci ) {
        suppressWarnings(
            devdiff_glm[i] <- delta_deviance_snp_glm( X[ i, ], LFs, trait )
        )
    }
    expect_equal( devdiff, devdiff_glm )
})

test_that("gcat.stat works", {
    expect_silent(
        devdiff <- gcat.stat( X, LFs, trait )
    )
    expect_equal( length( devdiff ), m_loci )
    expect_true( is.numeric( devdiff ) )
    # delta deviances can be NA if LFA/glm.fit fail to converge
    expect_true( !anyNA( devdiff ) )
    # theoretically this is true, but in practice it depends on LFA fitting these models well, so testing this is not appropriate here (fails sometimes)
    ## expect_true( all( devdiff >= 0 ) )

    # loop through `glm` version
    devdiff_glm <- vector( 'numeric', m_loci )
    for ( i in 1 : m_loci ) {
        suppressWarnings(
            devdiff_glm[i] <- delta_deviance_snp_glm( X[ i, ], LFs, trait )
        )
    }
    expect_equal( devdiff, devdiff_glm )
    
    # test version with adjustments
    expect_silent(
        devdiff <- gcat.stat( X, LFs, trait, adjustment )
    )
    expect_equal( length( devdiff ), m_loci )
    expect_true( is.numeric( devdiff ) )
    # NOTE: the binary adjustments we use here cause bad fits that even `glm` doesn't handle well (it appears to cause ill-defined problems sometimes), so we don't bother testing against `glm` here
})

test_that( "gcat works", {
    expect_silent(
        pvals <- gcat(X, LFs, trait)
    )
    expect_equal( length(pvals), m_loci )
    expect_true( is.numeric(pvals) )
    expect_true( all(pvals >= 0, na.rm = TRUE) )
    expect_true( all(pvals <= 1, na.rm = TRUE) )

    # repeat with adjustments
    expect_silent(
        pvals <- gcat(X, LFs, trait, adjustment = adjustment)
    )
    expect_equal( length(pvals), m_loci )
    expect_true( is.numeric(pvals) )
    expect_true( all(pvals >= 0, na.rm = TRUE) )
    expect_true( all(pvals <= 1, na.rm = TRUE) )
})

### BEDMatrix tests

# require external packages for this...

if (
    suppressMessages(suppressWarnings(require(BEDMatrix))) &&
    suppressMessages(suppressWarnings(require(genio)))
) {
    context('gcat_BEDMatrix')
    
    # write the same data we simulated onto a temporary file
    file_bed <- tempfile('delete-me-random-test') # output name without extensions!
    genio::write_plink( file_bed, X )

    # load as a BEDMatrix object
    X_BEDMatrix <- suppressMessages(suppressWarnings( BEDMatrix( file_bed ) ))

    test_that("delta_deviance_lf works with BEDMatrix", {
        LF1 <- cbind( LFs, trait )
        devdiff_basic <- delta_deviance_lf( X, LFs, LF1 )
        expect_silent(
            devdiff_BM <- delta_deviance_lf( X_BEDMatrix, LFs, LF1 )
        )
        expect_equal( devdiff_basic, devdiff_BM )
    })

    test_that("gcat.stat works with BEDMatrix", {
        devdiff_basic <- gcat.stat( X, LFs, trait )
        expect_silent(
            devdiff_BM <- gcat.stat( X_BEDMatrix, LFs, trait )
        )
        expect_equal( devdiff_basic, devdiff_BM )
        
        # test version with adjustments
        devdiff_basic <- gcat.stat( X, LFs, trait, adjustment )
        expect_silent(
            devdiff_BM <- gcat.stat( X_BEDMatrix, LFs, trait, adjustment )
        )
        expect_equal( devdiff_basic, devdiff_BM )
    })

    test_that( "gcat works with BEDMatrix", {
        pvals_basic <- gcat(X, LFs, trait)
        expect_silent(
            pvals_BM <- gcat(X_BEDMatrix, LFs, trait)
        )
        expect_equal( pvals_basic, pvals_BM )

        # repeat with adjustments
        pvals_basic <- gcat(X, LFs, trait, adjustment = adjustment)
        expect_silent(
            pvals_BM <- gcat(X_BEDMatrix, LFs, trait, adjustment = adjustment)
        )
        expect_equal( pvals_basic, pvals_BM )
    })

    # delete temporary data when done
    genio::delete_files_plink( file_bed )
}
