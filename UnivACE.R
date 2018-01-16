#---------------- UNIVARIATE TWIN MODELING IN lavaan --------------

library(lavaan)
library(MASS)

####### Data Simulation

# simulate data for monozygotic and dizygotic twins

n.MZ <- 250 # number of MZ twins
n.DZ <- 250 # number of DZ twins

mzcor <- .7
dzcor <- .4

S.MZ <- matrix( c( 1, mzcor, mzcor, 1 ), 2 ) # covariance matrix MZ twins
S.DZ <- matrix( c( 1, dzcor, dzcor, 1 ), 2 ) # covariance matrix DZ twins

MZdata <- mvrnorm( n.MZ, rep( 0, 2 ), S.MZ, emp = TRUE ) # generate MZ twins
DZdata <- mvrnorm( n.DZ, rep( 0, 2 ), S.DZ, emp = TRUE ) # generate DZ twins

# tie MZ and DZ data together in a dataframe
# with zygosity as grouping variable

twindata <- as.data.frame( rbind( cbind( MZdata, 1 ), cbind( DZdata, 2 ) ) )
colnames(twindata) <- c( "xt1", "xt2", "zyg" ) 

######## Model( ACE )

ACEmodel <- '
            # twin 1
            At1 =~ equal( c( h, h ) ) * xt1
            Ct1 =~ equal( c( c, c ) ) * xt1
            Et1 =~ equal( c( e, e ) ) * xt1

            # twin 2
            At2 =~ equal( c( h, h ) ) * xt2
            Ct2 =~ equal( c( c, c ) ) * xt2
            Et2 =~ equal( c( e, e ) ) * xt2

            # genetic and environmental similarity
            At1 ~~ c( 1,.5 ) * At2
            Ct1 ~~ c( 1, 1 ) * Ct2

            # error variances set at zero 
            # only ACE components explain variance
            xt1 ~~ 0 * xt1
            xt2 ~~ 0 * xt2

            # intercepts (all equal)
            xt1 ~ equal(c ( mu, mu ) ) * 1
            xt2 ~ equal(c ( mu, mu ) ) * 1
    
            # standardize
            totvar := h^2 + c^2 + e^2

            H2 := h^2/totvar
            C2 := c^2/totvar
            E2 := e^2/totvar

            check := H2 + C2 + E2   
            '

###### Model fitting & Results

# fit the model
result <- sem( ACEmodel, twindata, orthogonal = TRUE, group = "zyg", std.lv = TRUE )

# show results
summary( result, standardized = TRUE )
