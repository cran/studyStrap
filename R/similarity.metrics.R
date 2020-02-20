#' Study Strap similarity measures: Supporting function used as the default similarity measures in Study Strap, SSE, and CMSS algorithms. Compares similarity in covaraite profiles of 2 studies.
#'
#' @param dat1 A design matrix of the first study.
#' @param dat2 A design matrix of the second study to be compared to the first study.
#' @return Vector of similarity measures.
#' @examples
#' set.seed(1)
#'
#' ##########################
#' ##### Simulate Data ######
#' ##########################
#'
#' # create training dataset with 10 studies, 2 covariates
#' X <- matrix(rnorm(2000), ncol = 2)
#'
#' # true beta coefficients
#' B <- c(5, 10, 15)
#'
#' # outcome vector
#' y <- cbind(1, X) %*% B
#'
#' # study names
#' study <- sample.int(10, 1000, replace = TRUE)
#' data <- data.frame( Study = study,
#'                     Y = y,
#'                     V1 = X[,1],
#'                     V2 = X[,2] )
#'
#'
#' # create target study design matrix for
#' # covariate profile similarity weighting and
#' # accept/reject algorithm (covaraite-matched study strap)
#'
#' target <- matrix(rnorm(1000), ncol = 2) # design matrix only
#' colnames(target) <- c("V1", "V2")
#'
#' #############################
#' #### Similarity Measures ####
#' #############################
#' # compare the covariate profile of the entire training dataset with that of the target study.
#'
#' sim.vec <- sim.metrics(target, data[-c(1,2)])
#' @import CCA
#' @import MatrixCorrelation
#' @importFrom stats coef cor cov predict
#' @export

sim.metrics <- function(dat1, dat2){
    # Takes in 2 design matrices (outcome vector is not attached) to compare and returns vector of
    # similarity metrics that when combined for each iteration that can be used to calculate weights
    #   23 metrics

    # library(CCA)
    # library(MatrixCorrelation)
    similarity.vec <- c()


    # number of covariates
    covs <- ncol(dat1)
    n <- nrow(dat1)
    p <- ncol(dat1)
    q <- ncol(dat2)

    ############################################
    # Correlation Coefficient Similarity Metric
    ############################################
    # https://stats.idre.ucla.edu/r/dae/canonical-correlation-analysis/
    ## UCLA

    # correlation matrix (make sure rows match up)
    mc1 <- matcor(dat1[1:min(nrow(dat1), nrow(dat2)),], dat2[1:min(nrow(dat1), nrow(dat2)),]) # correlation matrix

    cfs <- diag( mc1[[3]][1:p, (q+1):(p+q)] ) # correlation coefficients between corresponding values

    similarity.metric.corr.abs <- mean(abs(cfs)) # the higher the better
    similarity.metric.corr.sq <- sum(mc1[[3]][1:p, (q+1):(p+q)]) # the higher the better
    similarity.metric.corr.sq.abs <- sum(abs(mc1[[3]][1:p, (q+1):(p+q)]))
    # add to vector of similarity metrics
    similarity.vec <- c(similarity.vec, similarity.metric.corr.abs,
                        similarity.metric.corr.sq, similarity.metric.corr.sq.abs)


    ############################################
    # Canonical Correlations Similarity Metric
    ############################################
    ############################################
    # Covariance of Canonical Variables
    ############################################
    # http://users.stat.umn.edu/~helwig/notes/cancor-Notes.pdf
    # Nathaniel E. Helwig
    # Univ of Minnesota

    # standardize data
    Xs <- scale(dat1[1:min(nrow(dat1), nrow(dat2)),])
    Ys <- scale(dat2[1:min(nrow(dat1), nrow(dat2)),])

    # cca (the normal way)
    Sx <- cov(Xs)
    Sy <- cov(Ys)
    Sxy <- cov(Xs,Ys)
    Sxeig <- eigen(Sx, symmetric=TRUE)
    Sxisqrt <- Sxeig$vectors %*% diag(1/sqrt(Sxeig$values)) %*% t(Sxeig$vectors)
    Syeig <- eigen(Sy, symmetric=TRUE)
    Syisqrt <- Syeig$vectors %*% diag(1/sqrt(Syeig$values)) %*% t(Syeig$vectors)
    Xmat <- Sxisqrt %*% Sxy %*% solve(Sy) %*% t(Sxy) %*% Sxisqrt
    Ymat <- Syisqrt %*% t(Sxy) %*% solve(Sx) %*% Sxy %*% Syisqrt
    Xeig <- eigen(Xmat, symmetric=TRUE)
    Yeig <- eigen(Ymat, symmetric=TRUE)


    #-------------------------------------------
    # Canonical Correlations Similarity Metric
    #-------------------------------------------
    rho <- sqrt(Xeig$values)

    similarity.metric.can <- mean(rho)
    similarity.metric.can.abs <- mean(abs(rho))
    similarity.metric.can.sq <- mean((rho)^2)

    # add to vector of similarity metrics
    similarity.vec <- c(similarity.vec, similarity.metric.can,
                        similarity.metric.can.sq)
    #-------------------------------------------

    #-------------------------------------------
    # Covariance of Canonical Variables
    #-------------------------------------------

    #
    Ahat <- Sxisqrt %*% Xeig$vectors
    Bhat <- Syisqrt %*% Yeig$vectors
    Ainv <- solve(Ahat)
    Binv <- solve(Bhat)

    # canonical variables
    U <- Xs %*% Ahat
    V <- Ys %*% Bhat

    # covariance of canonical variables (U and V)
    rhomat <- cbind(diag(rho), matrix(0, p, q-p))
    similarity.metric.can.UV.rho <- sum( (cov(U, V) - rhomat)^2 )
    similarity.metric.can.UV.cov <- sum( (Sxy - crossprod(Ainv, rhomat) %*% Binv)^2 )
    sim4 <- sum( (cov(U, V) - rhomat) )
    sim5 <- sum( (Sxy - crossprod(Ainv, rhomat) %*% Binv) ) #So far the best
    sim6 <- sum(diag((cov(U, V) - rhomat)^2))
    sim7 <- sum(diag((Sxy - crossprod(Ainv, rhomat) %*% Binv)))
    sim8 <- sum(diag((Sxy - crossprod(Ainv, rhomat) %*% Binv))^2)

    mean.cor <- abs(cor(colMeans(dat1), colMeans(dat2)))

    similarity.vec <- c(similarity.vec, similarity.metric.can.UV.rho, similarity.metric.can.UV.cov, sim4, sim5, sim6, sim7, sim8, mean.cor)

    #-------------------------------------------
    # Matrix Correlations
    #-------------------------------------------

    # added from MatrixCorrelations package
    ## number of components is arbitrary--default choses the number of columns of each matrix to try
    comp1 <- p
    comp2 <- q
    similarity.vec <- c(similarity.vec, (allCorrelations(Xs, Ys,
                                            ncomp1 = comp1, ncomp2 = comp2, plot = FALSE)) )

    return(similarity.vec)
}
