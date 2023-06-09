% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity.metrics.R
\name{sim.metrics}
\alias{sim.metrics}
\title{Study Strap similarity measures: Supporting function used as the default similarity measures in Study Strap, SSE, and CMSS algorithms. Compares similarity in covaraite profiles of 2 studies.}
\usage{
sim.metrics(dat1, dat2)
}
\arguments{
\item{dat1}{A design matrix of the first study.}

\item{dat2}{A design matrix of the second study to be compared to the first study.}
}
\value{
Vector of similarity measures.
}
\description{
Study Strap similarity measures: Supporting function used as the default similarity measures in Study Strap, SSE, and CMSS algorithms. Compares similarity in covaraite profiles of 2 studies.
}
\examples{
set.seed(1)

##########################
##### Simulate Data ######
##########################

# create training dataset with 10 studies, 2 covariates
X <- matrix(rnorm(2000), ncol = 2)

# true beta coefficients
B <- c(5, 10, 15)

# outcome vector
y <- cbind(1, X) \%*\% B

# study names
study <- sample.int(10, 1000, replace = TRUE)
data <- data.frame( Study = study,
                    Y = y,
                    V1 = X[,1],
                    V2 = X[,2] )


# create target study design matrix for
# covariate profile similarity weighting and
# accept/reject algorithm (covaraite-matched study strap)

target <- matrix(rnorm(1000), ncol = 2) # design matrix only
colnames(target) <- c("V1", "V2")

#############################
#### Similarity Measures ####
#############################
# compare the covariate profile of the entire training dataset with that of the target study.

sim.vec <- sim.metrics(target, data[-c(1,2)])
}
