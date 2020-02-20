## ---- warning=FALSE, comment=FALSE--------------------------------------------
set.seed(1)
library(studyStrap)
# create half of training dataset from 1 distribution
X1 <- matrix(rnorm(2000), ncol = 2) # design matrix - 2 covariates
B1 <- c(5, 10, 15) # true beta coefficients
y1 <- cbind(1, X1) %*% B1

# create 2nd half of training dataset from another distribution
X2 <- matrix(rnorm(2000, 1,2), ncol = 2) # design matrix - 2 covariates
B2 <- c(10, 5, 0) # true beta coefficients
y2 <- cbind(1, X2) %*% B2

X <- rbind(X1, X2)
y <- c(y1, y2)

study <- sample.int(10, 2000, replace = TRUE) # 10 studies
data <- data.frame( Study = study, Y = y, V1 = X[,1], V2 = X[,2] )


# create target study design matrix for covariate profile similarity weighting and 
# accept/reject algorithm (Covariate-matched study strap)
target <- matrix(rnorm(1000, 3, 5), ncol = 2) # design matrix
colnames(target) <- c("V1", "V2")

## -----------------------------------------------------------------------------
head(data)

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# custom function
fn1 <- function(x1,x2){
    return( abs( cor( colMeans(x1), colMeans(x2) )) )
    } 

sseMod1 <- sse(formula = Y ~., 
               data = data, 
               target.study = target,
               ssl.method = list("pcr"), 
               ssl.tuneGrid = list(data.frame("ncomp" = 1)), 
               customFNs = list(fn1) )


## ---- warning=FALSE, comment=FALSE, warning=FALSE-----------------------------
preds <- studyStrap.predict(sseMod1, target)
head(preds)[1:3,]

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# custom function
fn1 <- function(x1,x2){
    return( abs( cor( colMeans(x1), colMeans(x2) )) )
    } 

sseMod2 <- sse(formula = Y ~., 
               data = data, 
               target.study = target,
               ssl.method = list("lm","pcr"), 
               ssl.tuneGrid = list(NA, data.frame("ncomp" = 2)), 
               customFNs = list(fn1) )


## ---- warning=FALSE, comment=FALSE, warning=FALSE-----------------------------
preds <- studyStrap.predict(sseMod2, target)
head(preds)[1:3,]

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
sseMod3 <- sse(formula = Y ~., 
               data = data,
               ssl.method = list("pcr"), 
               ssl.tuneGrid = list(NA, data.frame("ncomp" = 1)), 
               sim.mets = FALSE)

preds <- studyStrap.predict(sseMod3, target)
head(preds)[1:3,]


## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# 1 SSL
mrgMod1 <- merged(formula = Y ~.,
                  data = data,   
                  ssl.method = list("pcr"), 
                  ssl.tuneGrid = list( data.frame("ncomp" = 2))
                  )

# 2 SSLs
mrgMod2 <- merged(formula = Y ~.,
                  data = data,  
                  ssl.method = list("lm","pcr"), 
                  ssl.tuneGrid = list(NA, data.frame("ncomp" = 2))
                  )


## ---- warning=FALSE, comment=FALSE, warning=FALSE-----------------------------
preds <- studyStrap.predict(mrgMod2, target)
head(preds)

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# custom function
fn1 <- function(x1,x2){
    return( abs( cor( colMeans(x1), colMeans(x2) )) )
    } 

# 1 SSL
ssMod1 <- ss(formula = Y ~.,
             data = data,  
             target.study = target,
             bag.size = length(unique(data$Study)), 
             straps = 10, 
             stack = "standard",
             sim.covs = NA, 
             ssl.method = list("pcr"), 
             ssl.tuneGrid = list(data.frame("ncomp" = 2)), 
             sim.mets = TRUE,
             model = TRUE, 
             customFNs = list( fn1 ) )

# 2 SSLs
ssMod2 <- ss(formula = Y ~., 
             data = data, 
             target.study = target,
             bag.size = length(unique(data$Study)), 
             straps = 10, 
             stack = "standard",
             sim.covs = NA, 
             ssl.method = list("lm","pcr"), 
             ssl.tuneGrid = list(NA, data.frame("ncomp" = 2)), 
             sim.mets = TRUE,
             model = TRUE, 
             customFNs = list( fn1 ) )

## ---- warning=FALSE, comment=FALSE, warning=FALSE-----------------------------
preds <- studyStrap.predict(ssMod2, target)
head(preds)[1:3,]

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# custom function
fn1 <- function(x1,x2){
    return( abs( cor( colMeans(x1), colMeans(x2) )) )
    } 

ssMod3 <- ss(formula = Y ~., 
             data = data, 
             target.study = target,
             bag.size = length(unique(data$Study)), 
             straps = 10, 
             sim.covs = NA, ssl.method = list("pcr"), 
             ssl.tuneGrid = list(data.frame("ncomp" = 2)), 
             sim.mets = FALSE,
             customFNs = list( fn1 ) )

preds <- studyStrap.predict(ssMod3, target)
head(preds)[1:3,]

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------

ssMod4 <- ss(formula = Y ~., 
             data = data, 
             bag.size = length(unique(data$Study)), 
             straps = 10, 
             sim.covs = NA, ssl.method = list("pcr"), 
             ssl.tuneGrid = list(data.frame("ncomp" = 2)), 
             sim.mets = FALSE)

preds <- studyStrap.predict(ssMod4, target)
head(preds)[1:3,]

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# 1 SSL
arMod1 <-  cmss(formula = Y ~., 
                data = data, 
                target.study = target,
                converge.lim = 2,
                bag.size = length(unique(data$Study)), 
                max.straps = 50, 
                paths = 2, 
                ssl.method = list("pcr"), 
                ssl.tuneGrid = list(data.frame("ncomp" = 2))
                )

# 2 SSLs
arMod2 <-  cmss(formula = Y ~., 
                data = data, 
                target.study = target,
                converge.lim = 2,
                bag.size = length(unique(data$Study)), 
                max.straps = 50, 
                paths = 2, 
                ssl.method = list("lm","pcr"), 
                ssl.tuneGrid = list(NA, data.frame("ncomp" = 2))
                )

preds <- studyStrap.predict(arMod2, target)
head(preds)[1:3,]

## ---- warning=FALSE, comment=FALSE, message = FALSE---------------------------
# 1 SSL

# custom function for CPS
fn1 <- function(x1,x2){
    return( abs( cor( colMeans(x1), colMeans(x2) )) )
} 

# custom function for Accept/Reject step criteria
fn2 <- function(x1,x2){
    return( sum ( ( colMeans(x1) - colMeans(x2) )^2 ) )
    } 

arMod3 <-  cmss(formula = Y ~., 
                data = data, 
                target.study = target,
                converge.lim = 2,
                bag.size = length(unique(data$Study)), 
                max.straps = 50, 
                paths = 2, 
                ssl.method = list("pcr"), 
                ssl.tuneGrid = list(data.frame("ncomp" = 2)),
                sim.mets = FALSE,
                sim.fn = fn2,
                customFNs = list( fn1, fn2 ) 
                )

preds <- studyStrap.predict(arMod3, target)
head(preds)[1:3,]

## -----------------------------------------------------------------------------
sseMod1

## -----------------------------------------------------------------------------
sseMod1$models

## -----------------------------------------------------------------------------
names(sseMod1$modelInfo)

## -----------------------------------------------------------------------------
names(sseMod1$dataInfo)

