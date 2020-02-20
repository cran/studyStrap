#' fatTrim: Supporting function to reduce the size of models
#'
#' @param cmx A model object.
#' @return A model object.
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
#' ##########################
#' ##### Model Fitting #####
#' ##########################
#'
#' # Fit model with 1 Single-Study Learner (SSL): Linear Regression
#' mod1 <- lm(formula = Y ~., data = data)
#'
#'
#' ############################################
#' ##### Fat Trim to reduce model size #####
#' ############################################
#'
#' mod1.trim <- fatTrim(mod1)
#'
#' # compare sizes
#' object.size(mod1)
#' object.size(mod1.trim)
#' @export


fatTrim = function(cmx) {
    # modified from http://www.win-vector.com/blog/2014/05/trimming-the-fat-from-glm-models-in-r/
    cmx$y = c()
    cmx$model = c()
    cmx$residuals = c()
    cmx$fitted.values = c()
    cmx$effects = c()
    cmx$scores = c()
    cmx$loadings = c()
    cmx$weights = c()
    cmx$Yloadings = c()
    cmx$Xtotvar = c()
    cmx$trainingData = c()
    cmx$resample = c()
    cmx$results = c()
    #cmx$control = c()
    cmx$dots = c()
    cmx$times = c()


    attr(cmx$terms,".Environment") = c()
    attr(cmx$formula,".Environment") = c()

    cmx
}
