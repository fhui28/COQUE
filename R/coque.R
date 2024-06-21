#' @title Penalized composite quadratic estimator (COQUE) for multivariate longitudinal/independent cluster GLMMs.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The COQUE objective function being optimized is \code{-0.5 * t(psi - psihat) %^% G %*% (psi - psihat) - lambda * penalty}, where \code{psihat} are the composite likelihood estimates produced from the stacked univariate GLMM (see the [stackedGLMMs()] function), \code{G} is the sandwich or Godambe information matrix, and \code{penalty} is a set of broken adaptive ridge (fusion) penalties to achieve simultaneous variable selection and coefficient clustering on the fixed effect coefficients. Construction of the regularization path for COQUE is done via parallel computing.
#'
#' **No penalization of the random effects is performed.**
#'
#'
#' @param object An object of class \code{stackedGLMM}.
#' @param nlambda The number of tuning parameters values.
#' @param lambda_min_ratio The smallest value for the tuning parameter lambda, as a fraction of the maximum internally-derived lambda value i.e., a value such that penalization leads to zeros for all penalized estimates.
#' @param lambda A user-supplied tuning parameter sequence. Note this is usually not supplied, as it is standard to have the function itself compute its own lambda sequence based on \code{nlambda} and \code{lambda_min_ratio}.
#' @param min_df The minimum number of non-zero fixed-effect coefficients, excluding the intercept, allowed in the model. This is useful to supply if, when the function tries to find a appropriate lambda sequence, the user wants the maximum value of lambda to not necessarily shrink to all the penalized estimates to zero. Defaults to zero, and can take a maximum value given by \code{length(selection_indices)}.
#' @param num_cores The regularization path is constructed via parallel computing, and this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to \code{parallel::detectCores() - 2}.
#' @param control A list containing the following elements:
#' \itemize{
#' \item{maxit: }{the maximum number of iterations in the algorithm.}
#' \item{eps: }{the convergence criterion; the norm of the difference between estimates from successive iterations must be smaller than this value.}
#' \item{round_eps: }{a tolerance to round values to zero; this is technically not needed as COQUE and the broken adaptive ridge penalty will produce exactly zero estimates up to machine error, but is included anyway.}
#' \item{fallback_to_bread: }{whether to use replace the sandwich information matrix with the observed information (i.e., the bread) matrix if instabilities are detected in the former.}
#' }
#'
#'
#' @details
#' See associated manuscript (currently in revision).
#'
#'
#' @return An object of class \code{coque} with the following elements (and not necessarily in the order below):
#' \item{call:}{The function call.}
#' \item{lambda:}{The actual sequence of tuning parameter, lambda, values used.}
#' \item{psi_path:}{A sparse matrix showing the estimated parameters at each point on the regularization path. The number of rows is equal \code{length(psihat)}, while the number of columns is equal to \code{length(lambda)}. The parameters are ordered so that all the fixed effect coefficients appear first (running as covariates within response), followed by any nuisance parameters and random effects variance parameters.}
#' \item{coque_statistics_path:}{A vector showing the value of the composite quadratic approximation statistic at each point on the regularization path.}
#' \item{df_path:}{A vector showing the number of non-zero estimated values in the (reparametrized) \code{psihat} at each point on the regularization path. This is used for calculating information criteria. Note \code{df_path} counts all parameters, not just those up for selection. This is consistent with the fact that the composite quadratic estimator is constructed with respect to all the parameters in the composite likelihood. }
#' \item{ICs:}{A matrix of some information criteria and their values on the regularization path. The number of rows is equal to \code{length(lambda)}. The current information criteria calculated include:
#' Akaike information criterion (AIC) with model complexity penalty of 2;
#' Bayesian information criterion (BIC) with model complexity penalty of \code{log(num_clus)} where \code{num_clus} is the numbers of clusters in the multivariate GLMM;
#' Extended Bayesian information criterion (EBIC) with model complexity penalty of \code{log(num_clus) + log(penalized parameters)};
#' Extended Bayesian information criterion (EBIC2) with model complexity penalty of \code{log(num_clus) + 2 * log(penalized parameters)}.}
#' \item{time_taken:}{The time taken to construct the regularization path.}
#' \item{num_resp/num_clus...:}{Additional auxilary function that can be safely ignored.}
#' \item{additional_information:}{This can be safely ignored.}
#'
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>
#'
#'
#' @seealso [estimates.coque()] for extracting values from a COQUE fit, [stackedGLMMs()] for fitting stacked univariate generalized linear mixed models (note this function has to be run upstream of running the main \code{coque} function), and [generatedata_mglmm()] for simulating data from a longitudinal/independent cluster multivariate generalized linear mixed model.
#'
#'
#' @examples
#' \dontrun{
#' ##------------------
#' # Example 1: Simulated binary responses
#' # Inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
#' ##------------------
#' ## Basic set up and true parameters
#' set.seed(052024)
#' library(mvtnorm)
#' num_clus <- 200
#' num_resp <- 6
#' response_type <- "binomial"
#' dat <- data.frame(id = rep(1:num_clus, sample(10:15, size = num_clus, replace = TRUE)))

#' p <- 15
#' H <- abs(outer(1:p, 1:p, "-"))
#' X <- rmvnorm(nrow(dat), sigma = 0.5^H)
#' colnames(X) <- paste0("x", 1:ncol(X))
#' dat <- cbind(dat, X)
#' rm(X)
#' str(dat)

#' true_fixed_effects <- cbind(
#' runif(num_resp, -1, 1),
#' rnorm(num_resp,0.5,sd=0.5),
#' rnorm(num_resp,-0.5,sd=0.5),
#' replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 3 unique values, including zero
#' replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 2 unique values, including zero
#' replicate(5, numeric(num_resp)))
#' true_fixed_effects
#'
#' true_G <- rbind(c(1,0.8,0.4,-0.6,-0.3,-0.15),
#' c(0.8,1,0.6,-0.4,-0.2,-0.1),
#' c(0.4,0.6,1,-0.2,-0.1,-0.05),
#' c(-0.6,-0.4,-0.2,1,0.6,0.3),
#' c(-0.3,-0.2,-0.1,0.6,1,0.3),
#' c(-0.15,-0.1,-0.05,0.3,0.3,1))
#'
#' true_ranef_cov <- rbind(c(2,1,0,0,0),
#' c(1,2,0,0,0),
#' c(0,0,1,0,0),
#' c(0,0,0,0.5,0),
#' c(0,0,0,0,0.25))/2
#'
#' true_Sigma <- kronecker(true_G, true_ranef_cov)
#'
#' rm(H, p)
#'
#' ## Simulate longitudinal/independent cluster multivariate responses
#' simy <- generatedata_mglmm(formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat,
#' fixed_effects = true_fixed_effects,
#' Sigma = true_Sigma,
#' response_type = response_type,
#' trial_size = 1)
#'
#' str(simy)
#'
#'
#' ## Apply COQUE -- **Real application starts here**
#'
#' # 1. Fitting stacked univariate GLMMs
#' sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
#' response_type = response_type,
#' formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat,
#' trial_size = 1)
#'
#' str(sglmms)
#'
#' # 2. Pass to the COQUE fit
#' penalized_fit <- coque(object = sglmms)
#' penalized_fit$lambda
#' penalized_fit$df
#' penalized_fit$ICs
#'
#' parameter_estimates <- estimates(penalized_fit, lambda = "BIC", hybrid = TRUE)
#'
#' true_fixed_effects
#' parameter_estimates$hybrid # Compare with true_fixed_effects
#' parameter_estimates$random_effects_covariance
#'
#'
#' ##------------------
#' # Example 2: Simulated Poisson responses
#' # Inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
#' ##------------------
#' ## Basic set up and true parameters
#' set.seed(052024)
#' library(mvtnorm)
#' num_clus <- 100
#' num_resp <- 6
#' response_type <- "poisson"
#' dat <- data.frame(id = rep(1:num_clus, sample(10:15, size = num_clus, replace = TRUE)))
#'
#' p <- 15
#' H <- abs(outer(1:p, 1:p, "-"))
#' X <- rmvnorm(nrow(dat), sigma = 0.5^H)
#' colnames(X) <- paste0("x", 1:ncol(X))
#' dat <- cbind(dat, X)
#' rm(X)
#' str(dat)
#'
#' true_fixed_effects <- cbind(
#' runif(num_resp, -1, 1),
#' rnorm(num_resp,0.5,sd=0.5),
#' rnorm(num_resp,-0.5,sd=0.5),
#' replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 3 unique values, including zero
#' replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 2 unique values, including zero
#' replicate(5, numeric(num_resp)))
#' true_fixed_effects
#'
#' true_G <- rbind(c(1,0.8,0.4,-0.6,-0.3,-0.15),
#' c(0.8,1,0.6,-0.4,-0.2,-0.1),
#' c(0.4,0.6,1,-0.2,-0.1,-0.05),
#' c(-0.6,-0.4,-0.2,1,0.6,0.3),
#' c(-0.3,-0.2,-0.1,0.6,1,0.3),
#' c(-0.15,-0.1,-0.05,0.3,0.3,1))
#'
#' true_ranef_cov <- rbind(c(2,1,0,0,0),
#' c(1,2,0,0,0),
#' c(0,0,1,0,0),
#' c(0,0,0,0.5,0),
#' c(0,0,0,0,0.25))/2
#' true_Sigma <- kronecker(true_G, true_ranef_cov)
#'
#' rm(H, p)
#'
#' ## Simulate longitudinal/independent cluster multivariate responses
#' simy <- generatedata_mglmm(formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat,
#' fixed_effects = true_fixed_effects,
#' Sigma = true_Sigma,
#' response_type = response_type)
#'
#' str(simy)
#'
#'
#' ## Apply COQUE -- **Real application starts here**
#'
#' # 1. Fitting stacked univariate GLMMs
#' sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
#' response_type = response_type,
#' formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat)
#'
#' str(sglmms)
#'
#' # 2. Pass to the COQUE fit
#' penalized_fit <- coque(object = sglmms)
#' penalized_fit$lambda
#' penalized_fit$df
#' penalized_fit$ICs
#'
#' parameter_estimates <- estimates(penalized_fit, lambda = "EBIC", hybrid = TRUE)
#'
#' true_fixed_effects
#' parameter_estimates$hybrid # Compare with true_fixed_effects
#' parameter_estimates$random_effects_covariance
#'
#'
#' ##------------------
#' # Example 3: Simulated Tweedie responses
#' # Inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
#' ##------------------
#' ## Basic set up and true parameters
#' set.seed(052024)
#' library(mvtnorm)
#' num_clus <- 200
#' num_resp <- 6
#' response_type <- "tweedie"
#' dat <- data.frame(id = rep(1:num_clus, sample(10:15, size = num_clus, replace = TRUE)))
#'
#' p <- 15
#' H <- abs(outer(1:p, 1:p, "-"))
#' X <- rmvnorm(nrow(dat), sigma = 0.5^H)
#' colnames(X) <- paste0("x", 1:ncol(X))
#' dat <- cbind(dat, X)
#' rm(X)
#' str(dat)
#'
#' true_fixed_effects <- cbind(
#' runif(num_resp, -1, 1),
#' rnorm(num_resp,0.5,sd=0.5),
#' rnorm(num_resp,-0.5,sd=0.5),
#' replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 3 unique values, including zero
#' replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 2 unique values, including zero
#' replicate(5, numeric(num_resp)))
#' true_fixed_effects
#'
#' true_G <- rbind(c(1,0.8,0.4,-0.6,-0.3,-0.15),
#' c(0.8,1,0.6,-0.4,-0.2,-0.1),
#' c(0.4,0.6,1,-0.2,-0.1,-0.05),
#' c(-0.6,-0.4,-0.2,1,0.6,0.3),
#' c(-0.3,-0.2,-0.1,0.6,1,0.3),
#' c(-0.15,-0.1,-0.05,0.3,0.3,1))
#'
#' true_ranef_cov <- rbind(c(2,1,0,0,0),
#' c(1,2,0,0,0),
#' c(0,0,1,0,0),
#' c(0,0,0,0.5,0),
#' c(0,0,0,0,0.25))/2
#' true_Sigma <- kronecker(true_G, true_ranef_cov)
#' true_phi <- runif(num_resp,0,2)
#' true_power <- runif(num_resp, 1.2, 1.8)
#'
#' rm(H, p)
#'
#' ## Simulate longitudinal/independent cluster multivariate responses
#' simy <- generatedata_mglmm(formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat,
#' fixed_effects = true_fixed_effects,
#' Sigma = true_Sigma,
#' response_type = response_type,
#' phi = true_phi,
#' power = true_power)
#'
#' str(simy)
#'
#'
#' ## Apply COQUE -- **Real application starts here**
#'
#' # 1. Fitting stacked univariate GLMMs
#' sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
#' response_type = response_type,
#' formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat)
#'
#' str(sglmms)
#'
#' # 2. Pass to the COQUE fit
#' penalized_fit <- coque(object = sglmms)
#' penalized_fit$lambda
#' penalized_fit$df
#' penalized_fit$ICs
#'
#' parameter_estimates <- estimates(penalized_fit, lambda = "BIC", hybrid = TRUE)
#'
#' true_fixed_effects
#' parameter_estimates$hybrid # Compare with true_fixed_effects
#' parameter_estimates$random_effects_covariance
#'
#'
#' ##------------------
#' # Example 4: Simulated ordinal responses
#' # Inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
#' ##------------------
#' ## Basic set up and true parameters
#' set.seed(052024)
#' library(mvtnorm)
#' num_clus <- 200
#' num_resp <- 6
#' response_type <- "ordinal"
#' dat <- data.frame(id = rep(1:num_clus, sample(10:15, size = num_clus, replace = TRUE)))
#'
#' p <- 15
#' H <- abs(outer(1:p, 1:p, "-"))
#' X <- rmvnorm(nrow(dat), sigma = 0.5^H)
#' colnames(X) <- paste0("x", 1:ncol(X))
#' dat <- cbind(dat, X)
#' rm(X)
#' str(dat)
#'
#' true_fixed_effects <- cbind(
#' runif(num_resp, -1, 1),
#' rnorm(num_resp,0.5,sd=0.5),
#' rnorm(num_resp,-0.5,sd=0.5),
#' replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 3 unique values, including zero
#' replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)),
#' # 3 covariates where in each covariate there are at most 2 unique values, including zero
#' replicate(5, numeric(num_resp)))
#' true_fixed_effects
#'
#' true_G <- rbind(c(1,0.8,0.4,-0.6,-0.3,-0.15),
#' c(0.8,1,0.6,-0.4,-0.2,-0.1),
#' c(0.4,0.6,1,-0.2,-0.1,-0.05),
#' c(-0.6,-0.4,-0.2,1,0.6,0.3),
#' c(-0.3,-0.2,-0.1,0.6,1,0.3),
#' c(-0.15,-0.1,-0.05,0.3,0.3,1))
#'
#' true_ranef_cov <- rbind(c(2,1,0,0,0),
#' c(1,2,0,0,0),
#' c(0,0,1,0,0),
#' c(0,0,0,0.5,0),
#' c(0,0,0,0,0.25))/2
#' true_Sigma <- kronecker(true_G, true_ranef_cov)
#' true_cutoffs <- cbind(0, matrix(c(1,2,3), nrow = num_resp, ncol = 3, byrow = TRUE))
#'
#' rm(H, p)
#'
#' ## Simulate longitudinal/independent cluster multivariate responses
#' simy <- generatedata_mglmm(formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat,
#' fixed_effects = true_fixed_effects,
#' Sigma = true_Sigma,
#' response_type = response_type,
#' cutoffs = true_cutoffs)
#'
#' str(simy)
#'
#'
#' ## Apply COQUE -- **Real application starts here**
#'
#' # 1. Fitting stacked univariate GLMMs
#' sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
#' response_type = response_type,
#' formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat)
#'
#' str(sglmms)
#'
#' # 2. Pass to the COQUE fit
#' penalized_fit <- coque(object = sglmms)
#' penalized_fit$lambda
#' penalized_fit$df
#' penalized_fit$ICs
#'
#' parameter_estimates <- estimates(penalized_fit, lambda = "BIC", hybrid = TRUE)
#'
#' true_fixed_effects
#' parameter_estimates$hybrid # Compare with true_fixed_effects
#' parameter_estimates$random_effects_covariance
#' }
#'
#'
#' @export coque
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom mvtnorm rmvnorm
#' @importFrom parallel detectCores
#' @importFrom quadprog solve.QP
#' @importFrom TMB MakeADFun
#' @importFrom stats cor cov fitted nlminb dnorm pnorm qnorm rnorm dbinom pbinom rbinom dnbinom pnbinom rnbinom dbeta pbeta rbeta dexp pexp rexp dgamma pgamma rgamma dlogis plogis qlogis dpois ppois rpois runif dchisq pchisq qchisq qqnorm as.formula binomial formula Gamma logLik model.matrix na.omit optim nlminb residuals
#' @import Matrix
#'
#'
#' @md

coque <- function(object, nlambda = 100, lambda_min_ratio = 1e-6, lambda = NULL, min_df = 0, num_cores = NULL,
                    control = list(maxit = 1000, eps = 1e-6, round_eps = 1e-6, fallback_to_bread = TRUE)) {

    ##-------------------
    ## Checks and balances
    ##-------------------
    if(class(object) != "stackedGLMM")
        stop("object must be of class 'stackedGLMM'.")
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in parallel to fit stacked GLMMs.")

    if(is.null(control$maxit))
        control$maxit <- 100
    if(is.null(control$eps))
        control$eps <- 1e-6
    if(is.null(control$round_eps))
        control$round_eps <- 1e-6

    # if(min_df >= length(selection_indices))
    #     stop("min_df, the smallest number of non-zero penalized estimates in the model, must not be larger than length(selection_indices).")

    num_clus <- nrow(object[[1]]$point_estimates$random_effects)
    fixed_effects_names <- colnames(model.matrix(object$formula_X[[1]], data = object$data))

    ##-------------------
    ## Extract composite (stacked) estimates and reorder them to put fixed effects first
    ##-------------------
    psihat <- do.call(c, lapply(object[1:object$num_resp], function(x) x$fit$par))
    new_order <- grep("fixed_effects", names(psihat))
    new_order <- c(new_order, (1:length(psihat))[-new_order])
    psihat <- psihat[new_order]


    ##-------------------
    ## Construct bread and meat matrices
    ##-------------------
    bread_matrix <- bdiag(lapply(object[1:object$num_resp], function(x) x$sdreport$cov.fixed))
    bread_matrix <- solve(bread_matrix) / num_clus
    bread_matrix <- bread_matrix[new_order, new_order]

    meat_matrix <- cov(do.call(cbind, lapply(object[1:object$num_resp], function(x) x$score_vectors)))
    meat_matrix <- solve(nearPD(meat_matrix)$mat) # Near PD is needed as sometimes it is not invertible...
    meat_matrix <- meat_matrix[new_order, new_order]
    sandwich_matrix <- as.matrix(bread_matrix %*% meat_matrix %*% bread_matrix)
    if(inherits(try(solve(sandwich_matrix), silent = TRUE), "try-error")) {
        if(control$fallback_to_bread)
            sandwich_matrix <- bread_matrix
        if(!control$fallback_to_bread)
            stop("The composite quadratic estimator can not be constructed due to instabilities. Note a way to remedy this is to set control$fallback_to_bread = TRUE.")
        }

   # d =  data.frame(
   #     est = psihat[1:(16*6)],
   #     true = c(sglmms[[1]]$sdreport$sd[1:16], sglmms[[2]]$sdreport$sd[1:16], sglmms[[3]]$sdreport$sd[1:16], sglmms[[4]]$sdreport$sd[1:16], sglmms[[5]]$sdreport$sd[1:16], sglmms[[6]]$sdreport$sd[1:16]),
   #     bread = (solve(bread_matrix)/num_clus) %>% diag %>% sqrt %>% {.[1:(16*6)]},
   #     sandwich = (solve(sandwich_matrix)/num_clus) %>% diag %>% sqrt %>% {.[1:(16*6)]})


    ##-------------------
    ## Identify fixed effects and perform relevant transformation to set up broken adaptive ridge fusion
    ##-------------------
    fixed_effects_matrix <- psihat[grep("fixed_effects", names(psihat))]
    if((length(fixed_effects_matrix) %% object$num_resp) != 0)
        stop("The broken adaptive ridge fusion penalty is designed specifically for the case where the fixed effects structure is the same for all responses.")
    fixed_effects_matrix <- matrix(fixed_effects_matrix, nrow = object$num_resp, byrow = TRUE)

    RSmatrix <- .reparametrizationmatrix(coefficients = fixed_effects_matrix)
    bigRSmatrix <- bdiag(RSmatrix, diag(nrow = length(psihat) - length(fixed_effects_matrix)))
    rm(RSmatrix)

    reparametrized_psihat <- as.vector(bigRSmatrix %*% psihat)
    reparametrized_sandwich_matrix <- crossprod(solve(bigRSmatrix), sandwich_matrix) %*% solve(bigRSmatrix) # same as reparametrized_sandwich_matrix <- solve(t(bigRSmatrix)) %*% sandwich_matrix %*% solve(bigRSmatrix)

    reparametrized_fixed_effects_matrix <- matrix(reparametrized_psihat[1:length(fixed_effects_matrix)], nrow = object$num_resp, byrow = TRUE)
    additional_weights <- matrix(1/abs(colSums(reparametrized_fixed_effects_matrix[-1,])), nrow = object$num_resp, ncol = ncol(fixed_effects_matrix), byrow = TRUE)
    additional_weights <- as.vector(t(additional_weights))
    additional_weights <- c(additional_weights, rep(0, length(psihat) - length(fixed_effects_matrix))) #' Note that even if you pass non-zero weights to the intercept, the selection_indices argument takes care of this!
    #(RSmatrix %*% as.vector(t(fixed_effects_matrix))) %>% matrix(nrow = num_resp, byrow = TRUE)
    #(RSmatrix %*% as.vector(t(true_fixed_effects))) %>% matrix(nrow = num_resp, byrow = TRUE)

    weights_mat <- matrix(1, nrow = object$num_resp, ncol = ncol(fixed_effects_matrix))
    weights_mat[,1] <- 0 #' Do not penalize the intercept
    weights_mat <- as.vector(t(weights_mat))
    selection_indices <- which(weights_mat == 1)
    rm(weights_mat, reparametrized_fixed_effects_matrix)

    #(abs(reparametrized_psihat[1:length(fixed_effects_matrix)]/sqrt((diag(solve(reparametrized_sandwich_matrix)/num_clus)))[1:length(fixed_effects_matrix)]) > 1.96) %>% matrix(nrow = num_resp, byrow = TRUE)


    ##-----------------------
    ## Attempt to determine a minimum and maximum lambda, if required, for broken adaptive ridge fusion (BARF)
    ##-----------------------
    if(!is.null(lambda))
        lambdaseq <- lambda
    if(is.null(lambda)) {
        message("Finding an appropriate sequence of lambda values...")
        lambda_max <- 1e-1
        fit_togetmaxlambda <- .BARF_singlelambda(psihat = reparametrized_psihat,
                                                 sandwich_matrix = reparametrized_sandwich_matrix,
                                                 selection_indices = selection_indices,
                                                 additional_weights = additional_weights,
                                                 lambda = lambda_max,
                                                 num_clus = num_clus,
                                                 control = control)
        fit_togetmaxlambda <- solve(bigRSmatrix) %*% fit_togetmaxlambda

        while(sum(fit_togetmaxlambda[selection_indices] != 0) > min_df) {
            lambda_max <- lambda_max*1.1
            fit_togetmaxlambda <- .BARF_singlelambda(psihat = reparametrized_psihat,
                                                     sandwich_matrix = reparametrized_sandwich_matrix,
                                                     selection_indices = selection_indices,
                                                     additional_weights = additional_weights,
                                                     lambda = lambda_max,
                                                     num_clus = num_clus,
                                                     control = control)
            }

        lambdaseq <- lseq(lambda_max, lambda_max*lambda_min_ratio, length = nlambda, decreasing = TRUE)
        rm(fit_togetmaxlambda)
        }


    ##-----------------------
    #' ## Construct regularization path for broken adaptive ridge fusion (BARF)
    ##-----------------------
    tic <- proc.time()
    message("Constructing regularization path...")
    cwfit_fn <- function(l) {
        cwfit <- .BARF_singlelambda(psihat = reparametrized_psihat,
                                    sandwich_matrix = reparametrized_sandwich_matrix,
                                    selection_indices = selection_indices,
                                    additional_weights = additional_weights,
                                    lambda = lambdaseq[l],
                                    num_clus = num_clus,
                                    control = control)

        out <- list(
            cw_coefficients = cwfit,
            df_path = sum(cwfit != 0),
            coque_statistics = num_clus * as.vector(crossprod(cwfit - reparametrized_psihat, reparametrized_sandwich_matrix) %*% (cwfit - reparametrized_psihat)) #' Analogue of the multivariate Wald-statistic in [https://www.jstor.org/stable/24309261]
            )
        out$AIC <- out$coque_statistics + 2*out$df_path
        out$BIC <- out$coque_statistics + log(nrow(object$data))*out$df_path
        out$EBIC <- out$coque_statistics + (log(nrow(object$data)) + log(length(selection_indices)))*out$df_path
        out$EBIC2 <- out$coque_statistics + (log(nrow(object$data)) + 2 * log(length(selection_indices)))*out$df_path
        #out$ERIC <- out$coque_statistics - log(lambdaseq[l])*out$df_path
        #out$ERIC2 <- out$coque_statistics - 2*log(lambdaseq[l])*out$df_path

        return(out)
        }

    allfits <- foreach(l = 1:length(lambdaseq)) %dopar% cwfit_fn(l = l)
    toc <- proc.time()

    coefficients_path <- Matrix(sapply(allfits, function(x) x$cw_coefficients), sparse = TRUE)
    rownames(coefficients_path) <- names(reparametrized_psihat)
    df_path <- sapply(allfits, function(x) x$df_path)
    coque_statistics <- sapply(allfits, function(x) x$coque_statistics)
    AIC <- sapply(allfits, function(x) x$AIC)
    BIC <- sapply(allfits, function(x) x$BIC)
    EBIC <- sapply(allfits, function(x) x$EBIC)
    EBIC2 <- sapply(allfits, function(x) x$EBIC2)
    #ERIC <- sapply(allfits, function(x) x$ERIC)
    #ERIC2 <- sapply(allfits, function(x) x$ERIC2)


    out <- list(lambda = lambdaseq,
                psi_path = coefficients_path,
                coque_statistics_path = coque_statistics,
                df_path = df_path,
                ICs = data.frame(AIC, BIC, EBIC, EBIC2),
                time_taken = toc - tic)


    ##-----------------------
    ## Final touch up
    ##-----------------------
    out$psi_path <- solve(bigRSmatrix) %*% out$psi_path
    rownames(out$psi_path) <- names(psihat)

    out$call <- match.call()
    out$num_resp <- object$num_resp
    out$num_clus <- num_clus
    out$response_type <- object$response_type
    out$response_names <- names(object)[1:object$num_resp]
    out$fixed_effects_names <- fixed_effects_names

    out$additional_information <- list(reparametrized_psihat = reparametrized_psihat, reparametrized_sandwich_matrix = reparametrized_sandwich_matrix, reparametrized_psi_path = coefficients_path, bigRSmatrix = bigRSmatrix)
    class(out) <- "coque"
    return(out)
    }



#' @title Hidden internal function for applying the broken adaptive ridge (fusion) penalty applied to a composite form approximation -- single tuning parameter version.
#' @return A vector of penalized estimates with the same length as \code{psihat}, with potentially one or more elements shrunk to zero.
#' @noRd
#' @noMd

.BARF_singlelambda <- function(psihat, sandwich_matrix, selection_indices, additional_weights, lambda, num_clus, control) {

    ##-------------------------
    #' ## Attempting the analogue of equation 4 in [https://doi.org/10.1016/j.jmva.2018.08.007]
    ##-------------------------
    Dbar <- Diagonal(x = additional_weights, n = length(psihat))
    diag(Dbar)[-selection_indices] <- 0
    err <- Inf
    counter <- 0

    cw_psi <- psihat
    while(counter < control$maxit & err > control$eps) {
        GammaMatrix_psi <- Diagonal(n = length(psihat))
        diag(GammaMatrix_psi)[selection_indices] <- cw_psi[selection_indices]

        new_psi <- GammaMatrix_psi %*% solve(GammaMatrix_psi %*% sandwich_matrix %*% GammaMatrix_psi + lambda * Dbar) %*% GammaMatrix_psi %*% sandwich_matrix %*% psihat
        new_psi <- as.vector(new_psi)
        names(new_psi) <- names(psihat)

        err <- sum((new_psi - cw_psi)^2)
        cw_psi <- new_psi
        counter <- counter + 1
    }

    #' Note the convergence properties of the algorithm does basically ensure exact zeros can be achieved (up to machine error); see [https://doi.org/10.1155/2016/3456153]. However once can also round to zero if interested. In this case, I will round based on eps = 1e-6
    new_psi[abs(new_psi) < control$round_eps] <- 0
    return(new_psi)
}


