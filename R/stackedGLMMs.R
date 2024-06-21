#' @title Fit stacked univariate generalized linear mixed models (GLMMs).
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function can be used on its own, but in the context of this the package, it is primarily utilized upstream of applying [coque()] for simultaneous selection and clustering of the fixed effects. Fitting of the GLMMs is done via [TMB::MakeADFun()], which optimizes the Laplace-approximated marginal log-likelihood function using automatic differentiation.
#'
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial", "ordinal", "tweedie".
#' @param formula_X One-sided formula for the fixed effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{fixed_effects} argument.
#' @param formula_Z One-sided formula for the random effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{Sigma} argument.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data frame.
#' @param trial_size A vector of response-specific trial sizes, which is needed for binomial distributed responses.
#' @param start_params Starting values can be supplied. This must take the form of a list with length equal to the number of responses i.e., the number of columns in \code{y}. Each element in the list should then contain:
#' \itemize{
#' \item{fixed_effects: }{a vector of starting values for the fixed effect coefficients.}
#' \item{random_effects: }{a matrix of starting values for the random effect coefficients where the number of columns is equal to the number of columns of the model matrix implied by \code{formula_Z}.}
#' \item{sd: }{a vector of standard deviation parameters for the random effects covariance matrix.}
#' \item{unconstrained_cor_params: }{a vector of parameters corresponding to the lower-triangular elements of the Cholesky decomposition of the random effects correlation matrix. We refer the reader to \url{https://kaskr.github.io/adcomp/_book/Densities.html} to see how covariance matrices are parametrized in \code{TMB}.}
#' }
#' @param TMB_directories A list with two elements, identifying the directory where TMB C++ file exists (\code{cpp}), and the directory where the corresponding compiled files to be placed (\code{compile}). *Unless you really want to do some real mucking around, these should be left at their default i.e., the directory where the packages were installed locally. Please note a version of the C++ file will be copied to the \code{compile} directory.*
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{parallel::detectCores()-2}.
#' @param score_vectors Should vectors of the score values for each cluster within each responses be returned. Defaults to \code{TRUE} and this generally should *not* be turned off.
#' @param control A list controlling the some of the tuning parameters in the optimization process for fitting the stacked GLMMs, noting [stats::nlminb()] is used to perform the optimization. These include:
#' \itemize{
#' \item{iter.max: }{the maximum number of iterations allowed.}
#' \item{eval.max: }{the maximum number of evaluations of the objective function allowed.}
#' \item{trace: }{this can be set to a positive value so that the value of the objective function and the parameters is printed every \code{trace}'th iteration.}
#' \item{attempts: }{the number of attempts are made to fit each stacked univariate GLMM (if [stats::nlminb()] fails on the first try).}
#' }
#'
#'
#' @details
#' See associated manuscript (currently in revision).
#'
#'
#' @return A object for class "stackedGLMM" which contains the following:
#' \item{respxxx: }{The first \code{ncol(y)} elements contains the estimated stacked GLMMs corresponding to each response. Each element is itself a list containing: the fitted model output from TMB (\code{fit} and \code{sdreport}), a list of estimated parameters and random effects in a slightly tidier format (\code{point_estimates}), and a matrix of score vector values for each cluster (\code{score_vectors}).}
#' \item{call: }{The function call.}
#' \item{formula_X/formula_Z/data/response_type/num_resp: }{Arugment information that can be safely ignored.}
#'
#'
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>
#'
#'
#' @seealso [coque()] for the main function for constructing a COQUE fit, [estimates.coque()] for extracting values from a COQUE fit, and [generatedata_mglmm()] for simulating data from a longitudinal/independent cluster multivariate generalized linear mixed model#'
#'
#' @examples
#' \dontrun{
#' ##------------------
#' # Example 1: Simulated binary responses
#' # Inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
#' ##------------------
#' ## Basic set up and true parameters
#' set.seed(052024)
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
#' ## Fitting stacked univariate GLMMs -- **Real application starts here**
#'
#' sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
#' response_type = response_type,
#' formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat,
#' trial_size = 1)
#'
#' str(sglmms)
#'
#' sglmms$resp1$point_estimates # Example of point estimates from response 1
#'
#'
#' ##------------------
#' # Example 2: Simulated Poisson responses
#' # Inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
#' ##------------------
#' ## Basic set up and true parameters
#' set.seed(052024)
#' num_clus <- 200
#' num_resp <- 6
#' response_type <- "poisson"
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
#' response_type = response_type)
#'
#' str(simy)
#'
#'
#' ## Fitting stacked univariate GLMMs -- **Real application starts here**
#'
#' sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
#' response_type = response_type,
#' formula_X = ~ . - id,
#' formula_Z = ~ x1 + x2 + x3 + x4,
#' data = dat)
#'
#' str(sglmms)
#'
#' sglmms$resp3$point_estimates # Example of point estimates from response 3
#' }
#'
#'
#' @export stackedGLMMs
#' @importFrom TMB sdreport
#' @md


stackedGLMMs <- function(y, response_type, formula_X, formula_Z, data, trial_size = 1,
                         start_params = NULL,
                         TMB_directories = list(cpp = system.file("executables", package = "COQUE"), compile = system.file("executables", package = "COQUE")),
                         num_cores = NULL, score_vectors = TRUE,
                         control = list(iter.max = 1000, eval.max = 1000, trace = 0, attempts = 20)) {

    ##---------------------
    ## Checks and balances
    ##---------------------
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in parallel to fit stacked GLMMs.")

    if(is.null(control$iter.max))
        control$iter.max <- 1000
    if(is.null(control$eval.max))
        control$eval.max <- 1000
    if(is.null(control$trace))
        control$trace <- 0
    if(is.null(control$attempts))
        control$attempts <- 20

    response_type <- match.arg(response_type, choices = c("gaussian", "poisson", "binomial", "ordinal", "tweedie"))
    if(!(response_type %in% c("gaussian", "poisson", "binomial", "ordinal", "tweedie")))
        stop("Specified family not permitted. Sorry!")

    y <- as.matrix(y)
    num_resp <- ncol(y)

    if(!is.list(formula_X))
        formula_X <- replicate(num_resp, formula_X, simplify = FALSE)
    if(!is.list(formula_Z))
        formula_Z <- replicate(num_resp, formula_Z, simplify = FALSE)

    if(length(formula_X) != length(formula_Z))
        stop("formula_X and formula_X should either both be a single formula, or lists of the same length. In the case of the latter, each element corresponds to the formula to use for the fixed/random effects in the corresponding column of y.")

    full_dat <- na.omit(cbind(y, data))
    y <- full_dat[,1:num_resp]
    data <- full_dat[,-(1:num_resp)]
    data$id <- as.numeric(as.factor(data$id))
    rm(full_dat)


    ##---------------------
    ## Compile TMB C++ script -- doing it outside of the external functions below to save time and potential conflicts
    ##---------------------
    getDLL <- paste0("glmmMLE_", response_type)
    if(control$trace > 0) {
        message("Compiling TMB C++ file...")
        }

    file.copy(from = paste0(TMB_directories$cpp, "/", getDLL, ".cpp"), to = paste0(TMB_directories$compile, "/", getDLL, ".cpp"), overwrite = FALSE)

    origwd <- getwd()
    # Please see https://github.com/kaskr/adcomp/issues/321 for flags argument
    setwd(TMB_directories$compile)
    TMB::compile(paste0(getDLL, ".cpp"), flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
    dyn.load(paste0(TMB_directories$compile, "/", TMB::dynlib(getDLL)))
    setwd(origwd)


    ##---------------------
    ## Set up and fit the stacked GLMMs
    ##---------------------
    dofits <- function(select_response) {
        maketidbits <- .construct_tidbits(response_type = response_type,
                                          formula_X = formula_X[[select_response]],
                                          formula_Z = formula_Z[[select_response]],
                                          data = data,
                                          y = y[,select_response],
                                          start_params = start_params[[select_response]],
                                          trial_size = trial_size)

        glmm_MLE_obj <- MakeADFun(data = maketidbits$data,
                                  parameters = maketidbits$parameters,
                                  random = maketidbits$random,
                                  DLL = maketidbits$DLL,
                                  silent = TRUE)

        glmm_MLE_obj$lower_limits <- rep(-Inf, length(names(glmm_MLE_obj$par)))
        glmm_MLE_obj$lower_limits[grep("phi", names(glmm_MLE_obj$par))] <- 1e-3
        glmm_MLE_obj$lower_limits[grep("powerparam", names(glmm_MLE_obj$par))] <- 1+1e-3
        glmm_MLE_obj$lower_limits[grep("sd", names(glmm_MLE_obj$par))] <- 1e-3
        glmm_MLE_obj$lower_limits[grep("unconstrained_cor_params", names(glmm_MLE_obj$par))] <- -1e3
        if(response_type == "ordinal")
            glmm_MLE_obj$lower_limits[grep("difford", names(glmm_MLE_obj$lower_limits))] <- 1e-3

        glmm_MLE_obj$upper_limits <- rep(Inf, length(names(glmm_MLE_obj$par)))
        glmm_MLE_obj$upper_limits[grep("powerparam", names(glmm_MLE_obj$par))] <- 2-1e-3
        glmm_MLE_obj$upper_limits[grep("unconstrained_cor_params", names(glmm_MLE_obj$par))] <- 1e3

        glmmfit <- try(nlminb(start = glmm_MLE_obj$par,
                          objective = glmm_MLE_obj$fn,
                          gradient = glmm_MLE_obj$gr,
                          lower = glmm_MLE_obj$lower_limits,
                          upper = glmm_MLE_obj$upper_limits,
                          control = control),
                       silent = TRUE)

        attempt_counter <- 1
        while(inherits(glmmfit, "try-error") & attempt_counter < control$attempts) {
            maketidbits <- .construct_tidbits(response_type = response_type,
                                              formula_X = formula_X[[select_response]],
                                              formula_Z = formula_Z[[select_response]],
                                              data = data,
                                              y = y[,select_response],
                                              trial_size = trial_size)

            glmm_MLE_obj <- MakeADFun(data = maketidbits$data,
                                      parameters = maketidbits$parameters,
                                      random = maketidbits$random,
                                      DLL = maketidbits$DLL,
                                      silent = TRUE)


            glmm_MLE_obj$lower_limits <- rep(-Inf, length(names(glmm_MLE_obj$par)))
            glmm_MLE_obj$lower_limits[grep("phi", names(glmm_MLE_obj$par))] <- 1e-3
            glmm_MLE_obj$lower_limits[grep("powerparam", names(glmm_MLE_obj$par))] <- 1+1e-3
            glmm_MLE_obj$lower_limits[grep("sd", names(glmm_MLE_obj$par))] <- 1e-3
            glmm_MLE_obj$lower_limits[grep("unconstrained_cor_params", names(glmm_MLE_obj$par))] <- -1e3
            if(response_type == "ordinal")
                glmm_MLE_obj$lower_limits[grep("difford", names(glmm_MLE_obj$lower_limits))] <- 1e-3

            glmm_MLE_obj$upper_limits <- rep(Inf, length(names(glmm_MLE_obj$par)))
            glmm_MLE_obj$upper_limits[grep("powerparam", names(glmm_MLE_obj$par))] <- 2-1e-3
            glmm_MLE_obj$upper_limits[grep("unconstrained_cor_params", names(glmm_MLE_obj$par))] <- 1e3

            glmmfit <- try(nlminb(start = glmm_MLE_obj$par,
                                  objective = glmm_MLE_obj$fn,
                                  gradient = glmm_MLE_obj$gr,
                                  lower = glmm_MLE_obj$lower_limits,
                                  upper = glmm_MLE_obj$upper_limits,
                                  control = control),
                           silent = TRUE)

            attempt_counter <- attempt_counter + 1
            }


        #Hess <- numDeriv::hessian(func = glmm_MLE_obj$fn, x = glmmfit$par)
        #SD = sdreport( Obj, hessian.fixed=Hess )
        glmmfit_sdreport <- sdreport(glmm_MLE_obj)
        pt_estimates <- as.list(glmmfit_sdreport, what = "Estimate", report = TRUE)
        pt_estimates$random_effects_covariance <- .convert_to_covariance(sd = pt_estimates$sd, rho = pt_estimates$unconstrained_cor_params)

        out <- list(fit = glmmfit, sdreport = glmmfit_sdreport, point_estimates = pt_estimates)


        if(score_vectors) {
            num_clus <- length(unique(data$id))
            all_score_vectors <- NULL

            for(i in 1:num_clus) {
                maketidbits_singleclus <- .construct_tidbits(response_type = response_type,
                                                             formula_X = formula_X[[select_response]],
                                                             formula_Z = formula_Z[[select_response]],
                                                             data = data[data$id == i,],
                                                             y = y[data$id == i, select_response],
                                                             trial_size = trial_size,
                                                             max_levels = max(y[,select_response]),
                                                             num_ordinal_levels = max(y[,select_response]))

                glmm_singleclusterMLE_obj <- MakeADFun(data = maketidbits_singleclus$data,
                                                       parameters = maketidbits_singleclus$parameters,
                                                       random = maketidbits_singleclus$random,
                                                       DLL = maketidbits_singleclus$DLL,
                                                       silent = TRUE)

                all_score_vectors <- rbind(all_score_vectors, glmm_singleclusterMLE_obj$gr(x = glmmfit$par))
                }

            rownames(all_score_vectors) <- 1:num_clus
            colnames(all_score_vectors) <- paste0(colnames(y)[select_response], "_", names(glmmfit$par))
            out$score_vectors <- all_score_vectors
            }

        return(out)
        }


    ##---------------------
    ## Final output
    ##---------------------
    allfits <- foreach(l = 1:num_resp) %dopar% dofits(select_response = l)
    names(allfits) <- colnames(y)

    allfits$call <- match.call()
    allfits$formula_X <- formula_X
    allfits$formula_Z <- formula_Z
    allfits$data <- data
    allfits$response_type <- response_type
    allfits$num_resp <- num_resp
    class(allfits) <- "stackedGLMM"
    return(allfits)
    }




#' @title Hidden function to construct data, parameter, and random effect arguments needed to pass to TMB for fitting a GLMM
#' @noRd
#' @noMd

.construct_tidbits <- function(response_type = "gaussian", y, formula_X, formula_Z, data,
    start_params = NULL, trial_size = 1, max_levels = NULL, num_ordinal_levels = NULL,
    seed = NULL) {

    docompile <- FALSE

    MM_fixef <- as.matrix(model.matrix(formula_X, data = data))
    MM_ranef <- as.matrix(model.matrix(formula_Z, data = data))
    y <- as.vector(y)
    if(response_type == "ordinal") {
        if(!is.numeric(y))
            stop("For ordinal response types, y must take on numeric values. It is assumed that each level of the ordinal variable has been observed at least once in y.")
        }

    response_type <- match.arg(response_type, choices = c("gaussian", "poisson", "binomial", "ordinal", "tweedie"))
    if(!(response_type %in% c("gaussian", "poisson", "binomial", "ordinal", "tweedie")))
        stop("Specified family not permitted. Sorry!")

    ##---------------------
    ## Compile TMB C++ script -- This is automatically fixed to FALSE since the compilation should have occurred outside of this.
    ##---------------------
    if(docompile) {
        origwd <- getwd()
        filename <- paste0("glmmMLE_", response_type, ".cpp")
        setwd("cppfiles")
        TMB::compile(filename) # flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign"
        dyn.load(TMB::dynlib(paste0("glmmMLE_",response_type)))
        setwd(origwd)
        }

	if(!is.null(seed))
		set.seed(seed)


    ##---------------------
    ## Declare data
    ##---------------------
    data$id <- as.numeric(factor(data$id))
    tmb_data <- list(y = y, X = MM_fixef, Z = MM_ranef, clusid = data$id-1, trial_size = trial_size, num_clus = length(unique(data$id)))
	if(response_type == "ordinal") {
        if(is.null(num_ordinal_levels))
            tmb_data$num_ordinal_levels <- length(unique(y))
        if(!is.null(num_ordinal_levels))
            tmb_data$num_ordinal_levels <- num_ordinal_levels
        if(is.null(max_levels))
            tmb_data$max_levels <- max(y)
        if(!is.null(max_levels))
            tmb_data$max_levels <- max_levels
	    }

    ##---------------------
    ## Declare parameters and constraints
    ##---------------------
    if(is.null(start_params)) {
        tmb_parameters <- list(fixed_effects = rnorm(ncol(MM_fixef)), random_effects = matrix(0, nrow = tmb_data$num_clus, ncol = ncol(MM_ranef)))
        if(response_type %in% c("gaussian", "tweedie"))
            tmb_parameters$phi <-
        if(response_type %in% c("tweedie"))
        tmb_parameters$powerparam <- 1.6
        if(response_type == "ordinal")
            tmb_parameters$difford <- runif(tmb_data$num_ordinal_levels-2)
        tmb_parameters$sd <- rep(1, ncol(MM_ranef))
        tmb_parameters$unconstrained_cor_params <- rep(0, ncol(MM_ranef)*(ncol(MM_ranef)-1)/2)
        }
    if(!is.null(start_params)) {
        tmb_parameters <- list(fixed_effects = start_params$fixed_effects, random_effects = start_params$random_effects)
        if(response_type %in% c("gaussian", "tweedie"))
            tmb_parameters$phi <- start_params$phi
        if(response_type %in% c("tweedie"))
        tmb_parameters$powerparam <- 1.6
        if(response_type == "ordinal")
            tmb_parameters$difford <- runif(tmb_data$num_ordinal_levels-2)
        tmb_parameters$sd <- start_params$sd
        tmb_parameters$unconstrained_cor_params <- start_params$unconstrained_cor_params
        }

    lower_limits <- rep(-Inf, length(unlist(tmb_parameters)))
    names(lower_limits) <- names(unlist(tmb_parameters))
    lower_limits[grep("phi", names(lower_limits))] <- 1e-3
    lower_limits[grep("powerparam", names(lower_limits))] <- 1+1e-3
    lower_limits[grep("sd", names(lower_limits))] <- 1e-3
    lower_limits[grep("unconstrained_cor_params", names(lower_limits))] <- -1e3
    if(response_type == "ordinal")
        lower_limits[grep("difford", names(lower_limits))] <- 1e-3

    upper_limits <- rep(Inf, length(unlist(tmb_parameters)))
    names(upper_limits) <- names(unlist(tmb_parameters))
    upper_limits[grep("powerparam", names(upper_limits))] <- 2-1e-3
    upper_limits[grep("unconstrained_cor_params", names(upper_limits))] <- 1e3

    ##---------------------
    ## What is the random bit?
    ##---------------------
    tmb_random <- c("random_effects")

	set.seed(NULL)
	return(list(parameters = tmb_parameters,
	            data = tmb_data,
	            random = tmb_random,
	            lower_limits = lower_limits,
	            upper_limits = upper_limits,
	            DLL = paste0("glmmMLE_",response_type)))
    }


