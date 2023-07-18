#' @title Fit stacked univariate generalized linear mixed models (GLMMs)
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' This function can be used on its own, but in the context of this project is primarily used upstream of applying COQUE for simultaneous fixed effects selection and coefficient clustering. Fitting of the GLMMs is done via TMB, and optimizes the Laplace-approximated marginal log-likelihood function.
#'
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial", "ordinal", "tweedie".
#' @param formula_X One-sided formula for the fixed effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{fixed_effects} argument.
#' @param formula_Z One-sided formula for the random effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{Sigma} argument.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data fame.
#' @param trial_size A vector of response-specific trial sizes, which is needed for binomial distributed responses.
#' @param start_params Starting values can be supplied. This must take the form of a list with length equal to the number of columns in \code{y}. Each element in the list should then contain: \code{fixed_effects}, which is a vector starting values for the fixed effect coefficients; \code{random_effects}, which is a matrix of starting values for the random effect coefficients where the number of columns is equal to the number of columns of by the model matrix implied by \code{formula_Z}; \code{sd}, which is a vector of standard deviation parameters for the random effects covariance matrix; \code{unconstrained_cor_params}, which is a vector of parameters corresponding to the lower-triangular elements of the Cholesky decomposition of the random effects correlation matrix. We refer the reader to [https://kaskr.github.io/adcomp/_book/Densities.html] to see how covariance matrices are parametrized in TMB.
#' @param TMB_directory The directory where the TMB C++ files are located.
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{doParallel::detectCores()-2}.
#' @param score_vectors Should vectors of the score values for each cluster within each responses be returned. Defaults to \code{TRUE} and this generally should not be turned off.
#' @param control A list controlling the some of the tuning parameters in the optimization process for fitting the stacked GLMMs, noting (nlminb)[stats::nlminb()] is used to perform the optimization. These include: \code{iter.max}, which is the maximum number of iterations allowed; \code{eval.max} which is the maximum number of evaluations of the objective function allowed; \code{trace}, which can be set to a positive value so that the value of the objective function and the parameters is printed every trace'th iteration; \code{attempts}, which is the number of attempts are made to fit each stacked univariate GLMM (if \code{nlminb} fails on the first try).
#'
#'
#' @return A object for class "stackedGLMM" which contains the following:
#' \item{respxxx:}{The first \code{ncol(y)} elements contains the estimated stacked GLMMs corresponding to each response. Each element is itself a list containing: the fitted model output from TMB (\code{fit} and \code{sdreport}), a list of estimated parameters and random effects in a slightly tidier format (\code{point_estimates}), and a matrix of score vector values for each cluster (\code{socre_vectors}).}
#' \item{call:}{The function call.}
#' \item{formula_X/formula_Z/data/response_type/num_resp:}{Arugment information that can be safely ignored.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>


stackedGLMMs <- function(y, response_type, formula_X, formula_Z, data, trial_size = 1,
                         start_params = NULL, TMB_directory = NULL, num_cores = NULL, score_vectors = TRUE,
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
    origwd <- getwd()
    filename <- paste0("glmmMLE_", response_type, ".cpp")
    setwd(here(TMB_directory))
    TMB::compile(filename, flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")
    dyn.load(TMB::dynlib(here(TMB_directory, paste0("glmmMLE_",response_type))))
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
                                          trial_size = trial_size,
                                          docompile = FALSE
                                          )

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
                                              trial_size = trial_size,
                                              docompile = FALSE)

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
                                                             num_ordinal_levels = max(y[,select_response]),
                                                             docompile = FALSE)

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
    allfits <- foreach(k = 1:num_resp) %dopar% dofits(select_response = k)
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

.construct_tidbits <- function(response_type = "gaussian", y, formula_X, formula_Z, data,
    start_params = NULL, trial_size = 1, max_levels = NULL, num_ordinal_levels = NULL,
    docompile = FALSE, seed = NULL) {
    
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
    ## Compile TMB C++ script 
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


