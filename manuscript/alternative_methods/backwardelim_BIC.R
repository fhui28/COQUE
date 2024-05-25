#' @title Fit stacked univariate generalized linear mixed models (GLMMs) with fixed effects backward elimination using BIC.
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' This function applies [glmmTMB::glmmTMB()] or some variation in parallel across each response of the multivariate longitudinal data, where backward elimination is applied to the fixed effects via the Bayesian information criterion (BIC). Please see [glmmTMB::glmmTMB()] and [ordinal::clmm()] for more details and associated references. 
#' 
#' *Note for ordinal responses, it is recommend you use the \code{stacked_backwardelim_ordinalTMB} function below instead, as it is much faster.*
#' 
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data frame. The data **must** also be formatted such that the column names are \code{x1, x2, x3, ...}
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial", "ordinal", "tweedie". For ordinal responses though, it is recommend you use the \code{stacked_backwardelim_ordinalTMB} function below instead, as it is much faster.
#' @param num_fixef This argument determines the fixed effects component of the model. Effectively, the formula for the fixed effects component of the GLMM is determined as \code{~ x1 + x2 + x3 + ... + x_"num_fixef"}, where we note the \code{data} argument required the column names to be \code{x1, x2, x3, ...}.
#' @param formula_Z One-sided formula for the random effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{Sigma} argument.
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{doParallel::detectCores()-2}.
#' 
#' @return A list with the following:
#' \item{respxxx:}{The first \code{ncol(y)} elements contains the estimated penalized GLMMs corresponding to each response. Each element is itself a list containing only the element \code{coefficients}, which is the estimated coefficients using BIC.}
#' \item{call:}{The function call.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>



library(glmmTMB)
library(ordinal)

#' Note that num_fixed should exclude the intercept.
#' Note formula_Z must be compatible with the part used in glmmTMB and clmm e.g., (x1 + x2 + x3 + x4 | id)
stacked_backwardelim <- function(y, data, response_type = "gaussian", num_fixef, formula_Z, num_cores = NULL) {
    
    ##---------------------
    ## Checks and balances
    ##---------------------
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in parallel to fit stacked GLMMs.")
    
    y <- as.matrix(y)
    num_resp <- ncol(y)
    
    response_type <- match.arg(response_type, choices = c("gaussian", "binomial", "poisson", "ordinal", "tweedie"))
    if(!(response_type %in% c("gaussian", "binomial", "poisson", "ordinal", "tweedie"))) 
        stop("Current response type not permitted. Sorry!")

        
    ##---------------------
    ## Set up and fit the stacked GLMMs 
    ##---------------------
    dofits <- function(select_response) {
        cw_dat <- data.frame(resp = y[,select_response], data) 
        cw_dat$id <- as.factor(cw_dat$id)
        if(response_type == "ordinal")
            cw_dat$resp <- as.ordered(cw_dat$resp)
        
        make_fixform <- as.formula(paste0("resp ~ ", paste0("x", 1:num_fixef, collapse  = " + "), " + ", formula_Z, collapse = " + "))
        if(response_type != "ordinal")
            full_fit <- new_fit <- glmmTMB(make_fixform, data = cw_dat, family = response_type) 
        if(response_type == "ordinal")
            full_fit <- new_fit <- clmm(make_fixform, data = cw_dat) 
        best_BIC <- Inf
        useless_x <- NULL
        full_fit_BIC <- BIC(full_fit)
        if(!is.finite(full_fit_BIC))
            return(list(coefficients = fixef(full_fit)[[1]]))
        
        
        while(best_BIC > BIC(new_fit)) {
            best_BIC <- BIC(new_fit)
            dodrop <- try( drop1(new_fit, k = log(nrow(cw_dat))), silent = TRUE)
            
            if(inherits(dodrop, "try-error")) 
                break;
            if(best_BIC <= min(dodrop$AIC, na.rm = TRUE)) 
                break;
            
            useless_x <- c(useless_x, as.numeric(strsplit(rownames(dodrop)[which.min(dodrop$AIC)],"x")[[1]][2]))
            make_fixform <- as.formula(paste0("resp ~ ", paste0("x", (1:num_fixef)[-useless_x], collapse  = " + "), " + (x1 + x2 + x3 + x4|id)", collapse = "+"))
            if(response_type != "ordinal")
                new_fit <- glmmTMB(make_fixform, data = cw_dat, family = response_type) 
            if(response_type == "ordinal")
                new_fit <- clmm(make_fixform, data = cw_dat, family = response_type) 
            new_fit_BIC <- BIC(new_fit)
            if(!is.finite(new_fit_BIC))
                break;                
            }
        
        
        backelim_fit <- list(coefficients = numeric(length(fixef(full_fit)[[1]])))
        backelim_fit$coefficients[1] <- fixef(new_fit)[[1]][1]
        index_coefs <- as.numeric(sapply(strsplit(names(fixef(new_fit)[[1]][-1]),"x"), function(x) x[2]))+1
        backelim_fit$coefficients[index_coefs] <- fixef(new_fit)[[1]][-1]
        
        return(backelim_fit)
        }
    
    allfits <- foreach(k = 1:num_resp) %dopar% dofits(select_response = k) 
    
    
    ##---------------------
    ## Final output
    ##---------------------
    names(allfits) <- colnames(y)
    
    allfits$call <- match.call()
    return(allfits)
    }



#' @title Fit stacked univariate cumulative logit linear mixed models (GLMMs) with fixed effects backward elimination using BIC.
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' This function is preferred over the above function, which uses [ordinal::clmm()], as the fitting and hence backward elimination process using TMB is much faster.
#' 
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data frame. 
#' @param formula_X One-sided formula for the fixed effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{fixed_effects} argument.
#' @param formula_Z One-sided formula for the random effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{Sigma} argument.
#' @param TMB_directory The directory where the TMB C++ files are located.
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{doParallel::detectCores()-2}.
#' @param control A list controlling the some of the tuning parameters in the optimization process for fitting the stacked GLMMs, noting (nlminb)[stats::nlminb()] is used to perform the optimization. These include: \code{iter.max}, which is the maximum number of iterations allowed; \code{eval.max} which is the maximum number of evaluations of the objective function allowed; \code{trace}, which can be set to a positive value so that the value of the objective function and the parameters is printed every trace'th iteration; \code{attempts}, which is the number of attempts are made to fit each stacked univariate GLMM (if \code{nlminb} fails on the first try).
#'
#' @return A list with the following:
#' \item{respxxx:}{The first \code{ncol(y)} elements contains the estimated penalized GLMMs corresponding to each response. Each element is itself a list containing the two elements: \code{coefficients}, which is the estimated coefficients using BIC, and \code{point_estimates}, which a list of estimated parameters and random effects in a tidier format.}
#' \item{call:}{The function call.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>


stacked_backwardelim_ordinalTMB <- function(y, data, formula_X, formula_Z, TMB_directory = NULL, num_cores = NULL,
                                            control = list(iter.max = 1000, eval.max = 1000, trace = 0, attempts = 20)) {
    
    ##---------------------
    ## Checks and balances
    ##---------------------
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in parallel to fit stacked GLMMs.")
    
    y <- as.matrix(y)
    num_resp <- ncol(y)

    if(!is.list(formula_X))
        formula_X <- replicate(num_resp, formula_X, simplify = FALSE)
    if(!is.list(formula_Z))
        formula_Z <- replicate(num_resp, formula_Z, simplify = FALSE)
    
    
    ##---------------------
    ## Compile TMB C++ script -- doing it outside of the external functions below to save time and potential conflicts
    ##---------------------
    origwd <- getwd()
    filename <- paste0("glmmMLE_ordinal.cpp")
    setwd(here(TMB_directory))
    TMB::compile(filename, flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")
    dyn.load(TMB::dynlib(here(TMB_directory, "glmmMLE_ordinal")))
    setwd(origwd)
    
    
    ##---------------------
    ## Set up functions to fit a single ordinal GLMM, and drop1 
    ##---------------------
    single_ordinal_fit <- function(cw_y, cw_formula_X, cw_formula_Z, data) {
        maketidbits <- .construct_tidbits(response_type = "ordinal", 
                                          formula_X = cw_formula_X, 
                                          formula_Z = cw_formula_Z, 
                                          data = data, 
                                          y = cw_y, 
                                          trial_size = 1,
                                          docompile = FALSE
                                          )
        
        glmm_MLE_obj <- MakeADFun(data = maketidbits$data, 
                                  parameters = maketidbits$parameters,
                                  random = maketidbits$random, 
                                  DLL = maketidbits$DLL, 
                                  silent = TRUE)
        
        glmmfit <- try(nlminb(start = glmm_MLE_obj$par,
                              objective = glmm_MLE_obj$fn,
                              gradient = glmm_MLE_obj$gr,
                              lower = maketidbits$lower_limits,
                              control = control),
                       silent = TRUE)
        
        attempt_counter <- 1
        while(inherits(glmmfit, "try-error") & attempt_counter < control$attempts) {
            maketidbits <- .construct_tidbits(response_type = "ordinal", 
                                              formula_X = cw_formula_X, 
                                              formula_Z = cw_formula_Z, 
                                              data = data, 
                                              y = cw_y, 
                                              trial_size = trial_size,
                                              docompile = FALSE)
            
            glmm_MLE_obj <- MakeADFun(data = maketidbits$data, 
                                      parameters = maketidbits$parameters,
                                      random = maketidbits$random, 
                                      DLL = maketidbits$DLL, 
                                      silent = TRUE)
            
            glmmfit <- try(nlminb(start = glmm_MLE_obj$par,
                                  objective = glmm_MLE_obj$fn,
                                  gradient = glmm_MLE_obj$gr,
                                  lower = maketidbits$lower_limits,
                                  control = control),
                           silent = TRUE)
            
            attempt_counter <- attempt_counter + 1
            }
        
        if(!inherits(glmmfit, "try-error")) {
            glmmfit_sdreport <- sdreport(glmm_MLE_obj, ignore.parm.uncertainty = TRUE) 
            pt_estimates <- as.list(glmmfit_sdreport, what = "Estimate", report = TRUE)
            pt_estimates$random_effects_covariance <- .convert_to_covariance(sd = pt_estimates$sd, rho = pt_estimates$unconstrained_cor_params)
            cw_BIC <- 2 * glmmfit$objective + log(length(cw_y)) * length(glmmfit$par)
            
            out <- list(fit = glmmfit, point_estimates = pt_estimates, BIC = cw_BIC, #sdreport = glmmfit_sdreport, 
                        cw_formula_X = cw_formula_X, cw_formula_Z = cw_formula_Z)
            }
        
        if(inherits(glmmfit, "try-error")) {
            out <- list(BIC = Inf, cw_formula_X = cw_formula_X, cw_formula_Z = cw_formula_Z)
            }

        return(out)
        }
    
    drop1_ordinal_fit <- function(cw_fit, cw_y, data) {
        potentialDrops <- drop.scope(cw_fit$cw_formula_X, formula("~ 1"))
        if(length(potentialDrops) == 0)
            return()
        
        res <- foreach(theDrop = potentialDrops) %do% #' Do *not* run in parallel the outer stacked function itself is already in parallel
            single_ordinal_fit(cw_y = cw_y, 
                               cw_formula_X = update.formula(cw_fit$cw_formula_X, as.formula(paste("~ . -", theDrop))), 
                               cw_formula_Z = cw_fit$cw_formula_Z,
                               data = data)
        res$BIC <- sapply(res, function(x) x$BIC)
        
        return(res)
        }
    
    
    ##---------------------
    ## Go go backward elimination!
    ##---------------------
    dofits <- function(select_response) {
        full_fit <- new_fit <- single_ordinal_fit(cw_y = y[,select_response], 
                                                  cw_formula_X = formula_X[[select_response]], 
                                                  cw_formula_Z = formula_Z[[select_response]],
                                                  data = data) 
        best_BIC <- Inf
        useless_x <- NULL
        full_fit_BIC <- full_fit$BIC
        if(!is.finite(full_fit_BIC))
            return(list(coefficients = full_fit$point_estimates$fixed_effects, point_estimates = full_fit$point_estimates))
        
        while(best_BIC > new_fit$BIC) {
            best_BIC <- new_fit$BIC
            dodrop <- try(drop1_ordinal_fit(new_fit, cw_y = y[,select_response], data = data), silent = TRUE)

            if(is.null(dodrop)) # Intercept-only model reached 
                break;
            if(inherits(dodrop, "try-error")) 
                break;
            if(best_BIC <= min(dodrop$BIC, na.rm = TRUE)) 
                break;
            
            new_fit <- dodrop[[which.min(dodrop$BIC)]]
            new_fit_BIC <- new_fit$BIC
            if(!is.finite(new_fit_BIC))
                break;                
            }
        
        
        backelim_fit <- list(coefficients = numeric(length(full_fit$point_estimates$fixed_effects)))
        names(backelim_fit$coefficients) <- colnames(model.matrix(full_fit$cw_formula_X, data = data))
        index_coefs <- match(colnames(model.matrix(new_fit$cw_formula_X, data = data)), names(backelim_fit$coefficients))
        backelim_fit$coefficients[index_coefs] <- new_fit$point_estimates$fixed_effects
        backelim_fit$point_estimates <- new_fit$point_estimates
        
        return(backelim_fit)
        }
    
    allfits <- foreach(k = 1:num_resp) %dopar% dofits(select_response = k) 
    

    ##---------------------
    ## Final output
    ##---------------------
    names(allfits) <- colnames(y)
    
    allfits$call <- match.call()
    return(allfits)
    }

