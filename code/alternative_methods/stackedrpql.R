#' @title Fit stacked univariate penalized generalized linear mixed models (GLMMs) using the [rpql::rpql()] package
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' This function applies [rpql::rpql()] in parallel across each response of the multivariate longitudinal data, where the adaptive LASSO penalty is applied to perform fixed effects selection only and the second hybrid information criterion of the package is used to select the tuning parameter. Please see [rpql::rpql()], [rpql::rpqlseq()], and the corresponding (CRAN website)[https://cran.r-project.org/web/packages/rpql/index.html] for more details and associated references. 
#' 
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data frame. The data **must** also be formatted such that the column names are \code{x1, x2, x3, ...}
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial".
#' @param num_fixef This argument determines the fixed effects component of the model. Effectively, the formula for the fixed effects component of the GLMM is determined as \code{~ x1 + x2 + x3 + ... + x_"num_fixef"}, where we note the \code{data} argument required the column names to be \code{x1, x2, x3, ...}.
#' @param formula_Z One-sided formula for the random effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{Sigma} argument.
#' @param lambda_seq A tuning parameter sequence. Note a default sequence is used based on a set of 100 equally spaced tuning parameter values on the log-scale, although it is highly recommended that this be tested.
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{doParallel::detectCores()-2}.
#' 
#' @return A list with the following:
#' \item{respxxx:}{The first \code{ncol(y)} elements contains the estimated penalized GLMMs corresponding to each response. Each element has the same structure as the same structure as the output from usage of the [rpql::rpql()] function. }
#' \item{call:}{The function call.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>


library(rpql)
library(lme4)


stacked_rpql <- function(y, data, response_type = "gaussian", num_fixef, formula_Z,
                         lambda_seq = lseq(1e-6,1,length=100), num_cores = NULL) {
    
    ##---------------------
    ## Checks and balances
    ##---------------------
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in parallel to fit stacked rpql GLMMs.")
    
    y <- as.matrix(y)
    num_resp <- ncol(y)
    
    response_type <- match.arg(response_type, choices = c("gaussian", "binomial", "poisson"))
    if(!(response_type %in% c("gaussian", "binomial", "poisson"))) 
        stop("Current response type not permitted. Sorry!")

    if(response_type == "gaussian")
       get_family <- gaussian(link = "identity")
    if(response_type == "poisson")
       get_family <- poisson(link = "log")
    if(response_type == "binomial")
       get_family <- binomial(link = "logit")
    
    ##---------------------
    ## Set up and fit the stacked GLMMs 
    ##---------------------
    dofits <- function(select_response) {
       cw_dat <- data.frame(resp = y[,select_response], data) 
       cw_dat$id <- as.factor(cw_dat$id)
       make_fixform <- as.formula(paste0("resp ~ ", paste0("x", 1:num_fixef, collapse  = " + "), " + ", formula_Z, collapse = " + "))
       
       fit_satlme4 <- lme4::glmer(make_fixform, data = cw_dat, family = response_type) 
       fit_sat <- build.start.fit(fit_satlme4, gamma = 2)
       fit_sat$pen.weights$random$id <- fit_sat$pen.weights$random$id*0 # Do not penalized random effects
        
        
       XMM <- unname(model.matrix(fit_satlme4)) 
       ZMM <- getME(fit_satlme4,"mmList"); names(ZMM) <- "cluster"
       
       fit_rpql <- try(rpqlseq(y = cw_dat$resp,
                           X = XMM,
                           Z = ZMM,
                           id = list(cluster = as.numeric(cw_dat$id)),
                           family = get_family,
                           lambda = lambda_seq,
                           pen.type = "adl",
                           pen.weights = fit_sat$pen.weights,
                           start = fit_sat), silent = TRUE)
       
       if(inherits(fit_rpql, "try-error"))
          final_fit <- NA
       if(!inherits(fit_rpql, "try-error"))
          final_fit <- fit_rpql$best.fit[[5]]
       
       return(final_fit)
       }
    
    allfits <- foreach(k = 1:num_resp) %dopar% dofits(select_response = k) 
    
    ##---------------------
    ## Final output
    ##---------------------
    names(allfits) <- colnames(y)
    
    allfits$call <- match.call()
    return(allfits)
    }


