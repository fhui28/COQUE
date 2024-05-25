#' @title Fit stacked univariate penalized generalized linear mixed models (GLMMs) using the [glmmPen::glmmPen()] package
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' This function applies [glmmPen::glmmPen()] in parallel across each response of the multivariate longitudinal data, where the minimax concave plus (MCP) penalty is applied to perform fixed effects selection only and a hybrid Bayesian information criterion of the package is used to select the tuning parameter. Please see [glmmPen::glmmPen()] and the corresponding (CRAN website)[https://cran.r-project.org/web/packages/glmmPen/index.html] for more details and associated references. 
#' 
#' @formula A two-sided linear formula object describing both the fixed effects and random effects components of the GLMM, with \code{resp} specifically on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right. Random-effects terms are distinguished by vertical bars ("|") separating expression for design matrices from the grouping factor. Please see the \code{formula} argument in [glmmPen::glmmPen()] for more information.
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data frame.
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial".
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{doParallel::detectCores()-2}.
#' 
#' @return A list with the following:
#' \item{respxxx:}{The first \code{ncol(y)} elements contains the estimated penalized GLMMs corresponding to each response. Each element has the same structure as the same structure as the output from applying \code{fixef()} function to a [glmmPen::glmmPen()] fit. }
#' \item{call:}{The function call.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>, and the authors of [glmmPen::glmmPen()].


library(glmmPen)
library(ncvreg)

stacked_glmmPen <- function(formula, y, data, response_type = "gaussian", num_cores = NULL) {
    
    ##---------------------
    ## Checks and balances
    ##---------------------
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in parallel to fit stacked GLMMs.")
    
    data$id <- factor(data$id)
    y <- as.matrix(y)
    num_resp <- ncol(y)
    
    response_type <- match.arg(response_type, choices = c("gaussian", "poisson", "binomial"))
    if(!(response_type %in% c("gaussian", "poisson", "binomial"))) 
        stop("Specified family not permitted. Sorry!")

    if(response_type == "gaussian")
        get_family <- gaussian()
    if(response_type == "poisson")
        get_family <- poisson()
    if(response_type == "binomial")
        get_family <- binomial()

    ##---------------------
    ## Set up and fit the stacked GLMMs -- function adapted from [https://davidabugaber.com/blog/f/find-the-optimal-mixed-model-for-your-data-with-glmmlasso]
    ##---------------------
    dofits <- function(select_response) {
        cw_dat <- data.frame(resp = y[,select_response], data)

       glm1 <- try(glmmPen(formula = formula,
                           data = cw_dat,
                           family = get_family,
                           penalty = "MCP",
                           tuning_options = selectControl(lambda1_seq = 0, BIC_option = c("BICh", "BIC", "BICNgrp")),
                           progress = FALSE),
                   silent = TRUE)

        if(!inherits(glm1, "try-error")) {
            out <- fixef(glm1)
            }

        if(inherits(glm1, "try-error")) {
            out <- NA
            }

        return(out)
        }
        
    
    ##---------------------
    ## Final output
    ##---------------------
    allfits <- foreach(k = 1:num_resp) %dopar% dofits(select_response = k)
    names(allfits) <- colnames(y)
    gc()
    
    allfits$call <- match.call()
    return(allfits)
    }
    

