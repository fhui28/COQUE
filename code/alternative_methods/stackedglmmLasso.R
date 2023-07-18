#' @title Fit stacked univariate penalized generalized linear mixed models (GLMMs) using the [glmmLasso::glmmLasso()] package
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' This function applies [glmmLasso::glmmLasso()] in parallel across each response of the multivariate longitudinal data, where the LASSO penalty is applied to perform fixed effects selection only and the Bayesian information criterion (BIC) of the package is used to select the tuning parameter. Please see [glmmLasso::glmmLasso()] and the corresponding (CRAN website)[https://cran.r-project.org/web/packages/glmmLasso/index.html] for more details and associated references. 
#' 
#' @param fix A two-sided formula object describing the fixed-effects component of the model, with \code{resp} specifically on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right. Please see the \code{fix} argument in [glmmLasso::glmmLasso()] for more details.
#' @param rnd A two-sided formula object describing the random-effects part of the model, with the grouping factor on the left of a \code{~} operator and the random terms, separated by \code{+} operators, on the right. Please see the \code{rnd} argument in [glmmLasso::glmmLasso()] for more details.
#' @param y A matrix of multivariate longitudinal responses, where the number of columns equal to the number of responses.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data frame. 
#' @param lambda_seq A tuning parameter sequence. 
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial", "ordinal", "tweedie".
#' @param num_cores The number of cores to use for parallelization of the stacked GLMM fitting. Defaults to the \code{doParallel::detectCores()-2}.
#' @param switch.NR Should the algorithm swith to a Newton-Raphson update step, when reasonable? Default is \code{FALSE} as per [glmmLasso::glmmLasso()].
#' @param tweedie_var_power If Tweedie distributed responses are assumed, then this argument must be supplied to inform a fixed power parameter for each response. 
#' @param steps The maximum number of iterations. Please see the \code{steps} argument in [glmmLasso::glmmLassoControl()] for more details.
#' 
#' @return A list with the following:
#' \item{respxxx:}{The first \code{ncol(y)} elements contains the estimated penalized GLMMs corresponding to each response. Each element is itself a list with the following elements: \code{final_lambda}, which is the value of the tuning parameter chosen by BIC; \code{final_fixed_effects}, which is the vector of estimated coefficients using BIC; \code{BIC}, which is the vector of BIC values along the regularization path; \code{regularization_path}, which the matrix of estimated coefficients along the regularization path. The number of columns in \code{regularization_path} should be equal to \code{length(lambda_seq)}.}
#' \item{call:}{The function call.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>, and the author/s of the [glmmLasso::glmmLasso()] package.


library(glmmLasso)
library(statmod)

stacked_glmmLasso <- function(fix, rnd, y, data, lambda_seq, 
                              response_type = "gaussian", num_cores = NULL, switch.NR = FALSE,
                              tweedie_var_power = NULL, steps = 1000) {
    
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
    
    response_type <- match.arg(response_type, choices = c("gaussian", "poisson", "binomial", "ordinal", "tweedie"))
    if(!(response_type %in% c("gaussian", "poisson", "binomial", "ordinal", "tweedie"))) 
        stop("Specified family not permitted. Sorry!")

    if(response_type == "gaussian")
        get_family <- gaussian(link = "identity")
    if(response_type == "poisson")
        get_family <- poisson(link = "log")
    if(response_type == "binomial")
        get_family <- binomial(link = "logit")
    if(response_type == "ordinal")
        get_family <- glmmLasso::cumulative()
    
    use_control <- glmmLassoControl()
    use_control$steps <- steps
    
    ##---------------------
    ## Set up and fit the stacked GLMMs -- function adapted from [https://davidabugaber.com/blog/f/find-the-optimal-mixed-model-for-your-data-with-glmmlasso]
    ##---------------------
    dofits <- function(select_response) {
        cw_dat <- data.frame(resp = y[,select_response], data)
        BIC_vec <- rep(Inf, length(lambda_seq))
        Coeff_ma <- NULL
    
        for(j in 1:length(BIC_vec)) {
            #print(paste("Iteration ", j, sep=""))
        
            if(response_type != "tweedie")
               glm1 <- try(glmmLasso(fix = fix,
                                     rnd = rnd,
                                     data = cw_dat,
                                     family = get_family,
                                     lambda = lambda_seq[j],
                                     switch.NR = switch.NR,
                                     final.re = TRUE, # This tends to perform better in general
                                     control = use_control),
                           silent = TRUE)
            if(response_type == "tweedie")
               glm1 <- try(glmmLasso(fix = fix,
                                     rnd = rnd,
                                     data = cw_dat,
                                     family = tweedie(var.power = tweedie_var_power[select_response], link.power = 0),
                                     lambda = lambda_seq[j],
                                     switch.NR = switch.NR,
                                     final.re = TRUE, # This tends to perform better in general
                                     control = use_control),
                           silent = TRUE)
            
        if(!inherits(glm1, "try-error")) {
            BIC_vec[j] <- glm1$bic
            Coeff_ma <- cbind(Coeff_ma, glm1$coefficients)
            }
        if(inherits(glm1, "try-error")) {
            BIC_vec[j] <- Inf
            Coeff_ma <- cbind(Coeff_ma, rep(NA, nrow(Coeff_ma)))
            }
        }
    
        
        out <- list(final_lambda = lambda_seq[which.min(BIC_vec)], 
                    final_fixed_effects = Coeff_ma[,which.min(BIC_vec)],
                    BIC = BIC_vec,
                    regularization_path = Coeff_ma
                    )
                
        return(out)
        }
        
    
    ##---------------------
    ## Final output
    ##---------------------
    allfits <- foreach(k = 1:num_resp) %dopar% dofits(select_response = k)
    names(allfits) <- colnames(y)
    
    allfits$call <- match.call()
    return(allfits)
    }
    

# logLik.glmmLasso<-function(y,yhelp,mu,ranef.logLik=NULL,family,penal=FALSE, K = NULL, phi = 1, power = NULL) 
# {
#    fam <- family$family
#    
#    if(fam=="poisson")
#       loglik <- sum(y*log(mu)-mu-log(factorial(y)))
#    
#    if(fam=="binomial")
#       loglik <- sum(log(mu[y==1]))+sum(log((1-mu)[y==0]))
#    
#    if(fam=="Tweedie")
#       loglik <- sum(dtweedie(y, mu = mu, phi = phi, power = power, log = TRUE))
#    
#    if(fam=="gaussian")
#       loglik <- length(y)*-0.5*(log(2)+log(pi)+log(phi)) + sum(-0.5*(y-mu)^2/phi)
#    
#    if(fam=="acat")
#    {
#       mu_cat <- matrix(mu, byrow = TRUE, ncol = K)
#       mu_cat <- cbind(mu_cat,1-rowSums(mu_cat))
#       # yhelp <- matrix(y, byrow = TRUE, ncol = K)
#       # yhelp <- cbind(yhelp,1-rowSums(yhelp))
#       loglik <- sum(yhelp*log(mu_cat))
#    }
#    
#    if(fam=="cumulative")
#    {
#       mu_cat <- matrix(mu, byrow = TRUE, ncol = K)  
#       mu_help <- mu_cat
#       
#       
#       for (i in 2:K){
#          mu_cat[,i] <- mu_help[,i]-mu_help[,i-1]
#       }
#       
#       mu_cat <- cbind(mu_cat,1-mu_help[,K])
#       
#       loglik <- sum(yhelp*log(mu_cat))
#       #loglik 
#    }
#    
#    if(penal)
#       loglik <- loglik + ranef.logLik 
#    
#    return(loglik)
# }
