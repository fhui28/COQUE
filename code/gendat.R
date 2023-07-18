#' @title Simulate data from a independent cluster multivariate generalized linear mixed model (GLMM) 
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' @param formula_X Formula for the fixed effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{fixed_effects} argument.
#' @param formula_Z Formula for the random effects component of the model. All responses are assumed to be effected by the same formula; sparsity and clustering is instead controlled by the \code{Sigma} argument.
#' @param data A data frame from which containing relevant fixed and random effect covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data fame.
#' @param fixed_effects A matrix of fixed effect coefficients, where the number of rows is equal to the number of responses and the number columns must be equal to the number of columns of the model matrix implied by \code{formula_X}.
#' @param Sigma A random effects covariance matrix, where the number of rows and columns both must equal to the number of responses multiplied by the number of columns of the model matrix implied by \code{formula_Z}.
#' @param phi A vector of positive response-specific dispersion parameters, which is needed for Gaussian and Tweedie distributed responses.
#' @param power A vector of response-specific power parameters, which is needed for Tweedie distributed responses. All elements should be between 1 and 2.
#' @param trial_size A vector of response-specific trial sizes, which is needed for binomial distributed responses.
#' @param cutoffs A matrix of response-specific cutoffs, where the number of rows is equal to the number of responses, the number of columns determines the number of levels for the ordinal response, and the first column of the matrix must be a matrix of zeros for parameter identifiability reasons (since a fixed intercept is assumed to be included by default as part of \code{formula_X}). 
#' @param response_type The distribution of the responses for the model. Currently, it is assumed all responses follow the same distribution, which can be one of "gaussian", "poisson", "binomial", "ordinal", "tweedie".
#' @param max_count This argument can be used to try and control the size of the simulated responses. Basically, for any given response the function will try \code{attempts} times to (re)simulate until all the responses are below \code{max_count}. Currently this only appplies to Poisson and Tweedie distributed response. 
#' @param attempts The number of attempts used to try and re-simulate responses subject to \code{max_count}.
#' 
#' @return A list with the following elements: 
#' \item{y:}{The simulated data frame of multivariate longitudinal responses, where the number of columns equal to the number of responses plus one. The last column is just equal to \code{data$id}.}
#' \item{random_effects:}{A list of simulated random effects with length equal to the number of clusters. Each element in the list is a matrix with the number of rows equal to the number of responses and the number of columns equal to number of columns of the model matrix implied by \code{formula_Z}.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>




library(tweedie)

generatedata_mglmm <- function(formula_X, formula_Z, data, fixed_effects, 
                               Sigma, phi = NULL, power = NULL, trial_size = 1, cutoffs = NULL, response_type = "gaussian", 
                               max_count = Inf, attempts = 50) {

    num_resp <- nrow(fixed_effects)
    MM_fixef <- model.matrix(formula_X, data = data)
    MM_ranef <- model.matrix(formula_Z, data = data)
    if(ncol(MM_fixef) != ncol(fixed_effects)) 
        stop("The number of columns implied by the model matrix for formula_X does not match ncol(fixed_effects).")
    if(is.null(rownames(fixed_effects)))
        rownames(fixed_effects) <- paste0("resp", 1:num_resp)
    
    if(nrow(Sigma) != num_resp*ncol(MM_ranef)) 
        stop("Dimensions of Sigma and the number of columns implied by the model matrix for formula_Z do not match.")
    if(ncol(Sigma) != num_resp*ncol(MM_ranef)) 
        stop("Dimensions of Sigma and the number of columns implied by the model matrix for formula_Z do not match.")
    if(!is.null(cutoffs)) {
        if(nrow(cutoffs) != num_resp) 
            stop("Number of rows in cutoffs should be equal to number of rows in beta.")
        if(any(cutoffs[,1] != 0)) 
            stop("First column in cutoffs must be all equal to zero for identifiability reasons.")
        }
    
    response_type <- match.arg(response_type, choices = c("gaussian", "poisson", "binomial", "ordinal", "tweedie"))
    if(!(response_type %in% c("gaussian", "poisson", "binomial", "ordinal", "tweedie"))) 
        stop("Specified family not permitted. Sorry!")
    
    if(is.null(data$id))
        stop("data must contain an id column identifying the cluster id")
    data$id <- as.numeric(factor(data$id))
    

    ##---------------------
    ## Simulate random effects
    ##---------------------
    num_clus <- length(unique(data$id))
    true_ranef <- rmvnorm(num_clus, sigma = Sigma + diag(x = 1e-8, nrow = nrow(Sigma)))
    true_ranef[, which(diag(Sigma) == 0)] <- 0
    

    ##---------------------
    ## Simulate responses
    ##---------------------
    sim_resp <- matrix(NA, nrow = nrow(MM_fixef), ncol = num_resp)
    rownames(sim_resp) <- rownames(MM_fixef)
    colnames(sim_resp) <- rownames(fixed_effects)
    eta <- tcrossprod(MM_fixef, fixed_effects)
    for(i in 1:num_clus) {
        sel_cluster <- which(data$id == i)
        eta[sel_cluster,] <- eta[sel_cluster,] + tcrossprod(MM_ranef[sel_cluster,], matrix(true_ranef[i,], nrow = num_resp, ncol = ncol(MM_ranef), byrow = TRUE)) 
        for(k in 1:num_resp) {
            if(response_type == "gaussian") 
                sim_resp[sel_cluster, k] <- rnorm(length(sel_cluster), mean = eta[sel_cluster,k], sd = sqrt(phi[k]))
            if(response_type == "poisson") 
                sim_resp[sel_cluster, k] <- rpois(length(sel_cluster), lambda = exp(eta[sel_cluster,k]))
            if(response_type == "binomial") 
                sim_resp[sel_cluster, k] <- rbinom(length(sel_cluster), size = trial_size, prob = plogis(eta[sel_cluster,k]))
            if(response_type == "ordinal") 
                sim_resp[sel_cluster, k] <- .rorddata(cutoffs = cutoffs[k,], eta = eta[sel_cluster,k]) 
            if(response_type == "tweedie") 
               sim_resp[sel_cluster, k] <- rtweedie(length(sel_cluster), mu = exp(eta[sel_cluster,k]), phi = phi[k], power = power[k])
         }
        
        if(response_type %in% c("poisson", "tweedie")) {
            attempt_counter <- 0
            while(any(apply(sim_resp[sel_cluster,],2,max) > max_count) & attempt_counter <= attempts) {
                true_ranef[i,] <- as.vector(rmvnorm(1, sigma = Sigma + diag(x = 1e-8, nrow = nrow(Sigma))))
                true_ranef[i, which(diag(Sigma) == 0)] <- 0
                
                eta[sel_cluster,] <- tcrossprod(MM_fixef[sel_cluster,], fixed_effects) + tcrossprod(MM_ranef[sel_cluster,], matrix(true_ranef[i,], nrow = num_resp, ncol = ncol(MM_ranef), byrow = TRUE)) 
                
                for(k in 1:num_resp) {
                    if(response_type == "poisson") 
                        sim_resp[sel_cluster, k] <- rpois(length(sel_cluster), lambda = exp(eta[sel_cluster,k]))
                    if(response_type == "tweedie") 
                       sim_resp[sel_cluster, k] <- rtweedie(length(sel_cluster), mu = exp(eta[sel_cluster,k]), phi = phi[k], power = power[k])
                  }
                
                attempt_counter <- attempt_counter + 1
                }
            }
        
        
        
        }
    sim_resp <- as.data.frame(sim_resp)
    sim_resp$id <- data$id
    

    ##---------------------
    ## Finish up 
    ##---------------------
    true_ranef_list <- vector("list", num_clus) 
    for(i in 1:num_clus) {
        true_ranef_list[[i]] <- matrix(true_ranef[i,], nrow = num_resp, ncol = ncol(MM_ranef), byrow = TRUE)	
        rownames(true_ranef_list[[i]]) <- rownames(fixed_effects) 
        colnames(true_ranef_list[[i]]) <- colnames(MM_ranef) 
        }
    
    out <- list(y = sim_resp, random_effects = true_ranef_list)
    return(out)
    }



.rorddata <- function(cutoffs, eta, link = "logit") {
    out <- numeric(length(eta))
    Fs <- cbind(0, binomial(link = link)$linkinv(matrix(cutoffs, nrow = length(eta), ncol = length(cutoffs), byrow = TRUE) - eta),1)
    probs <- t(apply(Fs, 1, diff))
    for(i in 1:length(out)) 
        out[i] <- sample(1:ncol(probs), size = 1, prob = probs[i,])
    
    return(out)
    }


