#' @title Extract values from a penalized composite quadratic estimator (COQUE) fit.
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' @param object An object of class \code{coque}.
#' @param lambda The tuning parameter at which estimates are estimated. This can either be a string taking values e.g., \code{"BIC"} corresponding to taking the estimates based on minimizing of these information criterion, or a numeric value corresponding to specific tuning parameter value. Similar to (coef.glmnet)[glmnet::coef.glmnet()], in the case of the latter it uses linear interpolation to construct parameter estimates for values of \code{lambda} that do not coincide with those used in the fitting algorithm. 
#' @param hybrid Should a hybrid estimator be calculated also? If \code{TRUE}, then aside from the penalized estimator, a second hybrid estimator is obtained by fitting the composite quadratic approximation subject to the zero constraints implied by the penalized estimator i.e., parameters shrunk to zero are constrained to be zero.
#' 
#' @return A list with the following elements:
#' \item{fixed_effects_matrix}{The matrix of estimated fixed effect coefficients. Due to the use of the broken adaptive ridge fusion (BARF) penalty, some of the estimated coefficients will be shrunk to zero, but coefficients within a covariate many also take the same value i.e., are fused together.}
#' \item{random_effects_covariance}{A list of random effects covariance matrices for each responses.}
#' \item{psihat}{The full vector of parameter estimates at the specified tuning parameter value.}
#' \item{lambda}{The numerical value of the tuning parameter corresponding to the extracted estimates.}
#' \item{hybrid_fixed_effects_matrix}{If \code{hybrid = TRUE}, then the matrix of hybrid estimator of the fixed effect coefficients. Due to the use of the broken adaptive ridge fusion (BARF) penalty, some of the estimated coefficients will be shrunk to zero, but coefficients within a covariate many also take the same value i.e., are fused together.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>


estimates <- function(object, ...) {
    UseMethod("estimates")
}

estimates.coque <- function(object, lambda, hybrid = FALSE) {
    if(class(object) != "coque")
        stop("object must be of class 'coque'.")
    
    ##-------------------
    ## Extract psihat
    ##-------------------
    if(!is.numeric(lambda)) {
        lambda <- match.arg(lambda, choices = c("AIC", "BIC", "EBIC", "EBIC2", "ERIC", "ERIC2"))
        sel_lambda <- which.min(object$ICs[,lambda])
        
        psihat <- object$psi_path[,sel_lambda]
        final_lambda <- object$lambda[sel_lambda]
    }    
    
    if(is.numeric(lambda)) {
        lambda_list <- .lambda_interp(lambda_seq = object$lambda, s = lambda)
        
        psihat <- object$psi_path[, lambda_list$left, drop = FALSE] %*% Diagonal(x = lambda_list$frac) + object$psi_path[, lambda_list$right, drop = FALSE] %*% Diagonal(x = 1 - lambda_list$frac)
        psihat <- as.vector(psihat)
        names(psihat) <- rownames(object$psi_path)
        final_lambda <- lambda
    }
    
    
    ##-------------------
    ## Format
    ##-------------------
    fixed_effects_matrix <- matrix(psihat[grep("fixed_effects",names(psihat))], nrow = object$num_resp, byrow = TRUE)
    rownames(fixed_effects_matrix) <- object$response_names
    colnames(fixed_effects_matrix) <- object$fixed_effects_names
    
    remaining_parameters <- psihat[-grep("fixed_effects",names(psihat))]
    all_sd <- remaining_parameters[grep("sd", names(remaining_parameters))]
    all_rho <- remaining_parameters[grep("unconstrained_cor_params", names(remaining_parameters))]
    num_ranef <- length(all_sd)/object$num_resp
    
    random_effects_covariance <- vector("list", object$num_resp)
    for(k in 1:object$num_resp) {
        random_effects_covariance[[k]] <- .convert_to_covariance(
            sd = all_sd[(k*num_ranef - num_ranef + 1):(k*num_ranef)],
            rho = all_rho[(k*((num_ranef*(num_ranef-1)/2)) - (num_ranef*(num_ranef-1)/2) + 1):(k*((num_ranef*(num_ranef-1)/2)))]
        )
    }
    
    names(random_effects_covariance) <- object$response_names
    rm(all_sd, all_rho, num_ranef)
    
    
    ##-------------------
    ## Compute hybrid estimator
    ##-------------------
    if(!hybrid) {
        hybrid_fixed_effects_matrix <- NULL
    }
    if(hybrid) {
        getpenalized_estimates <- as.vector(object$additional_information$bigRSmatrix %*% psihat) 
        A <- Matrix(0, nrow = length(getpenalized_estimates), ncol = sum(getpenalized_estimates == 0))
        for(l0 in 1:ncol(A))
            A[which(getpenalized_estimates == 0)[l0], l0] <- 1
        
        get_hybrid_estimate <- solve.QP(Dmat = Matrix::nearPD(object$additional_information$reparametrized_sandwich_matrix)$mat, 
                                        dvec = Matrix::nearPD(object$additional_information$reparametrized_sandwich_matrix)$mat %*% object$additional_information$reparametrized_psihat, 
                                        Amat = A, 
                                        meq = ncol(A), 
                                        factorized = FALSE)
        
        # get_hybrid_estimate$unconstrained.solution - penalized_fit$additional_information$reparametrized_psihat
        get_hybrid_estimate$solution[which(getpenalized_estimates == 0)] <- 0
        get_hybrid_estimate$solution <- as.vector(solve(object$additional_information$bigRSmatrix) %*% get_hybrid_estimate$solution)
        names(get_hybrid_estimate$solution) <- rownames(object$psi_path) 
        
        hybrid_fixed_effects_matrix <- matrix(get_hybrid_estimate$solution[grep("fixed_effects",names(get_hybrid_estimate$solution))], 
                                              nrow = object$num_resp, byrow = TRUE)
        rownames(hybrid_fixed_effects_matrix) <- object$response_names
        colnames(hybrid_fixed_effects_matrix) <- object$fixed_effects_names
    }
    
    
    out <- list(fixed_effects_matrix = fixed_effects_matrix, random_effects_covariance = random_effects_covariance, 
                psihat = psihat, lambda = final_lambda, hybrid_fixed_effects_matrix = hybrid_fixed_effects_matrix)
    return(out)
}

