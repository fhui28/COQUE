
#' @title A hidden function to use linear interpolation to construct parameter estimates for values of the tuning parameter that do not coincide with those used in the fitting algorithm.
#'
#' @param lambda_seq is the index sequence that is produced by the model. This *must* be a strictly decreasing sequence.
#' @param s is the new vector at which evaluations are required.
#'
#' @return A list containing the left and right indices to take, and the corresponding fraction. The new parameter values are interpolated between the two using the fraction.
#'
#' @note Modified from and acknowledgments goes to the \code{glmnet} team!
#' @note Note: lambda_seq decreases. you take: \code{sfrac*left+(1-sfrac*right)}
#' @noRd
#' @noMd

.lambda_interp <- function(lambda_seq, s) {

    if(length(lambda_seq) == 1){
        nums <- length(s)
        left <- rep(1,nums)
        right <- left
        sfrac <-rep(1,nums)
        }

    else {
        k <- length(lambda_seq)
        sfrac <- (lambda_seq[1] - s)/(lambda_seq[1] - lambda_seq[k])
        lambda_seq <- (lambda_seq[1] - lambda_seq)/(lambda_seq[1] - lambda_seq[k])

        sfrac[sfrac < min(lambda_seq)] <- min(lambda_seq)
        sfrac[sfrac > max(lambda_seq)] <- max(lambda_seq)

        coord <- stats::approx(lambda_seq, seq(lambda_seq), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)

        sfrac <- (sfrac - lambda_seq[right])/(lambda_seq[left] - lambda_seq[right])
        sfrac[left == right] <- 1
        sfrac[abs(lambda_seq[left] - lambda_seq[right]) < .Machine$double.eps] <- 1
        }

    return(list(left = left,right = right,frac = sfrac))
    }


#' @title Hidden internal function to calculate the covariance matrix given the standard deviations and the lower triangular elements in the Cholesky decomposition of the correlation matrix (the diagonal matrix of standard deviations) and the lower triangular elements of R (a correlation matrix).

#' @description The covariance matrix is given by Sigma = D R D, where the diagonal elements of D are given by sd, and R is constructed by forming a lower triangular matrix L with diagonals set to one and the remaining elements filled using rho column-wise, and then setting R = cov2cor(crossprod(L)).
#' @noRd
#' @noMd

.convert_to_covariance <- function(sd, rho) {
    num_ranef <- length(sd)

    L <- diag(num_ranef)
    L[lower.tri(L)] <- rho
    C <- cov2cor(tcrossprod(L))
    diag(C) <- 1

    out <- diag(x = sd) %*% C %*% diag(x = sd)
    return(out)
    }



#' @title Hidden internal function that creates the reparametrization matrix
#'
#' @description Basically, when the outputted matrix is right multiplied by the vector of regression coefficients (which are ordered as responses within covariates), it produces the reparametrized coefficients ready for homogeneity pursuit.
#' @note Function adapted from and credit goes to the authors of the \code{metafuse} package
#' @noRd
#' @noMd

.reparametrizationmatrix <- function(coefficients) {
    Smat <- .orderingmatrix(coefficients = coefficients)
    Rmat <- .adjustedifferencematrix(coefficients = coefficients, S = Smat)

    return(Rmat %*% Smat)
    }


#' @title Hidden internal function that creates an ordering matrix.
#'
#' @description Basically, when the outputted matrix is right multiplied by the vector of regression coefficients (which are ordered as responses within covariates), it then orders coefficients within each covariate.
#' @note Function adapted from and credit goes to the authors of the \code{metafuse} package
#' @noRd
#' @noMd

.orderingmatrix <- function(coefficients) {
    num_resp <- nrow(coefficients)
    num_X <- ncol(coefficients)

    S_list <- list()
    for(i in 1:num_X) {
        get_order <- order(coefficients[,i])
        Stemp <- matrix(0, nrow = num_resp, ncol = num_resp)
        for(j in 1:num_resp) {
            Stemp[j, get_order[j]] <- 1
        }

        S_list[[i]] <- Stemp
    }

    out <- bdiag(S_list)

    # Since out currently runs on a per-covariate basis, then reorder to run it on a per-response basis
    response_index <- rep(1:num_resp, num_X)
    out <- out[order(response_index),order(response_index)]
    return(out)
    }


#' @title Hidden internal function that creates an differences matrix.
#'
#' @description Basically, when the outputted matrix is right multiplied Smat, and then right multiplied by the vector of regression coefficients (which are ordered as responses within covariates), it produces the reparametrized coefficients
#' @note Function adapted from and credit goes to the authors of the \code{metafuse} package
#' @noRd
#' @noMd

.adjustedifferencematrix <- function(coefficients, S) {
    num_resp <- nrow(coefficients)
    num_X <- ncol(coefficients)
    new_coefficients <- matrix(S %*% as.vector(t(coefficients)), nrow = num_resp, byrow = TRUE)

    R_list <- c()
    for(i in 1:num_X) {
        Rtemp <- diag(x = c(0, rep(1, num_resp - 1)))
        Rtemp[1, which.min(abs(new_coefficients[,i]))] <- 1
        for(k in 2:num_resp)
            Rtemp[k, k-1] <- -1

        R_list[[i]] <- Rtemp
        }

    out <- bdiag(R_list)

    # Since out currently runs on a per-covariate basis, then reorder to run it on a per-response basis
    response_index <- rep(1:num_resp, num_X)
    out <- out[order(response_index),order(response_index)]
    return(out)
    }

