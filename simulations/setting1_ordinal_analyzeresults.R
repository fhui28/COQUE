#' ---
#' title: Simulation setting for comparing COQUE to some current methods for fixed effects selection in multivariate independent cluster GLMMs -- analyzing results
#' abstract: Ordinal responses. Note the code was also set up to run as an array job on a HPC.
#' author: Francis KC Hui
#' ---

rm(list = ls())
library(tidyverse)
library(colorspace)
library(patchwork)
library(Matrix)
library(quadprog)
library(mvtnorm)
library(ROCR)
here::i_am("simulations/setting1_ordinal_analyzeresults.R")
library(here)

#' Code for BARF and composite quadratic estimator
source(here("code", "coque.R"))
source(here("code","estimates.R"))
source(here("code","stackedGLMMs.R"))
source(here("code","gendat.R"))
source(here("code","utils.R"))


##------------------
#' # Simulate data -- inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
##------------------
set.seed(072023)
num_resp <- 6
response_type <- "ordinal"
p <- 15
H <- abs(outer(1:p, 1:p, "-"))
num_clus <- 100

dat <- data.frame(id = rep(1:num_clus, sample(5:10, size = num_clus, replace = TRUE)))
X <- rmvnorm(nrow(dat), sigma = 0.5^H)
colnames(X) <- paste0("x", 1:ncol(X))
dat <- cbind(dat, X)
rm(X)
str(dat)

true_fixed_effects_n100 <- cbind(
    runif(num_resp, -1, 1),
    rnorm(num_resp,0.5,sd=0.5),
    rnorm(num_resp,-0.5,sd=0.5),
    replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
    replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 2 unique values, including zero
    replicate(5, numeric(num_resp)))

true_distinct_values_nointercept_n100 <- apply(true_fixed_effects_n100[,-1], 2, function(x) length(unique(x))) 



set.seed(072023)
num_resp <- 6
response_type <- "ordinal"
p <- 15
H <- abs(outer(1:p, 1:p, "-"))
num_clus <- 200

dat <- data.frame(id = rep(1:num_clus, sample(10:15, size = num_clus, replace = TRUE)))
X <- rmvnorm(nrow(dat), sigma = 0.5^H)
colnames(X) <- paste0("x", 1:ncol(X))
dat <- cbind(dat, X)
rm(X)
str(dat)

true_fixed_effects_n200 <- cbind(
    runif(num_resp, -1, 1),
    rnorm(num_resp,0.5,sd=0.5),
    rnorm(num_resp,-0.5,sd=0.5),
    replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
    replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 2 unique values, including zero
    replicate(5, numeric(num_resp)))

true_distinct_values_nointercept_n200 <- apply(true_fixed_effects_n200[,-1], 2, function(x) length(unique(x))) 


true_fixed_effects <- abind::abind(true_fixed_effects_n100, true_fixed_effects_n200, along = 3)
true_distinct_values_nointercept <- abind::abind(true_distinct_values_nointercept_n100, true_distinct_values_nointercept_n200, along = 2)



true_G <- rbind(c(1,0.8,0.4,-0.6,-0.3,-0.15),
                c(0.8,1,0.6,-0.4,-0.2,-0.1),
                c(0.4,0.6,1,-0.2,-0.1,-0.05),
                c(-0.6,-0.4,-0.2,1,0.6,0.3),
                c(-0.3,-0.2,-0.1,0.6,1,0.3),
                c(-0.15,-0.1,-0.05,0.3,0.3,1))
true_ranef_cov <- rbind(c(2,1,0,0,0),
                        c(1,2,0,0,0),
                        c(0,0,1,0,0),
                        c(0,0,0,0.5,0),
                        c(0,0,0,0,0.25))/2
true_Sigma <- kronecker(true_G, true_ranef_cov)
true_phi <- runif(num_resp,0,2)
true_power <- runif(num_resp, 1.2, 1.8)

rm(H, p, dat, true_distinct_values_nointercept_n100, true_distinct_values_nointercept_n200, true_fixed_effects_n100, true_fixed_effects_n200)
rm(true_G, true_Sigma, true_ranef_cov, num_clus, true_phi, true_power)



##------------------
#' # Process results
##------------------
N_seq <- c(100,200)
num_datasets <- 500
methods <- c("glmmLasso", "backward_elimination", "COQUE")
all_time_taken <- all_RMSE <- all_MSE_uniquevalues <- all_sensitivity <- 
    all_specificity <- all_accuracy <- all_F1 <- array(NA, 
                                                       dim = c(2, num_datasets, length(methods)), 
                                                       dimnames = list(N = N_seq, datasets = 1:num_datasets, methods = methods))


for(k1 in 1:2) { for(k0 in 1:num_datasets) {
    filename <- paste0("setting1_ordinal_n", N_seq[k1], "_dataset", k0, ".RData")
    message("Onto dataset ", k0)
    if(file.exists(filename))
        load(filename)
    if(!file.exists(filename))
        next;
    
    
    #' Pass over any erroneous fits (testing suggesting this was very little e.g., 2-4 in 500 datasets)
    if(inherits(try(sapply(glmmlasso_fit[1:num_resp], function(x) x$final_fixed_effects), silent = TRUE), "try-error"))
        next;
    if(inherits(try(sapply(glmm_backwardelim[1:num_resp], function(x) x$coefficients), silent = TRUE), "try-error"))
        next;
    if(inherits(penalized_fit, "try-error"))
        next;
    
    
    #' Computation time
    all_time_taken[k1,k0,] <- c(
        glmmlasso_fit$time_taken_full[3],
        glmm_backwardelim$time_taken_full[3],
        penalized_fit$time_taken_full[3]
        )
    
    
    #' Metrics for each of the estimators
    cw_coefficients <- t(sapply(glmmlasso_fit[1:num_resp], function(x) x$final_fixed_effects[-(2:4)])) # Remove all but 1 of the cutoff parameters 
    makepred_obj <- ROCR::prediction(predictions = as.vector(1*(cw_coefficients != 0)), labels = as.vector(1*(true_fixed_effects[,,k1] != 0)))
    cw_uniquevalues_nointercept <- apply(cw_coefficients[,-1], 2, function(x) length(unique(x))) 
    all_sensitivity[k1,k0,1] <- ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2]
    all_specificity[k1,k0,1] <- ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2]
    all_accuracy[k1,k0,1] <- ROCR::performance(makepred_obj, measure = "acc")@y.values[[1]][2]
    all_F1[k1,k0,1] <- ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2]
    all_RMSE[k1,k0,1] <- sqrt(mean(cw_coefficients[,-1] - true_fixed_effects[,,k1][,-1])^2)
    all_MSE_uniquevalues[k1,k0,1] <- mean((cw_uniquevalues_nointercept - true_distinct_values_nointercept[,k1])^2)
    
    
    cw_coefficients <- t(sapply(glmm_backwardelim[1:num_resp], function(x) x$coefficients))
    makepred_obj <- ROCR::prediction(predictions = as.vector(1*(cw_coefficients != 0)), labels = as.vector(1*(true_fixed_effects[,,k1] != 0)))
    cw_uniquevalues_nointercept <- apply(cw_coefficients[,-1], 2, function(x) length(unique(x))) 
    all_sensitivity[k1,k0,2] <- ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2]
    all_specificity[k1,k0,2] <- ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2]
    all_accuracy[k1,k0,2] <- ROCR::performance(makepred_obj, measure = "acc")@y.values[[1]][2]
    all_F1[k1,k0,2] <- ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2]
    all_RMSE[k1,k0,2] <- sqrt(mean(cw_coefficients[,-1] - true_fixed_effects[,,k1][,-1])^2)
    all_MSE_uniquevalues[k1,k0,2] <- mean((cw_uniquevalues_nointercept - true_distinct_values_nointercept[,k1])^2)

    
    cw_coefficients <- estimates(penalized_fit, lambda = "BIC", hybrid = TRUE)$hybrid
    makepred_obj <- ROCR::prediction(predictions = as.vector(1*(cw_coefficients != 0)), labels = as.vector(1*(true_fixed_effects[,,k1] != 0)))
    cw_uniquevalues_nointercept <- apply(cw_coefficients[,-1], 2, function(x) length(unique(x))) 
    all_sensitivity[k1,k0,3] <- ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2]
    all_specificity[k1,k0,3] <- ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2]
    all_accuracy[k1,k0,3] <- ROCR::performance(makepred_obj, measure = "acc")@y.values[[1]][2]
    all_F1[k1,k0,3] <- ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2]
    all_RMSE[k1,k0,3] <- sqrt(mean(cw_coefficients[,-1] - true_fixed_effects[,,k1][,-1])^2)
    all_MSE_uniquevalues[k1,k0,3] <- mean((cw_uniquevalues_nointercept - true_distinct_values_nointercept[,k1])^2)

    rm(list = ls(pattern = "glmm"))
    rm(penalized_fit, simy, cw_uniquevalues_nointercept, cw_coefficients)
    } }



##-----------------------
#' # Examine results
##-----------------------
apply(all_RMSE[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_RMSE[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), sd, na.rm = TRUE) %>% round(3)

apply(all_F1[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_F1[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), sd, na.rm = TRUE) %>% round(3)

apply(all_MSE_uniquevalues[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_MSE_uniquevalues[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), sd, na.rm = TRUE) %>% round(3)

apply(all_time_taken[,,c("COQUE", "glmmLasso", "backward_elimination")]/60, c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_time_taken[,,c("COQUE", "glmmLasso", "backward_elimination")]/60, c(1,3), sd, na.rm = TRUE) %>% round(3)


apply(all_sensitivity[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_sensitivity[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), sd, na.rm = TRUE) %>% round(3)

apply(all_specificity[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_specificity[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), sd, na.rm = TRUE) %>% round(3)

apply(all_accuracy[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_accuracy[,,c("COQUE", "glmmLasso", "backward_elimination")], c(1,3), sd, na.rm = TRUE) %>% round(3)

