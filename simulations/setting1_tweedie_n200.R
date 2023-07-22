#' ---
#' title: Simulation setting for comparing COQUE to some current methods for fixed effects selection in multivariate independent cluster GLMMs
#' abstract: Tweedie responses. Note the code was also set up to run as an array job on a HPC.
#' author: Francis KC Hui
#' ---

rm(list = ls())
here::i_am("simulations/setting1_tweedie_n200.R")
library(here)

#' Code for COQUE
source(here("code", "coque.R"))
source(here("code","estimates.R"))
source(here("code","stackedGLMMs.R"))
source(here("code","gendat.R"))
source(here("code","utils.R"))

#' Code for alternative methods
#source(here("code","alternative_methods","stackedglmmLasso.R"))
source(here("code", "alternative_methods", "backwardelim_BIC.R"))
#source(here("code", "alternative_methods", "stackedglmmPen.R"))
#source(here("code", "alternative_methods", "stackedrpql.R"))


##------------------
#' # Simulate data -- inspired by [http://dx.doi.org/10.1016/j.csda.2017.09.004]
##------------------
set.seed(072023)
num_clus <- 200
num_resp <- 6
response_type <- "tweedie"
dat <- data.frame(id = rep(1:num_clus, sample(10:15, size = num_clus, replace = TRUE)))

p <- 15
H <- abs(outer(1:p, 1:p, "-"))
X <- rmvnorm(nrow(dat), sigma = 0.5^H)
colnames(X) <- paste0("x", 1:ncol(X))
dat <- cbind(dat, X)
rm(X)
str(dat)

true_fixed_effects <- cbind(
   runif(num_resp, -1, 1),
   rnorm(num_resp,0.5,sd=0.5),
   rnorm(num_resp,-0.5,sd=0.5),
   replicate(4, sample(c(-0.25,0,0.25),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
   replicate(4, sample(c(-0.125,0),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 2 unique values, including zero
   replicate(5, numeric(num_resp)))
true_fixed_effects


#' Adapted from HPGFE manuscript
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

rm(H, p)


dosims <- function(NAI) {
    set.seed(NAI)
    
    simy <- generatedata_mglmm(formula_X = ~ . - id,
                               formula_Z = ~ x1 + x2 + x3 + x4,
                               data = dat,
                               fixed_effects = true_fixed_effects,
                               Sigma = true_Sigma,
                               response_type = response_type,
                               phi = true_phi,
                               power = true_power,
                               trial_size = 1)


    ##------------------
    #' # Apply COQUE
    ##------------------
    tic <- proc.time()
    sglmms <- stackedGLMMs(y = simy$y[,-ncol(simy$y)],
                           response_type = response_type,
                           formula_X = ~ . - id,
                           formula_Z = ~ x1 + x2 + x3 + x4,
                           data = dat,
                           trial_size = 1,
                           num_cores = num_resp,
                           TMB_directory = here("code","cppfiles")
                           )

    penalized_fit <- try(coque(object = sglmms, num_cores = num_resp), silent = TRUE)
    toc <- proc.time()
    penalized_fit$time_taken_full <- toc - tic

  
    ##------------------
    #' # Stacked penalized GLMMs using the [glmmLasso::glmmLasso()] package.
    #' Note it is possible to apply this, but in simulations we found this to be extremely computationally burdensome and so not consider it.
    ##------------------
    # tic <- proc.time()
    # glmmlasso_fit <- try(stacked_glmmLasso(fix = resp ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15,
    #                                        rnd = list(id = ~ 1 + x1 + x2 + x3 + x4),
    #                                        y = simy$y[,-ncol(simy$y)],
    #                                        data = dat,
    #                                        num_cores = num_resp,
    #                                        response_type = response_type,
    #                                        tweedie_var_power = true_power, 
    #                                        lambda = 10^seq(-3,5, length = 100),
    #                                        steps = 500), 
    #                      silent = TRUE)
    # toc <- proc.time()
    # if(!inherits(glmmlasso_fit, "try-error"))
    #     glmmlasso_fit$time_taken_full <- toc - tic

    ##------------------
    #' # Stacked GLMMs with fixed effects backward elimination approach using BIC. 
    ##------------------
    tic <- proc.time()
    glmm_backwardelim <- try(stacked_backwardelim(y = simy$y[,-ncol(simy$y)],
                                                  data = dat,
                                                  response_type = response_type,
                                                  num_fixef = ncol(true_fixed_effects)-1,
                                                  formula_Z = "(x1 + x2 + x3 + x4 | id)",
                                                  num_cores = num_resp),
                             silent = TRUE)
    toc <- proc.time()
    if(!inherits(glmm_backwardelim, "try-error"))
        glmm_backwardelim$time_taken_full <- toc - tic

    ##------------------
    #' ## Save results
    ##------------------
    save(penalized_fit, 
         #glmmlasso_fit, 
         glmm_backwardelim, 
         simy,
         file = paste0("setting1_tweedie_n", num_clus, "_dataset", NAI, ".RData"))
  }


dosims(NAI = NAI) # NAI here would be the index of the array job in HPC



##------------------------                
sessionInfo()
##------------------------                
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.9 (Final)
# 
# Matrix products: default
# BLAS:   /usr/local/R/4.2.2/lib64/R/lib/libRblas.so
# LAPACK: /usr/local/R/4.2.2/lib64/R/lib/libRlapack.so
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#  [1] rpql_0.8           ncvreg_3.14.1      glmmPen_1.5.3.4    Rcpp_1.0.10       
#  [5] bigmemory_4.6.1    lme4_1.1-33        ordinal_2022.11-16 glmmTMB_1.1.7     
#  [9] statmod_1.5.0      glmmLasso_1.6.2    tweedie_2.3.5      doParallel_1.0.17 
# [13] iterators_1.0.14   foreach_1.5.2      TMB_1.9.4          quadprog_1.5-8    
# [17] mvtnorm_1.2-2      Matrix_1.5-4.1     here_1.0.1         colorspace_2.1-0  
# 
# loaded via a namespace (and not attached):
#  [1] splines_4.2.2       StanHeaders_2.26.27 RcppParallel_5.1.7 
#  [4] ucminf_1.2.0        stats4_4.2.2        numDeriv_2016.8-1.1
#  [7] pillar_1.9.0        lattice_0.21-8      glue_1.6.2         
# [10] uuid_1.1-0          minqa_1.2.5         sandwich_3.0-2     
# [13] plyr_1.8.8          pkgconfig_2.0.3     rstan_2.21.8       
# [16] xtable_1.8-4        scales_1.2.1        processx_3.8.1     
# [19] emmeans_1.8.6       tibble_3.2.1        generics_0.1.3     
# [22] ggplot2_3.4.2       TH.data_1.1-2       cli_3.6.1          
# [25] survival_3.5-5      magrittr_2.0.3      crayon_1.5.2       
# [28] estimability_1.4.1  ps_1.7.5            fansi_1.0.4        
# [31] nlme_3.1-162        MASS_7.3-60         pkgbuild_1.4.1     
# [34] loo_2.6.0           prettyunits_1.1.1   tools_4.2.2        
# [37] coxme_2.2-18.1      matrixStats_1.0.0   lifecycle_1.0.3    
# [40] multcomp_1.4-24     stringr_1.5.0       munsell_0.5.0      
# [43] gamlss.dist_6.0-5   callr_3.7.3         compiler_4.2.2     
# [46] rlang_1.1.1         grid_4.2.2          nloptr_2.0.3       
# [49] bigmemory.sri_0.1.6 boot_1.3-28.1       gtable_0.3.3       
# [52] codetools_0.2-19    inline_0.3.19       reshape2_1.4.4     
# [55] R6_2.5.1            gridExtra_2.3       rstantools_2.3.1   
# [58] zoo_1.8-12          dplyr_1.1.2         bdsmatrix_1.3-6    
# [61] utf8_1.2.3          rprojroot_2.0.3     stringi_1.7.12     
# [64] vctrs_0.6.3         tidyselect_1.2.0    coda_0.19-4        

