//-----------------------------
#include <math.h> 
#include <TMB.hpp> 
//#include <boost/accumulators/statistics/density.hpp>

using namespace density;
template<class Type>

Type objective_function<Type>::operator() () {  
     // Declare data
     DATA_VECTOR(y);
     DATA_MATRIX(X);
     DATA_MATRIX(Z);
     DATA_IVECTOR(clusid);
     DATA_IVECTOR(trial_size); // For binomial responses
     DATA_INTEGER(num_clus);
     int num_fixef = X.cols();
     int num_ranef = Z.cols();
     int N = y.size();


	 // Declare parameters
     PARAMETER_VECTOR(fixed_effects);
     PARAMETER_MATRIX(random_effects);
     PARAMETER_VECTOR(sd);
     PARAMETER_VECTOR(unconstrained_cor_params);

     // Negative log-likelihood
     Type nll = Type(0.0);
	 vector<Type> eta = X * fixed_effects;
     for(int i=0; i<N; i++) { for(int l1=0; l1<num_ranef; l1++) { eta(i) += Z(i,l1) * random_effects(clusid(i),l1);} }
     for(int i=0; i<N; i++) { nll -= dbinom(y(i), Type(trial_size(0)), invlogit(eta(i)), true); }

     for(int i=0; i<num_clus; i++) {
          nll += VECSCALE(UNSTRUCTURED_CORR(unconstrained_cor_params), sd)(random_effects.row(i));
          }


     matrix<Type> random_effects_covariance(num_ranef,num_ranef);
     random_effects_covariance = UNSTRUCTURED_CORR(unconstrained_cor_params).cov();

     ADREPORT(fixed_effects);
     ADREPORT(random_effects);
     ADREPORT(sd);
     ADREPORT(unconstrained_cor_params);
     
     return nll;
	 }
