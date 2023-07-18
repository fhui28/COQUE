//-----------------------------
#include <math.h> 
#include <TMB.hpp> 
//#include <boost/accumulators/statistics/density.hpp>

using namespace density;
template<class Type>

Type objective_function<Type>::operator() () {  
     // Declare data
     DATA_IVECTOR(y);
     DATA_MATRIX(X);
     DATA_MATRIX(Z);
     DATA_IVECTOR(clusid);
     DATA_IVECTOR(trial_size); // For binomial responses
     DATA_INTEGER(num_clus);
     DATA_IVECTOR(num_ordinal_levels);
     DATA_INTEGER(max_levels);
     int num_fixef = X.cols();
     int num_ranef = Z.cols();
     int N = y.size();


	 // Declare parameters
     PARAMETER_VECTOR(fixed_effects);
     PARAMETER_MATRIX(random_effects);
     PARAMETER_VECTOR(difford); //Number of cutoffs - 1 due to first cutoff set at zero
     vector<Type> ordinal_cutoffs(max_levels-1);
     ordinal_cutoffs(0) = 0;
     for(int l1=1; l1<(num_ordinal_levels(0)-1); l1++) {
         ordinal_cutoffs(l1) = ordinal_cutoffs(l1-1) + difford(l1-1); 
         }
     PARAMETER_VECTOR(sd);
     PARAMETER_VECTOR(unconstrained_cor_params);

     // Negative log-likelihood
     Type nll = Type(0.0);
	 vector<Type> eta = X * fixed_effects;
     for(int i=0; i<N; i++) { for(int l1=0; l1<num_ranef; l1++) { eta(i) += Z(i,l1) * random_effects(clusid(i),l1);} }
     matrix<Type> prob_ord(N, num_ordinal_levels(0));
     for(int i=0; i<N; i++) {
          prob_ord(i,0) = invlogit(ordinal_cutoffs(0)-eta(i));
          for(int l=1; l<(num_ordinal_levels(0)-1); l++) {
               prob_ord(i,l) = invlogit(ordinal_cutoffs(l) - eta(i)) - invlogit(ordinal_cutoffs(l-1) - eta(i));
               }
          prob_ord(i,num_ordinal_levels(0)-1) = 1 - invlogit(ordinal_cutoffs(num_ordinal_levels(0)-2) - eta(i));
          nll -= log(prob_ord(i,y(i)-1));
          }

     for(int i=0; i<num_clus; i++) {
          nll += VECSCALE(UNSTRUCTURED_CORR(unconstrained_cor_params), sd)(random_effects.row(i));
          }


     matrix<Type> random_effects_covariance(num_ranef,num_ranef);
     random_effects_covariance = UNSTRUCTURED_CORR(unconstrained_cor_params).cov();

     ADREPORT(fixed_effects);
     ADREPORT(difford);
     ADREPORT(random_effects);
     ADREPORT(sd);
     ADREPORT(unconstrained_cor_params);
     
     return nll;
	 }
