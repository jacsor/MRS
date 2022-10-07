#include "RcppArmadillo.h"
#include "helpers.h"
#include "recursion.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
Rcpp::List fitMRScpp( arma::mat X, 
                      arma::vec G, 
                      int n_groups, 
                      arma::vec init_state,
                      arma::mat Omega,  
                      int K = 5, 
                      double alpha = 0.5, 
                      double beta = 1.0, 
                      double gamma = 0.3,
                      double delta = 0.3,
                      double eta = 0.3,
                      bool return_global_null = true,
                      bool return_tree = true,
                      int n_post_samples = 0,
                      int baseline = 0,
                      int min_n_node = 0
                 )
{

  /* **************************************** */
  /* center the data and transform in binary  */
  /* **************************************** */
  
  // number of dimensions
  int p = X.n_cols;
  // total number of observations
  int n_tot = X.n_rows;
  // total number of observations in each subgroup (used in dANOVA)
  Col<int> n_subgroups(n_groups); n_subgroups.fill(1);     
  // sub-gorups label for each observation (used in dANOVA)
  vec H(n_tot); H.fill(1);
  // prior on parameter nu (used in dANOVA)
  arma::vec nu_vec({1}); 
  vec a = 1.0 / (Omega.col(1) - Omega.col(0));
  vec b = - Omega.col(0) % a;
  // Initialize matrix of observations normalized to the p-dimensional hypercube
  Mat<unsigned int> X_binary(n_tot,p);    
  // Map each observation to the p-dimensional hypercube
  // Transform each observation in p-dimensional binary vector (up to resolution K+1)
  for(int i = 0; i < n_tot; i++ )
  {
    for(int d = 0; d < p; d++ )
      X_binary(i,d) = convert_to_inverse_base_2(a(d)*X(i,d)+b(d), K+1 );
  }
    
  /* ***************** */
  /* Compute posterior */
  /* ***************** */
  class_tree my_tree( X_binary, 
                      G, 
                      H,
                      init_state, 
                      n_groups, 
                      n_subgroups,
                      K, 
                      nu_vec,
                      alpha, 
                      beta, 
                      gamma, 
                      delta,
                      eta, 
                      return_global_null, 
                      return_tree,
                      n_post_samples,
                      baseline,
                      min_n_node,
                      0,  // Method is not used for MRS, use 0 as place holder
                      1   // n_grid_theta is not used for MRS, use 1 as place holder
                     );
  my_tree.update();
  
  double loglike = my_tree.get_marginal_loglikelihood();

  double post_glob_null = NA_REAL;  
  double prior_glob_null = my_tree.get_prior_global_null();
  if(return_global_null == true)
    post_glob_null = my_tree.get_posterior_global_null();
    
  Rcpp::List data = Rcpp::List::create(
    Rcpp::Named( "X") = X,
    Rcpp::Named( "G") = G,
    Rcpp::Named( "Groups") = n_groups,
    Rcpp::Named( "Omega" ) = Omega
  );
    
  Rcpp::List other_info = Rcpp::List::create(
    Rcpp::Named("K") = K,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("delta") = delta
  );
  
  Rcpp::List post_sample_trees(n_post_samples);
    
  if (n_post_samples > 0) {
    
    my_tree.sample_tree();
    vector<unsigned short> levels_sample;
    vector<double> alt_probs_sample;
    vector<vec> effect_sizes_sample; 
    vector<int> directions_sample;
    vector< vector<double> > sides_sample;
    vector< Col< unsigned > > data_points_sample;
    vector<unsigned short> node_idx_sample;
    
    for (int sample_id = 0; sample_id < n_post_samples; sample_id++) {
    
      levels_sample = my_tree.get_level_sample_nodes(sample_id);
      alt_probs_sample = my_tree.get_alt_prob_sample_nodes(sample_id);
      effect_sizes_sample = my_tree.get_effect_size_sample_nodes(sample_id); 
      directions_sample = my_tree.get_direction_sample_nodes(sample_id);
      sides_sample = my_tree.get_sides_sample_nodes(a, b,sample_id);
      data_points_sample =  my_tree.get_data_points_sample_nodes(sample_id);
      node_idx_sample = my_tree.get_idx_sample_nodes(sample_id);
      
      post_sample_trees[sample_id] = Rcpp::List::create( 
        Rcpp::Named( "Levels") = levels_sample,
        Rcpp::Named( "AltProbs") = alt_probs_sample,
        Rcpp::Named( "EffectSizes") = effect_sizes_sample,
        Rcpp::Named( "Directions") = directions_sample,
        Rcpp::Named( "Regions") = sides_sample,
        Rcpp::Named( "DataPoints") = data_points_sample,
        Rcpp::Named( "Ids") = node_idx_sample 
      );
    }
  }
  
  else {
    post_sample_trees = Rcpp::List::create();
  }
  
  if(return_tree == true)
  {
    my_tree.compute_varphi_post();  
    my_tree.representative_tree();     
    vector<unsigned short> levels = my_tree.get_level_nodes();
    vector<double> alt_probs = my_tree.get_alt_prob_nodes();
    vector<vec> effect_sizes = my_tree.get_effect_size_nodes(); 
    vector<int> directions = my_tree.get_direction_nodes();
    vector< vector<double> > sides = my_tree.get_sides_nodes(a, b);   
    vector< Col< unsigned > > data_points =  my_tree.get_data_points_nodes();
    vector<unsigned short> node_idx = my_tree.get_idx_nodes();
    
    my_tree.clear();
        
    Rcpp::List rep_tree = Rcpp::List::create( 
      Rcpp::Named( "Levels") = levels,
      Rcpp::Named( "AltProbs") = alt_probs,
      Rcpp::Named( "EffectSizes") = effect_sizes,
      Rcpp::Named( "Directions") = directions,
      Rcpp::Named( "Regions") = sides,
      Rcpp::Named( "DataPoints") = data_points,
      Rcpp::Named( "Ids") = node_idx 
    );
    
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "PostSamples" ) = post_sample_trees,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   
  }
  else
  {
    my_tree.clear();
    Rcpp::List rep_tree = Rcpp::List::create();
        
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "PostSamples" ) = post_sample_trees,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   

  }
  


    
}





// [[Rcpp::export]]
Rcpp::List fitMRSNESTEDcpp( arma::mat X, 
                            arma::vec G, 
                            arma::vec H, 
                            int n_groups, 
                            arma::Col<int> n_subgroups,
                            arma::vec init_state,
                            arma::mat Omega, 
                            arma::vec nu_vec,
                            int K = 5,                             
                            double alpha = 0.5, 
                            double beta = 1.0, 
                            double gamma = 0.07, 
                            double delta = 0.4,
                            double eta = 0,
                            bool return_global_null = true,
                            bool return_tree = true,
                            int n_post_samples = 0,
                            int baseline = 0,
                            int method = 0,  // method of integration over theta (0 is Newton-Raphson, 1 is Riemann quadrature)
                            int n_grid_theta = 4 // the number of theta values used in the Riemann quadrature either when Riemann quadrature is chosen or when NR has failed.
                          )
{

  /* **************************************** */
  /* center the data and transform in binary  */
  /* **************************************** */
  
  // number of dimensions
  int p = X.n_cols;
  // total number of observations
  int n_tot = X.n_rows;
  vec a = 1.0 / (Omega.col(1) - Omega.col(0));
  vec b = - Omega.col(0) % a;
  // Initialize matrix of observations normalized to the p-dimensional hypercube
  Mat<unsigned int> X_binary(n_tot,p);    
  // Map each observation to the p-dimensional hypercube
  // Transform each observation in p-dimensional binary vector (up to resolution K+1)
  for(int i = 0; i < n_tot; i++ )
  {
    for(int d = 0; d < p; d++ )
      X_binary(i,d) = convert_to_inverse_base_2(a(d)*X(i,d)+b(d), K+1 );
  }
    
  /* ***************** */
  /* Compute posterior */
  /* ***************** */
  class_tree my_tree( X_binary, 
                      G, 
                      H,
                      init_state, 
                      n_groups, 
                      n_subgroups,
                      K, 
                      nu_vec,
                      alpha, 
                      beta, 
                      gamma, 
                      delta,
                      eta, 
                      return_global_null, 
                      return_tree,
                      n_post_samples,
                      baseline,
                      0, // n_min_node is not used for ANDOVA so use 0 as place holder in the constructor
                      method, 
                      n_grid_theta); 
  my_tree.update();
  
  double loglike = my_tree.get_marginal_loglikelihood();

  double post_glob_null = NA_REAL;  
  double prior_glob_null = my_tree.get_prior_global_null();
  if(return_global_null == true)
    post_glob_null = my_tree.get_posterior_global_null();
    
  Rcpp::List data = Rcpp::List::create(
    Rcpp::Named( "X") = X,
    Rcpp::Named( "G") = G,
    Rcpp::Named( "H") = H,
    Rcpp::Named( "Groups") = n_groups,
    Rcpp::Named( "Replicates") = n_subgroups,
    Rcpp::Named( "Omega" ) = Omega
  );
    
  Rcpp::List other_info = Rcpp::List::create(
    Rcpp::Named("K") = K,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("delta") = delta,
    Rcpp::Named("nu_vec") = nu_vec
  );
  
  Rcpp::List post_sample_trees(n_post_samples);
  
  if (n_post_samples > 0) {
    
    my_tree.sample_tree();
    vector<unsigned short> levels_sample;
    vector<double> alt_probs_sample;
    vector<vec> effect_sizes_sample; 
    vector<int> directions_sample;
    vector< vector<double> > sides_sample;
    vector< Col< unsigned > > data_points_sample;
    vector<unsigned short> node_idx_sample;
    
    for (int sample_id = 0; sample_id < n_post_samples; sample_id++) {
      
      levels_sample = my_tree.get_level_sample_nodes(sample_id);
      alt_probs_sample = my_tree.get_alt_prob_sample_nodes(sample_id);
      effect_sizes_sample = my_tree.get_effect_size_sample_nodes(sample_id); 
      directions_sample = my_tree.get_direction_sample_nodes(sample_id);
      sides_sample = my_tree.get_sides_sample_nodes(a, b,sample_id);
      data_points_sample =  my_tree.get_data_points_sample_nodes(sample_id);
      node_idx_sample = my_tree.get_idx_sample_nodes(sample_id);
      
      post_sample_trees[sample_id] = Rcpp::List::create( 
        Rcpp::Named( "Levels") = levels_sample,
        Rcpp::Named( "AltProbs") = alt_probs_sample,
        Rcpp::Named( "EffectSizes") = effect_sizes_sample,
        Rcpp::Named( "Directions") = directions_sample,
        Rcpp::Named( "Regions") = sides_sample,
        Rcpp::Named( "DataPoints") = data_points_sample,
        Rcpp::Named( "Ids") = node_idx_sample 
      );
    }
  }
  
  else {
    post_sample_trees = Rcpp::List::create();
  }
  
  
  if(return_tree == true)
  {
    my_tree.compute_varphi_post();  
    my_tree.representative_tree();     
    vector<unsigned short> levels = my_tree.get_level_nodes();
    vector<double> alt_probs = my_tree.get_alt_prob_nodes();
    vector<vec> effect_sizes = my_tree.get_effect_size_nodes(); 
    vector<int> directions = my_tree.get_direction_nodes();
    vector< vector<double> > sides = my_tree.get_sides_nodes(a, b);   
    vector< Col< unsigned > > data_points =  my_tree.get_data_points_nodes();
    vector<unsigned short> node_idx = my_tree.get_idx_nodes();
    
    my_tree.clear();
        
    Rcpp::List rep_tree = Rcpp::List::create( 
      Rcpp::Named( "Levels") = levels,
      Rcpp::Named( "AltProbs") = alt_probs,
      Rcpp::Named( "EffectSizes") = effect_sizes,
      Rcpp::Named( "Directions") = directions,
      Rcpp::Named( "Regions") = sides,
      Rcpp::Named( "DataPoints") = data_points,
      Rcpp::Named( "Ids") = node_idx 
    );
    
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "PostSamples" ) = post_sample_trees,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   
  }
  else
  {
    my_tree.clear();
    Rcpp::List rep_tree = Rcpp::List::create();
        
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "PostSamples" ) = post_sample_trees,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   

  }

    
}

