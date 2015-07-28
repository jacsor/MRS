#ifndef HELPERS_H
#define HELPERS_H

#define MAXVAR 15  // down no more than 14 levels
#include "RcppArmadillo.h"


using namespace Rcpp;
using namespace arma;


union INDEX_TYPE_t 
{
  //var[MAXVAR] gives the bits that represent the left and right children
  unsigned short var[MAXVAR+1];   
  unsigned long long index;       // probably not used
};
typedef INDEX_TYPE_t INDEX_TYPE;


struct side_type
{
  unsigned short var;
	unsigned short extremes[2];
};


struct cube_type
{
  std::vector<side_type> sides;
  ushort level;
  double alt_prob;
  vec effect_size;         
  int direction;
  ushort node_idx;  
  arma::Col< unsigned > data_points;
};

typedef std::vector<cube_type> result_cubes_type;



double log_exp_x_plus_exp_y(double x, double y);

unsigned long long pow2(int k); 

uint convert_to_inverse_base_2(double x, int k);

unsigned long long int Choose(int n, int k);

INDEX_TYPE init_index(int level); 

unsigned long long int get_node_index(INDEX_TYPE& I,
                                      int level, 
                                      int dim);
                                      
std::pair<bool, INDEX_TYPE> make_parent_index( INDEX_TYPE& I, 
                                                unsigned short part_dim, 
                                                int level, 
                                                ushort which );

INDEX_TYPE make_child_index(  INDEX_TYPE& I, 
                              unsigned short part_dim, 
                              int level, 
                              ushort which); 
                              
int sum_elem(int * my_array , int num_elem);

INDEX_TYPE get_next_node( INDEX_TYPE& I, 
                          int p, 
                          int level); 

#endif
