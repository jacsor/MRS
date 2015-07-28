
#include <RcppArmadillo.h>
#include "helpers.h" 



using namespace Rcpp;
using namespace arma;
using namespace std;


const double log2pi = std::log(2.0 * M_PI);

double log_exp_x_plus_exp_y(double x, double y) 
{
  double result;
  if( ( exp(x) == 0 ) && ( exp(y) == 0 ) )
    result = x;
  else if ( x - y >= 100 ) result = x;
  else if ( x - y <= -100 ) result = y;
  else {
    if (x > y) {
      result = y + log( 1 + exp(x-y) );
    }
    else result = x + log( 1 + exp(y-x) );
  }
  return result;
}

unsigned long long pow2(int k) 
{
  unsigned long long res = (unsigned long long) 1 << k;
  return res;
}

uint convert_to_inverse_base_2(double x, int k)
{
  uint x_base_2 = (uint) floor (x * pow2(k) );
  uint x_base_inverse_2 = 0;

  for (int i = 0; i < k; i++) 
  {
      x_base_inverse_2 |= ((x_base_2 >> (k - 1 -i)) & 1) << i;
  }

  return x_base_inverse_2;
}

unsigned long long int Choose(int n, int k) 
{
    unsigned long long int c = 1;
    unsigned long long int d = 1;
    for (int i = 0; i <k; i++) 
    {
        c *= (n-i); d *= (i+1);
    }
    return c  / d;
}

INDEX_TYPE init_index(int level) 
{   
    INDEX_TYPE init;
    for (int i=0; i < level; i++) {init.var[i] = i+1;}    
    for (int i=level; i <= MAXVAR; i++) {init.var[i] = 0;}    
    return init;
}

unsigned long long int get_node_index(INDEX_TYPE& I, int level, int dim)
{
  unsigned long long  r = 0;
  unsigned long long numerator = 1;
  unsigned long long denominator = 1;
  
  for (int i = 0; i < level; i++) 
  {
    numerator = 1;
    denominator *= (i+1);     
    for (int j = 1; j <= i+1; j++) 
      numerator *= I.var[i]-j;
    
    r += numerator / denominator;   
  }
  return dim*(r*pow2(level) + (unsigned long long) I.var[MAXVAR]);
    
}


std::pair<bool, INDEX_TYPE> make_parent_index( INDEX_TYPE& I, 
                                          unsigned short part_dim, 
                                          int level, 
                                          ushort which )
{
    INDEX_TYPE parent_index = I;
    INDEX_TYPE child_index = I;
    unsigned short data = part_dim+1; 
    int i = 0;
    int j = 0;
    int x_curr;
    int parent_index_var_prev;
    bool first_time = true;
    bool parent_exists = true;
    
    if (level == 0) 
      parent_exists = false;
    else
    {
      i = 0;
      x_curr = parent_index.var[i]; // current dimension
      parent_index_var_prev = x_curr - 1;
    }
    
    if(x_curr != data)
    {
      while (parent_index.var[i] > 0 && data > x_curr ) 
      {
        x_curr += parent_index.var[i] - parent_index_var_prev - 1;
        parent_index_var_prev = parent_index.var[i];
        i++;
      }
      i = i - 1;
    }

    if ( data != x_curr ) 
      parent_exists = false;
    else
    {
      while( parent_index.var[i] > 0 )
      {
        i++;
        x_curr += parent_index.var[i] - parent_index_var_prev - 1;
        parent_index_var_prev = parent_index.var[i];
        if( parent_index_var_prev == 0 )
          parent_index.var[i - 1] = 0;
        else if( x_curr > data )  
          parent_index.var[i - 1] = parent_index.var[i] - 1;
          
        if( first_time && parent_index.var[i-1]!= child_index.var[i-1])
        {
          j = i - 1;
          first_time = false;
        }
      }
      
    
      if( ( ( parent_index.var[MAXVAR] >> j) & (ushort)1 ) != (ushort)which )
        parent_exists = false;
      else
        parent_index.var[MAXVAR] = 
          ( I.var[MAXVAR] & ~((~((ushort)0)) << j ) ) | ((I.var[MAXVAR] >> 1) & ((~((ushort)0)) << j ));
          
    }

    return std::make_pair(parent_exists, parent_index);

}


INDEX_TYPE make_child_index(  INDEX_TYPE& I, 
                              unsigned short part_dim, 
                              int level, 
                              ushort which) {
    INDEX_TYPE child_index = I;
    unsigned short data = part_dim+1; 
    int i;
    int j;
    int x_curr;
    int child_index_var_prev;

    if (level == 0) 
    {
      x_curr=1;
      i=0;
      child_index_var_prev = 0;
    }
    else 
    {
       x_curr = child_index.var[0]; // current dimension
       child_index_var_prev = child_index.var[0];
       i = 1;
    }

    while (i<MAXVAR) 
    {
        while (child_index.var[i] >0 && data >= x_curr ) 
        {
          x_curr += child_index.var[i] - child_index_var_prev - 1;
	        child_index_var_prev = child_index.var[i];
	        i++;
        }

        if (child_index.var[i] == 0 && data >= x_curr) 
        {
	        child_index.var[i] = data - x_curr + 1 + child_index_var_prev;
	        j=i;
	        i=MAXVAR;
        } 
        else 
        { // this corresponds to the first i such that data < x_curr
	        for (int h = level; h >= i; h--) 
          {
	            child_index.var[h] = child_index.var[h-1]+1;
	        }
	    
    	    child_index.var[i-1] = child_index.var[i] - (x_curr - data + 1); 	    
	        j=i-1; 
	        i=MAXVAR;
	    }
    }
    // update the bits of the child
    child_index.var[MAXVAR] = 
        ((I.var[MAXVAR] << 1) & ((~((ushort)0)) << (j+1) )) | ((ushort) which << j) | ( I.var[MAXVAR] & ~((~((ushort)0)) << j) );
    return child_index;
}



int sum_elem(int *  my_array , int num_elem)
{
  int output = 0;
  for(int i=0; i < num_elem; i++)
  {
    output += my_array[i];
  }
  return output;
}


INDEX_TYPE get_next_node(INDEX_TYPE& I, int p, int level) 
{
    INDEX_TYPE node = I;
    int i = level-1; int j = p + level -2;
    while (i>=0 && node.var[i] == j+1) {i--;j--;}
    if (i < 0) 
    { //reach the end of nodes
        for (int h = 0; h <= MAXVAR; h++) 
        {
          node.var[h]=0; //invalid node
        }
    } 
    else 
    {
        node.var[i] += 1;
        for (j=i+1;j<level;j++) 
        {
            node.var[j] = node.var[i]+ j-i;
        }
    }
    node.var[MAXVAR] = 0;

    return node;
}



