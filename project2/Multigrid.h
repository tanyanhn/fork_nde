#ifndef __MGGS_MULTIGRID__
#define __MGGS_MULTIGRID__
#include <iostream>
#include <cmath>
#include <map>

class function{
 public:
  double f(double _x){
     double tmp1 = sin(_x);
     double tmp2 = cos(_x);
     double tmp3 = exp(tmp1);
     return -tmp1*tmp3+tmp2*tmp2*tmp3;
  }
};

class full_weighting{
public:
  static double* action(double* _vec, int& n){
    if ( n < 8 || n%2 != 0){
      std::cerr << "No matching n!" << std::endl;
      return NULL;
    }
    else{
      double* vec = new double[n/2 - 1];
      for (int i = 0; i < n/2 - 1 ; i++)
	vec[i] = 0.25*(_vec[2*i]+2*_vec[2*i+1]+_vec[2*i+2]);
      n = n/2;
      return vec;
    }  
  }
};

class linear{
public:
  static double* action(double* _vec, int& n){
    double* vec = new double[2*n - 1];
    vec[0] = 0.5*_vec[0];
    vec[2*n-2] = 0.5*_vec[n-2];
    for (int i = 1; i < 2*n -2 ; i++){
      if (i%2 == 1)
	vec[i] = _vec[(i-1)/2];
      else
	vec[i] = 0.5*(_vec[i/2 - 1]+_vec[i/2]);
    }
    n = n*2;
    return vec;
  }
};

template <class RestrictionPolicy, class InterpolationPolicy>
class Multigrid: public RestrictionPolicy, public InterpolationPolicy{
private:
  function* u;
  std::pair<double,double> boundary;
  double* initial;
  std::pair<int,double> criteria;
public:
  Multigrid();
  Multigrid(std::pair<double,double> _boundary, double* _initial, std::pair<int,double> _criteria);
  ~Multigrid();
  void load_boundary(std::pair<double,double> _boundary);
  std::pair<double,double> get_boundary() const;
  void load_initial(double* _initial);
  double* get_initial() const;
  void load_criteria(std::pair<int,double> _criteria);
  std::pair<int,double> get_criteria() const;
  void test_restriction(int &n);
  void test_interpolation(int &n);
};


template <class RestrictionPolicy, class InterpolationPolicy>
Multigrid<RestrictionPolicy,InterpolationPolicy>::Multigrid(){
  u = new function;
  initial = NULL;
}

template <class RestrictionPolicy, class InterpolationPolicy>
Multigrid<RestrictionPolicy,InterpolationPolicy>::Multigrid(std::pair<double,double> _boundary, double* _initial, std::pair<int,double> _criteria){
  u = new function;
  boundary = _boundary;
  initial = _initial;
  criteria = _criteria;
}

template <class RestrictionPolicy, class InterpolationPolicy>
Multigrid<RestrictionPolicy,InterpolationPolicy>::~Multigrid(){
  delete u;
}

template <class RestrictionPolicy, class InterpolationPolicy>
void Multigrid<RestrictionPolicy,InterpolationPolicy>::load_boundary(std::pair<double,double> _boundary){
  boundary = _boundary;
}

template <class RestrictionPolicy, class InterpolationPolicy>
std::pair<double,double> Multigrid<RestrictionPolicy,InterpolationPolicy>::get_boundary() const{
  return boundary;
}

template <class RestrictionPolicy, class InterpolationPolicy>
void Multigrid<RestrictionPolicy,InterpolationPolicy>::load_initial(double* _initial){
  initial = _initial;
}

template <class RestrictionPolicy, class InterpolationPolicy>
double* Multigrid<RestrictionPolicy,InterpolationPolicy>::get_initial() const{
  return initial;
}

template <class RestrictionPolicy, class InterpolationPolicy>
void Multigrid<RestrictionPolicy,InterpolationPolicy>::load_criteria(std::pair<int,double> _criteria){
  criteria = _criteria;
}

template <class RestrictionPolicy, class InterpolationPolicy>
std::pair<int,double> Multigrid<RestrictionPolicy,InterpolationPolicy>::get_criteria() const{
  return criteria;
}


template <class RestrictionPolicy, class InterpolationPolicy>
void Multigrid<RestrictionPolicy,InterpolationPolicy>::test_restriction(int &n){
  initial = RestrictionPolicy().action(initial,n);
}

template <class RestrictionPolicy, class InterpolationPolicy>
void Multigrid<RestrictionPolicy,InterpolationPolicy>::test_interpolation(int &n){
  initial = InterpolationPolicy().action(initial,n);
}








#else
//do nothing
#endif
