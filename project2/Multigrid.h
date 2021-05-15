#ifndef __MGGS_MULTIGRID__
#define __MGGS_MULTIGRID__
#include <iostream>
#include <cmath>
#include <map>
#include "cblas.h"

class function{
 public:
  double action(double _x){
     double tmp1 = sin(_x);
     double tmp2 = cos(_x);
     double tmp3 = exp(tmp1);
     return tmp1*tmp3-tmp2*tmp2*tmp3;
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

class injection{
public:
  static double* action(double* _vec, int& n){
    if ( n < 8 || n%2 != 0){
      std::cerr << "No matching n!" << std::endl;
      return NULL;
    }
    else{
      double* vec = new double[n/2 - 1];
      for (int i = 0; i < n/2 - 1 ; i++)
	vec[i] = _vec[2*i+1];
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
    for (int i = 1; i < 2*n -2 ; i++){
      if (i%2 == 1)
	vec[i] = _vec[(i-1)/2];
      else
	vec[i] = 0.5*(_vec[i/2 - 1]+_vec[i/2]);
    }
    vec[2*n-2] = 0.5*_vec[n-2];
    n = n*2;
    return vec;
  }
};

class quadratic{
public:
  static double* action(double* _vec, int& n){
    double* vec = new double[2*n - 1];
    vec[0] = 1.5*_vec[0]-0.5*_vec[1];
    for (int i = 1; i < 2*n -4 ; i++){
      if (i%2 == 1)
	vec[i] = _vec[(i-1)/2];
      else
	vec[i] = 0.125*(15*_vec[i/2 - 1]-10*_vec[i/2]+3*_vec[i/2 + 1]);
    }
    vec[2*n-2] = 1.5*_vec[n-2]-0.5*_vec[n-3];
    vec[2*n-3] = _vec[n-2];
    vec[2*n-4] = 0.125*(15*_vec[n-2]-10*_vec[n-3]+3*_vec[n-4]);
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
  double* lefthand(int _n);
  double* righthand(int _n);
  double* weighted_Jacobi(double* _A, double* _f, double* _init, int _n, double weight, int times);
  //double* V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2);
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

template <class RestrictionPolicy, class InterpolationPolicy>
double* Multigrid<RestrictionPolicy,InterpolationPolicy>::lefthand(int _n){
  double h = 1.0/_n;
  double tmp1 = -1/(h*h);
  double tmp2 = -2*tmp1;
  double* A = new double[(_n-1)*(_n-1)];
  for (int i = 0 ; i < (_n-1)*(_n-1) ; i++)
    A[i] = 0;
  A[0] = tmp2;A[1] = tmp1;
  for (int i = 1 ; i < _n - 2 ; i++){
    A[i*_n-1] = tmp1;
    A[i*_n] = tmp2;
    A[i*_n+1] = tmp1;
  }
  A[(_n-1)*(_n-1)-2] = tmp1;
  A[(_n-1)*(_n-1)-1] = tmp2;
  return A;
}

template <class RestrictionPolicy, class InterpolationPolicy>
double* Multigrid<RestrictionPolicy,InterpolationPolicy>::righthand(int _n){
  double h = 1.0/_n;
  double tmp1 = boundary.first/(h*h);
  double tmp2 = boundary.second/(h*h);
  double* f = new double[_n-1];
  f[0] = u->action(h)+tmp1;
  for (int i = 1 ; i < _n-2 ; i++)
    f[i] = u->action((i+1)*h);
  f[_n-2] = u->action((_n-1)*h)+tmp2;
  return f;
}


template <class RestrictionPolicy, class InterpolationPolicy>
double* Multigrid<RestrictionPolicy,InterpolationPolicy>::weighted_Jacobi(double* _A, double* _f, double* _init, int _n, double weight, int times)
{
  double h = 1.0/_n;
  double tmp = 0.5*weight*h*h;
  double* T = new double[(_n-1)*(_n-1)];
  for (int i = 0 ; i < _n - 1 ; i++){
    for (int j = 0 ; j < _n - 1 ; j++){
      if (i == j)
	T[i*(_n-1)+j] = 1-tmp*_A[i*(_n-1)+j];
      else
	T[i*(_n-1)+j] = -tmp*_A[i*(_n-1)+j];
    }
  }
  double* cc = new double[_n-1];
  for (int i = 0 ; i < _n-1 ; i++)
    cc[i] = tmp*_f[i];
  double* c = new double[_n-1];
  double* result = _init;
  int count = 0; 
  while (count < times){
    cblas_dgemv(CblasRowMajor,CblasNoTrans,_n-1,_n-1,1.0,T,_n-1,result,1,0,c,1);
    for ( int i = 0 ; i < _n-1 ; i++)
      result[i]=c[i]+ cc[i];
    count++;
  }
  return result;
}

#else
//do nothing
#endif
