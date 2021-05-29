#ifndef __MGGS_MULTIGRID__
#define __MGGS_MULTIGRID__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "cblas.h"

class function1{
 public:
  double ref(double _x){
    return exp(sin(_x));
  }
  
  double action(double _x){
     double tmp1 = sin(_x);
     double tmp2 = cos(_x);
     double tmp3 = exp(tmp1);
     return tmp1*tmp3-tmp2*tmp2*tmp3;
  }
};

class function2{
  public:
  double ref(double _x){
    return sin(M_PI*_x);
  }
  
  double action(double _x){
     return M_PI*M_PI*sin(M_PI*_x);
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
    vec[0] = 0.0625*(5*_vec[0]+3*_vec[1]);
    vec[1] = _vec[0];
    vec[2] = 0.0625*(5*_vec[0]+5*_vec[1]+3*_vec[2]);
    for (int i = 3; i < 2*n -4 ; i++){
      if (i%2 == 1)
	vec[i] = _vec[(i-1)/2];
      else
	vec[i] = 0.0625*(3*_vec[i/2 - 2]+5*_vec[i/2 - 1]+5*_vec[i/2]+3*_vec[i/2 + 1]);
    }
    vec[2*n-4] = 0.0625*(5*_vec[n-2]+5*_vec[n-3]+3*_vec[n-4]);
    vec[2*n-3] = _vec[n-2];
    vec[2*n-2] = 0.0625*(5*_vec[n-2]+3*_vec[n-3]);
    n = n*2;
    return vec;
  }
};

class Grid{
 protected:
  std::pair<double,double> boundary;
  double* initial;
  std::pair<int,double> criteria;
 public:
  Grid();
  virtual void load_boundary(std::pair<double,double> _boundary) = 0;
  virtual std::pair<double,double> get_boundary() const = 0;
  virtual void load_initial(double* _initial) = 0;
  virtual double* get_initial() const = 0;
  virtual void load_criteria(std::pair<int,double> _criteria) = 0;
  virtual std::pair<int,double> get_criteria() const = 0;
  virtual double* righthand(int _n) = 0;
  virtual double* V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2) = 0;
  virtual double* fm_cycle(int &_n, double* _f, int _t1, int _t2) = 0;
  virtual double analysis_V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2) = 0;
  virtual double analysis_fm_cycle(int &_n, double* _f, int _t1, int _t2) = 0;
  virtual void show_result(double* _result, int _n, const std::string str) = 0;
};

Grid::Grid(){};


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
class Multigrid:public Grid, public RestrictionPolicy, public InterpolationPolicy, public Function{
public:
  Multigrid();
  Multigrid(std::pair<double,double> _boundary, double* _initial, std::pair<int,double> _criteria);
  ~Multigrid();
  virtual void load_boundary(std::pair<double,double> _boundary);
  virtual std::pair<double,double> get_boundary() const;
  virtual void load_initial(double* _initial);
  virtual double* get_initial() const;
  virtual void load_criteria(std::pair<int,double> _criteria);
  virtual std::pair<int,double> get_criteria() const;
  void test_restriction(int &n);
  void test_interpolation(int &n);
  double* lefthand(int _n);
  virtual double* righthand(int _n);
  double* weighted_Jacobi(double* _A, double* _f, double* _init, int _n, double weight, int times);
  double* residual(double* _A, double* _f, double* _u, int _n);
  double* ref_solution(int _n);
  double max_norm(double* _u, int _n);
  double two_norm(double* _u, int _n);
  double* onestep_V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2);
  double* onestepr_V_cycle(double* _v, double* _r, int& _n, double* _f, int _t1, int _t2);
  int n_iteration_V_cycle(int _n, int _t1, int _t2);
  virtual double* V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2);
  virtual double* fm_cycle(int &_n, double* _f, int _t1, int _t2);
  int n_iteration_fm_cycle(int _n, int _t1, int _t2);
  virtual double analysis_V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2);
  virtual double analysis_fm_cycle(int &_n, double* _f, int _t1, int _t2);
  virtual void show_result(double* _result, int _n, const std::string str);
};


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
  Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::Multigrid(){
  initial = NULL;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::Multigrid(std::pair<double,double> _boundary, double* _initial, std::pair<int,double> _criteria){
  boundary = _boundary;
  initial = _initial;
  criteria = _criteria;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::~Multigrid(){
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
void Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::load_boundary(std::pair<double,double> _boundary){
  boundary = _boundary;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
std::pair<double,double> Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::get_boundary() const{
  return boundary;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
void Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::load_initial(double* _initial){
  initial = _initial;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::get_initial() const{
  return initial;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
void Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::load_criteria(std::pair<int,double> _criteria){
  criteria = _criteria;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
std::pair<int,double> Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::get_criteria() const{
  return criteria;
}


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
void Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::test_restriction(int &n){
  initial = RestrictionPolicy().action(initial,n);
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
void Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::test_interpolation(int &n){
  initial = InterpolationPolicy().action(initial,n);
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::lefthand(int _n){
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

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::righthand(int _n){
  double h = 1.0/_n;
  double tmp1 = boundary.first/(h*h);
  double tmp2 = boundary.second/(h*h);
  double* f = new double[_n-1];
  f[0] = Function().action(h)+tmp1;
  for (int i = 1 ; i < _n-2 ; i++)
    f[i] = Function().action((i+1)*h);
  f[_n-2] = Function().action((_n-1)*h)+tmp2;
  return f;
}


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::weighted_Jacobi(double* _A, double* _f, double* _init, int _n, double weight, int times)
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
  double* result = new double[_n-1];
  for (int i = 0 ; i < _n-1 ; i++)
    result[i] = _init[i];
  int count = 0; 
  while (count < times){
    cblas_dgemv(CblasRowMajor,CblasNoTrans,_n-1,_n-1,1.0,T,_n-1,result,1,0,c,1);
    for ( int i = 0 ; i < _n-1 ; i++)
      result[i]=c[i]+ cc[i];
    count++;
  }
  return result;
}


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::residual(double* _A, double* _f, double* _u, int _n){
  double* result = new double[_n-1];
  cblas_dgemv(CblasRowMajor,CblasNoTrans,_n-1,_n-1,-1.0,_A,_n-1,_u,1,0,result,1);
  for (int i = 0 ; i < _n-1 ; i++)
    result[i] += _f[i];
  return result;
}
  

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::ref_solution(int _n){
  double h = 1.0/_n;
  double* solution = new double[_n-1];
  for (int i = 0 ; i < _n-1 ; i++)
    solution[i] = Function().ref((i+1)*h);
  return solution;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::max_norm(double* _u, int _n){
  int index = (int)cblas_idamax(_n-1,_u,1);
  return fabs(_u[index]);
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::two_norm(double* _u, int _n){
  double result = cblas_dnrm2(_n-1,_u,1);
  return result;
}


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::onestep_V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2){
  double weight = 2.0/3;
  double* A = this->lefthand(_n);
  double* v = weighted_Jacobi(A,_f,_v,_n,weight,_t1);
  if (_n > 4){
    double* f2 = RestrictionPolicy().action(this->residual(A,_f,v,_n),_n);
    double* _v2 = new double[_n-1];
    for (int i = 0 ; i < _n-1 ; i++)
      _v2[i] = 0;
    double* v2 = this->onestep_V_cycle(_v2,_n,f2,_t1,_t2);
    double* v2c = InterpolationPolicy().action(v2,_n);
    for (int i = 0 ; i < _n-1 ; i++)
      v[i] += v2c[i];
  }
  double* result = weighted_Jacobi(A,_f,v,_n,weight,_t2);
  return result;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::onestepr_V_cycle(double* _v, double* _r, int& _n, double* _f, int _t1, int _t2){
  double weight = 2.0/3;
  double* A = this->lefthand(_n);
  double* f2 = RestrictionPolicy().action(_r,_n);
  double* _v2 = new double[_n-1];
  for (int i = 0 ; i < _n-1 ; i++)
    _v2[i] = 0;
  double* v2 = this->onestep_V_cycle(_v2,_n,f2,_t1,_t2);
  double* v2c = InterpolationPolicy().action(v2,_n);
  for (int i = 0 ; i < _n-1 ; i++)
    _v[i] += v2c[i];
  double* result = weighted_Jacobi(A,_f,_v,_n,weight,_t2);
  return result;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
int Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::n_iteration_V_cycle(int _n, int _t1, int _t2){
  int max = 0;
  while(_n >= 4){
    max += _t1+_t2;
    _n = _n/2;
  }
  return max;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2){
  if (criteria.first == 0){
    int MAX = (int)(criteria.second + 0.01);
    int count = 1;
    double* A = this->lefthand(_n);
    _v = this->onestep_V_cycle(_v,_n,_f,_t1,_t2);
    double* res = this->residual(A,_f,_v,_n);
    while (count < MAX){
      _v = this->onestepr_V_cycle(_v,res,_n,_f,_t1,_t2);
      count++;
      res = this->residual(A,_f,_v,_n);
    }
    return _v;
  }
  else if (criteria.first == 1){
    double tol = criteria.second;;
    double* A = this->lefthand(_n);
    double f_norm = this ->two_norm(_f,_n);
    int count = 1;
    _v = this->onestep_V_cycle(_v,_n,_f,_t1,_t2);
    double* res = this->residual(A,_f,_v,_n);
    double res_norm = this->two_norm(res,_n);
    while (res_norm/f_norm > tol && count < 20){
      _v = this->onestepr_V_cycle(_v,res,_n,_f,_t1,_t2);
      count++;
      res = this->residual(A,_f,_v,_n);
      res_norm = this->two_norm(res,_n);
    }
    return _v;
  }
  else{
    std::cerr << "Wrong criteria!" << std::endl;
    return NULL;
  }
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double* Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::fm_cycle(int &_n, double* _f, int _t1, int _t2){
  double* v;
  if (_n > 4){
    double* f2 = RestrictionPolicy().action(_f,_n);
    double* v2 = this->fm_cycle(_n,f2,_t1,_t2);
    v = InterpolationPolicy().action(v2,_n);
  }
  else{
    v = new double[_n-1];
    for (int i = 0 ; i < _n-1 ; i++)
      v[i] = 0;
  }
  double* result = this->V_cycle(v,_n,_f,_t1,_t2);
  return result;
}


template <class RestrictionPolicy, class InterpolationPolicy, class Function>
int Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::n_iteration_fm_cycle(int _n, int _t1, int _t2){
  int max = 0;
  while(_n >= 4){
    max += this->n_iteration_V_cycle(_n,_t1,_t2);
    _n = _n/2;
  }
  return max;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::analysis_V_cycle(double* _v, int& _n, double* _f, int _t1, int _t2){
  std::ofstream os;
  std::string ss = "V_cycle_" + std::to_string(_n) + ".m";
  const char *s = ss.c_str();
  os.open(s);
  double f_norm = this ->two_norm(_f,_n);
  double res_norm;
  os << "res=[\n";
  if (criteria.first == 0){
    int MAX = (int)(criteria.second + 0.01);
    int count = 1;
    double* A = this->lefthand(_n);
    _v = this->onestep_V_cycle(_v,_n,_f,_t1,_t2);
    double* res = this->residual(A,_f,_v,_n);
    res_norm = this->two_norm(res,_n);
    std::cout << "After " << count << "V-cycle, residual is:" << std::endl;
    for (int i = 0 ; i < _n-1 ; i++)
      std::cout << res[i] << std::endl;
    os << res_norm << ",\n";
    while (count < MAX){
      _v = this->onestepr_V_cycle(_v,res,_n,_f,_t1,_t2);
      count++;
      res = this->residual(A,_f,_v,_n);
      res_norm = this->two_norm(res,_n);
      std::cout << "After " << count << "V-cycle, residual is:" << std::endl;
      for (int i = 0 ; i < _n-1 ; i++)
	std::cout << res[i] << std::endl;
      os << res_norm << ",\n";
    }
    double* ref_solution = this->ref_solution(_n);
    double* err_vector = new double[_n-1];
    for (int i = 0 ; i < _n-1 ; i++)
      err_vector[i] = ref_solution[i] - _v[i];
    double err = this->max_norm(err_vector, _n);
    os << "];\n";
    os << "[n,m] = size(res);\n";
    os << "t = res(2:n)./res(1:n-1);\n";
    os << "rate = sum(t)/(n-1);\n";
    os << "disp(['Reduction rate of residuals =',num2str(rate)]);\n";
    os << "semilogy(1:" << count << ",res,'*-');\n";
    os.close();
    return err;
  }
  else if (criteria.first == 1){
    double tol = criteria.second;
    double* ref_solution = this->ref_solution(_n);
    //double ref_norm = this ->two_norm(ref_solution,_n);
    double* A = this->lefthand(_n);
    _v = this->onestep_V_cycle(_v,_n,_f,_t1,_t2);
    double* res = this->residual(A,_f,_v,_n);
    res_norm = this->two_norm(res,_n);
    double* err_vector = new double[_n-1];
    int count = 1;
    std::cout << "After " << count << " times V-cycle, residual is:" << std::endl;
    for (int i = 0 ; i < _n-1 ; i++)
      std::cout << res[i] << std::endl;
    os << res_norm << ",\n";
    while (res_norm/f_norm > tol && count < 20){
      _v = this->onestepr_V_cycle(_v,res,_n,_f,_t1,_t2);
      count++;
      res = this->residual(A,_f,_v,_n);
      res_norm = this->two_norm(res,_n);
      std::cout << "After " << count << " times V-cycle, residual is:" << std::endl;
      for (int i = 0 ; i < _n-1 ; i++)
        std::cout << res[i] << std::endl;
      os << res_norm << ",\n";
    }
    for (int i = 0 ; i < _n-1 ; i++)
	err_vector[i] = ref_solution[i] - _v[i];
    double err = this->max_norm(err_vector, _n);
    os << "];\n";
    os << "[n,m] = size(res);\n";
    os << "t = res(2:n)./res(1:n-1);\n";
    os << "rate = sum(t)/(n-1);\n";
    os << "disp(['Reduction rate of residuals =',num2str(rate)]);\n";
    os << "semilogy(1:" << count << ",res,'*-');\n";
    os.close();
    if (count == 20)
      std::cout << "Fail to achieve the preset accuracy." << std::endl;
    return err;
  }
  else{
    std::cerr << "Wrong criteria!" << std::endl;
    return 0;
  }
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
double Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::analysis_fm_cycle(int &_n, double* _f, int _t1, int _t2){
    double* result;
    result = this->fm_cycle(_n,_f,_t1,_t2);
    double* ref_solution = this->ref_solution(_n);
    double* err_vector = new double[_n-1];
    for (int i = 0 ; i < _n-1 ; i++)
      err_vector[i] = ref_solution[i] - result[i];
    double err = this->max_norm(err_vector, _n);
    return err;
}

template <class RestrictionPolicy, class InterpolationPolicy, class Function>
void Multigrid<RestrictionPolicy,InterpolationPolicy,Function>::show_result(double* _result, int _n, const std::string str){
  double* sol = this->ref_solution(_n);
  std::ofstream os;
  std::string ss = str + "_show_" + std::to_string(_n) + ".m";
  const char *s = ss.c_str();
  os.open(s);
  os << "y = [\n";
  for (int i = 0 ; i < _n-1 ; i++)
    os << _result[i] << "," << sol[i] <<";\n";
  os << "];\n";
  os << "x = 1:" << _n-1 << ";\n";
  os << "x = x/" << _n << ";\n";
  os << "result = y(:,1);\n";
  os << "solution = y(:,2);\n";
  os << "hold on;\n";
  os << "plot(x,result,'*');\n";
  os << "plot(x,solution,'-');\n";
  os << "legend('result','exact sulution');\n";
  os.close();
}

#else
//do nothing
#endif
