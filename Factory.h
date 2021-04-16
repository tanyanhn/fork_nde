#ifndef __MGGS_FACTORY__
#define __MGGS_FACTORY__
#include "Info.h"
#include "cal.h"

template<class T>
class TimeIntegratorFactory{
 public:
 TimeIntegratorFactory();
 ~TimeIntegratorFactory();
 T* get_pointer();
 private:
 T* pInstance_;
};


template<class CalPolicy, int acc>
class Method: public CalPolicy{
 private:
  Info_Table* ITable;
 public:
  Method();
  ~Method();
  Info_Table* get_table();
  Vector4d one_step(Vector4d* u, const double dt, const double mu);
  Vector4d n_steps(double& time, Vector4d u0, const double dt, const double mu, int N);
  double err_Initial(double& time, Vector4d u0, const double dt, const double mu, const double T);
  double err_Richardson(double tol, double& time, Vector4d u0, double dt, const double mu, int N);
  double Grid_Refine1(Vector4d u0, double dt, const double mu, const double T);
  double Grid_Refine2(double tol,Vector4d u0, double dt, const double mu, int N);
};


#else
//do nothing
#endif
