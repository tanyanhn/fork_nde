#ifndef __MGGS_FACTORY__
#define __MGGS_FACTORY__
#include "Info.h"
#include "cal.h"

template<class T>
class TimeIntegratorFactory{
 public:
 static T& Instance();
 private:
 TimeIntegratorFactory();
 TimeIntegratorFactory(const TimeIntegratorFactory&);
 ~TimeIntegratorFactory();
 static TimeIntegratorFactory* pInstance_;
};


template<class CalPolicy>
class TimeIntegrator: public CalPolicy{
 private:
  Info_Table* ITable;
 public:
  TimeIntegrator();
  ~TimeIntegrator();
  Info_Table* get_table();
  Vector4d one_step(Vector4d* u, const double dt, const double mu, int _acc);
  Vector4d n_steps(double& time, Vector4d u0, const double dt, const double mu, const int _acc, int N);
  double err_Initial(double& time, Vector4d u0, const double dt, const double mu, const int _acc, const double T);
};



#else
//do nothing
#endif
