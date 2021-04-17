#ifndef __MGGS_FACTORY__
#define __MGGS_FACTORY__
#include <map>
#include "Info.h"
#include "cal.h"

class TimeIntegrator{
 public:
  TimeIntegrator();
  ~TimeIntegrator();
};


template<class CalPolicy, int acc>
class Method: public CalPolicy, public TimeIntegrator{
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


class TimeIntegratorFactory{
 public :
  typedef TimeIntegrator* (*IntegratorCreator)();
  static TimeIntegratorFactory* Instance();
  
 public:
  bool RegisterIntegrator(const std::pair<const std::string,int> &MethodId, IntegratorCreator createFn);
  bool UnregisterIntegrator(const std::pair<const std::string,int> &MethodId);
  TimeIntegrator* CreateIntegrator(const std::pair<const std::string,int> &MethodId) const;

 private:
  TimeIntegratorFactory();
  TimeIntegratorFactory(const TimeIntegratorFactory&);
  static TimeIntegratorFactory* pInstance_;

 private:
  typedef std::map<const std::pair<const std::string,int>, IntegratorCreator> CallbackMap;
  CallbackMap callbacks_;
};


#else
//do nothing
#endif
