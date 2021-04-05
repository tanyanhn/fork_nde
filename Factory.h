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


template<template <class> class CalPolicy>
class TimeIntegrator{
 private:
  Info_Table ITable;
 public:
  TimeIntegrator();
  ~TimeIntegrator();
  Vector4d one_step(Vector4d* u, const double dt, const double mu, int _acc);
  Vector4d n_steps(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N);
  double err_initial(double& time, Vector4d u0, const double dt, const double mu, double T);
};
