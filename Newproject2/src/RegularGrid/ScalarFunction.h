#ifndef _SCALARFUNCTION_H_
#define _SCALARFUNCTION_H_

#include "Core/Vec.h"
#include <cmath>

// #define M_PI 3.141592653589793238

template <int Dim> class ScalarFunction;

template <int Dim>
class ScalarFunction{
public:
  virtual const Real operator()(const Vec<Real,Dim>& pt) const = 0;
};


class D2Func1f: public ScalarFunction<2>{
public:
  const Real operator()(const Vec<Real,2>& pt) const {
    Real x = pt[0];
    Real y = pt[1];
    return (x*x+y*y)*exp(x*y);
  }
};

class D2Func1F: public ScalarFunction<2>{
public:
  const Real operator()(const Vec<Real,2>& pt) const {
    Real x = pt[0];
    Real y = pt[1];
    return exp(x*y);
  }
};

class D2FuncTest: public ScalarFunction<2>{
public:
  const Real operator()(const Vec<Real,2>& pt) const {
    Real x = pt[0]*100000;
    Real y = pt[1]*10000;
    return x + y;
  }
};

class D2Func1Fx: public ScalarFunction<2>{
public:
  const Real operator()(const Vec<Real,2>& pt) const {
    Real x = pt[0];
    Real y = pt[1];
    return y*exp(x*y);
  }
};

class D2Func1Fy: public ScalarFunction<2>{
public:
  const Real operator()(const Vec<Real,2>& pt) const {
    Real x = pt[0];
    Real y = pt[1];
    return x*exp(x*y);
  }
};

class D2Func0: public ScalarFunction<2>{
public:
  const Real operator()(const Vec<Real,2>& pt) const {
    return 0.0;
  }
};




class D1Func1f: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
    Real x = pt[0];
    return -M_PI*M_PI*sin(M_PI*x);
  }
};

class D1Func1F: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
    Real x = pt[0];
    return sin(M_PI*x);
  }
};

class D1Func1Fx: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
     Real x = pt[0];
     return M_PI*cos(M_PI*x);
  }
};

class D1Func0: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
    return 0.0;
  }
};





#endif // _SCALARFUNCTION_H_
