#ifndef _SCALARFUNCTION_H_
#define _SCALARFUNCTION_H_

#include "Core/Vec.h"
#include <cmath>

#define PI 3.141592653589793238

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
    return -PI*PI*sin(PI*x);
  }
};

class D1Func1F: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
    Real x = pt[0];
    return sin(PI*x);
  }
};

class D1Func1Fx: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
     Real x = pt[0];
     return PI*cos(PI*x);
  }
};

class D1Func0: public ScalarFunction<1>{
public:
  const Real operator()(const Vec<Real,1>& pt) const {
    return 0.0;
  }
};





#endif // _SCALARFUNCTION_H_
