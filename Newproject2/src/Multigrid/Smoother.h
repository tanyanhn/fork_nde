#ifndef _SMOOTHER_H_
#define _SMOOTHER_H_

#include "RegularGrid/RectDomain.h"
#include "Core/Tensor.h"
#include "Core/TensorSlice.h"

template <int Dim> class Smoother;

template <int Dim>
class Smoother{
public:
  Smoother(const RectDomain<Dim>& adomain):domain(adomain){}

  virtual void apply(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                     rhs, Tensor<Real,Dim>& res) const = 0; 
protected:
  RectDomain<Dim> domain;
};



#endif // _SMOOTHER_H_
