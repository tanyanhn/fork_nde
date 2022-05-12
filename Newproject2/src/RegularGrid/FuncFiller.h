#ifndef _FUNCFILLER_H_
#define _FUNCFILLER_H_

#include "RegularGrid/ScalarFunction.h"
#include "RegularGrid/RectDomain.h"
#include "Core/Tensor.h"

template <int Dim> class FuncFiller;

template <int Dim>
class FuncFiller{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;
  
  FuncFiller(const RectDomain<Dim>& adomain);

  void fill(Tensor<Real,Dim>& res, const ScalarFunction<Dim>* pfunc) const;

protected:
  RectDomain<Dim> domain;
};

template <int Dim>
FuncFiller<Dim>::FuncFiller(const RectDomain<Dim>& adomain){
  assert(adomain.getCentering() == NodeCentered);
  domain = adomain;
}

template <>
void FuncFiller<2>::fill(Tensor<Real,2>& res, const ScalarFunction<2>*
                         pfunc) const{
  const Box<2>& bx = domain;
  const Vec<int,2>& lo = bx.lo();
  const Vec<Real,2>& dx = domain.spacing();
  loop_box_2(bx,i,j){
    iVec Node{i,j};
    rVec rNode = (Node-lo)*dx;
    res(i,j) = (*pfunc)(rNode);
  }
}

template <>
void FuncFiller<1>::fill(Tensor<Real,1>& res, const ScalarFunction<1>*
                         pfunc) const{
  const Box<1>& bx = domain;
  const Vec<int,1>& lo = bx.lo();
  const Vec<Real,1>& dx = domain.spacing();
  loop_box_1(bx,i){
    Vec<int,1> Node{i};
    Vec<Real,1> rNode = (Node-lo)*dx;
    res(i) = (*pfunc)(rNode);
  }
}

#endif //_FUNCFILLER_H_
