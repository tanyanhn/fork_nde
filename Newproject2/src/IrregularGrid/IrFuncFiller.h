#ifndef _IRFUNCFILLER_H_
#define _IRFUNCFILLER_H_

#include "IrregularGrid/IrRectDomain.h"
#include "Core/Tensor.h"

template <int Dim> class IrFuncFiller;

template <>
class IrFuncFiller<2>{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;
  
  IrFuncFiller(const IrRectDomain<2>& aDomain);

  template <class TFunc>
  void fillFunc(Tensor<Real,2>& res, const TFunc& Func);

protected:
  IrRectDomain<2> Domain;
};


IrFuncFiller<2>::IrFuncFiller(const IrRectDomain<2>& aDomain):Domain(aDomain){}

template <class TFunc>
void IrFuncFiller<2>::fillFunc(Tensor<Real,2>& res, const TFunc&
                               Func){
  const RectDomain<2>& rDomain = Domain.getDomain();
  iVec lo = rDomain.lo();
  const rVec& dx = rDomain.spacing();
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (Domain.isInside(Node)){
      rVec rNode = (Node-lo)*dx;
      res(i,j) = Func(rNode);
    }
  }
}


#endif //_IRFUNCFILLER_H_
