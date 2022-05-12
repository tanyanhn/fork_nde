#ifndef _IRRECTDOMAIN_H_
#define _IRRECTDOMAIN_H_

#include "RegularGrid/RectDomain.h"
#include "RegularGrid/Func.h"
#include "Core/Tensor.h"

#include <vector>

template <int Dim> class IrRectDomain;

template <>
class IrRectDomain<2>{
public:

  using iVec = Vec<int,2>;

  using rVec = Vec<Real,2>;

  template <class T>
  using Vector = std::vector<T>;

  IrRectDomain() = default;

  IrRectDomain(const RectDomain<2>& aDomain, const Func<1>* paFunc, const
               Real atol);

  const RectDomain<2>& getDomain() const;

  const Func<1>* getplbFunc() const;

  bool isInside(const iVec& pt) const;

  bool isRePoint(const iVec& pt) const;

  Vec<Real,2> getProjPoint(const iVec& pt) const;
protected:
  bool isBelowCurve(const rVec& pt) const;
protected:
  RectDomain<2> Domain;
  const Func<1>* plbFunc;
  Real tol;
};

IrRectDomain<2>::IrRectDomain(const RectDomain<2>& aDomain, const Func<1>*
                              paFunc, const Real
                              atol):Domain(aDomain),plbFunc(paFunc),tol(atol){
  assert(aDomain.getCentering() == NodeCentered);
  assert(aDomain.getNumGhost() == 1);
}

const RectDomain<2>& IrRectDomain<2>::getDomain() const{
  return Domain;
}

const Func<1>* IrRectDomain<2>::getplbFunc() const{
  return plbFunc;
}

bool IrRectDomain<2>::isInside(const iVec& pt) const{
  if (!Domain.contain(pt))
    return 0;
  const iVec& lo = Domain.lo();
  //const iVec& hi = Domain.hi();
  const rVec& dx = Domain.spacing();
  rVec rpt = (pt - lo) * dx;//set 0.0 by default
  if (isBelowCurve(rpt))
    return 0;
  return 1;
}

bool IrRectDomain<2>::isRePoint(const iVec& pt) const{
  const iVec& lo = Domain.lo();
  const rVec& dx = Domain.spacing();
  if (isBelowCurve((pt - iVec{1,0} - lo)*dx) || isBelowCurve((pt +
                                                             iVec{1,0} - lo)*dx) ||
      isBelowCurve((pt - iVec{0,1} - lo)*dx) || isBelowCurve((pt +
                                                              iVec{0,1} - lo)*dx))
    return 0;
  return 1;
}

bool IrRectDomain<2>::isBelowCurve(const rVec& pt) const{
  Real y = (*plbFunc)(pt[0]);
  if (y - pt[1] > tol)
    return 1;
  return 0;
}

Vec<Real,2> IrRectDomain<2>::getProjPoint(const iVec& pt) const{
  const int num = 40;
  const iVec& lo = Domain.lo();
  const rVec& dx = Domain.spacing();
  rVec rpt = (pt - lo) * dx;
  Real x1 = rpt[0]-dx[0];
  Real x2 = rpt[0]+dx[0];
  Real h = (x2 - x1)/num;
  Real Dot = 10;
  int flag = 0;
  for (int i = 0 ; i <= num ; i++){
    Real x = x1 + h*i;
    rVec v1 = rpt - rVec{x,(*plbFunc)(x)};
    rVec v2 = rVec{x,(*plbFunc).Fx(x)};
    Real tmp = dot(normalize(v1),normalize(v2));
    if (fabs(tmp) < Dot){
      Dot = fabs(tmp);
      flag = i;
    }
  }
  Real x = x1 + h*flag;
  rVec respt = rVec{x,(*plbFunc)(x)};
  rVec frac = rVec{(respt-rpt)/dx};
  return frac;
}



#endif // _IRRECTDOMAIN_H_
