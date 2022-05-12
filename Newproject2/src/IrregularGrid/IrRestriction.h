#ifndef _IRRESTRICTION_H_
#define _IRRESTRICTION_H_

#include "IrregularGrid/IrRectDomain.h"
#include "Core/Tensor.h"

template <int Dim> class IrRestriction;

template <>
class IrRestriction<2>{
public:
  using iVec = Vec<int,2>;
  
  IrRestriction(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
  acoarseDomain);
  virtual ~IrRestriction() = default;
  virtual void applyRestrict(const Tensor<Real,2>& data,
                             Tensor<Real,2>& res) = 0;
protected:
  IrRectDomain<2> fineDomain;
  IrRectDomain<2> coarseDomain;
};

IrRestriction<2>::IrRestriction(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
                                acoarseDomain):fineDomain(afineDomain),coarseDomain(acoarseDomain){}

template <int Dim> class Injection;

template <>
class Injection<2>: public IrRestriction<2>{
public:
  using BaseClass = IrRestriction<2>;

  Injection(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
            acoarseDomain):IrRestriction(afineDomain,acoarseDomain){}
  
  void applyRestrict(const Tensor<Real,2>& data,
                     Tensor<Real,2>& res){
    Box<2> bx = res.box();
    const Vec<int,2>& lo = bx.lo();
    const Vec<int,2>& hi = bx.hi();
    loop_box_2(bx,i,j){
      if (i != lo[1] && i != hi[1] && j != lo[0] && j != hi[0] && coarseDomain.isInside(iVec{i,j}))
        res(i,j) = data(2*i,2*j);
    }
  }
};

template <int Dim> class FullWeightingRestriction;

template <>
class FullWeightingRestriction<2>: public IrRestriction<2>{
public:
  using BaseClass = IrRestriction<2>;

  FullWeightingRestriction(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
                           acoarseDomain):IrRestriction(afineDomain,acoarseDomain){}
  
  void applyRestrict(const Tensor<Real,2>& data,
                     Tensor<Real,2>& res){
    Box<2> bx = res.box();
    const Vec<int,2>& lo = bx.lo();
    const Vec<int,2>& hi = bx.hi();
    loop_box_2(bx,i,j){
      if (i != lo[1] && i != hi[1] && j != lo[0] && j != hi[0] && coarseDomain.isInside(iVec{i,j})){
        iVec Node{2*i,2*j};
        if (!fineDomain.isInside(Node+iVec{0,-1})){
          if (!fineDomain.isInside(Node+iVec{-1,0})){
            res(i,j) = (4*data(2*i,2*j)+4*data(2*i+1,2*j)+4*data(2*i,2*j+1)+4*data(2*i+1,2*j+1))/16.0;
          }
          else if (!fineDomain.isInside(Node+iVec{1,0})){
            res(i,j) = (4*data(2*i,2*j)+4*data(2*i-1,2*j)+4*data(2*i,2*j+1)+4*data(2*i-1,2*j+1))/16.0;
          }
          else{
            res(i,j) = (4*data(2*i,2*j)+2*data(2*i+1,2*j)+4*data(2*i,2*j+1)+2*data(2*i-1,2*j)+2*data(2*i+1,2*j+1)+2*data(2*i-1,2*j+1))/16.0;
          }
        }
        else if (!fineDomain.isInside(Node+iVec{1,0})){
          res(i,j) = (4*data(2*i,2*j)+4*data(2*i-1,2*j)+2*data(2*i,2*j+1)+2*data(2*i-1,2*j+1)+2*data(2*i,2*j-1)+2*data(2*i-1,2*j-1))/16.0;
        }
        else if (!fineDomain.isInside(Node+iVec{-1,0})){
          res(i,j) = (4*data(2*i,2*j)+4*data(2*i+1,2*j)+2*data(2*i,2*j+1)+2*data(2*i+1,2*j+1)+2*data(2*i,2*j-1)+2*data(2*i+1,2*j-1))/16.0;
        }
        else{
          res(i,j) = (4*data(2*i,2*j)+2*data(2*i+1,2*j)+2*data(2*i,2*j+1)+2*data(2*i-1,2*j)+2*data(2*i,2*j-1)+data(2*i+1,2*j+1)+data(2*i-1,2*j+1)+data(2*i+1,2*j-1)+data(2*i-1,2*j-1))/16.0;
        }
      }  
    }  
  }
};



#endif // _IRRESTRICTION_H_
