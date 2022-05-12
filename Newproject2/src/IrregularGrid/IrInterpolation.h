#ifndef _IRINTERPOLATION_H_
#define _IRINTERPOLATION_H_

#include "IrregularGrid/IrRectDomain.h"
#include "Core/Tensor.h"

template <int Dim> class IrInterpolation;

template <>
class IrInterpolation<2>{
public:
  using iVec = Vec<int,2>;
  
  IrInterpolation(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
                  acoarseDomain);
  virtual ~IrInterpolation() = default;
  virtual void applyInterpolate(const Tensor<Real,2>& data,
                                Tensor<Real,2>& res) = 0;
protected:
  IrRectDomain<2> fineDomain;
  IrRectDomain<2> coarseDomain;
};

IrInterpolation<2>::IrInterpolation(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
                                    acoarseDomain):fineDomain(afineDomain),coarseDomain(acoarseDomain){}

template <int Dim> class LinearInterpolation;

template <>
class LinearInterpolation<2>: public IrInterpolation<2>{
public:
  using BaseClass = IrInterpolation<2>;

  LinearInterpolation(const IrRectDomain<2>& afineDomain, const IrRectDomain<2>&
                      acoarseDomain):IrInterpolation(afineDomain,acoarseDomain){}
  
  void applyInterpolate(const Tensor<Real,2>& data,
                     Tensor<Real,2>& res){
    Box<2> bx = res.box();
    const Vec<int,2>& lo = bx.lo();
    const Vec<int,2>& hi = bx.hi();
    loop_box_2(bx,i,j){
      if (i != lo[1] && i != hi[1] && j != lo[0] && j != hi[0] && fineDomain.isInside(iVec{i,j})){
        if (i%2 == 0 && j%2 == 0)
          res(i,j) = data(i/2,j/2);
        else if (i%2 == 0 && j%2 != 0)
          res(i,j) = (data(i/2,j/2)+data(i/2,j/2+1))/2.0;
        else if (i%2 != 0 && j%2 == 0)
          res(i,j) = (data(i/2,j/2)+data(i/2+1,j/2))/2.0;
        else
          res(i,j) = (data(i/2,j/2)+data(i/2+1,j/2)+data(i/2,j/2+1)+data(i/2+1,j/2+1))/4.0;
      }
    }
  }
};


#endif // _IRINTERPOLATION_H_
