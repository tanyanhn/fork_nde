#ifndef _RESTRICTOR_H_
#define _RESTRICTOR_H_

#include "RegularGrid/RectDomain.h"
#include "Core/Tensor.h"

template <int Dim> class Restrictor;

template <int Dim>
class Restrictor{
public:
  Restrictor(const RectDomain<Dim>& afineDomain, const RectDomain<Dim>&
  acoarseDomain);
  virtual ~Restrictor() = default;
  virtual void apply(const Tensor<Real,Dim>& data,
                     Tensor<Real,Dim>& res) const = 0;
protected:
  RectDomain<Dim> fineDomain;
  RectDomain<Dim> coarseDomain;
};

template <int Dim>
Restrictor<Dim>::Restrictor(const RectDomain<Dim>& afineDomain, const RectDomain<Dim>&
                              acoarseDomain){
  assert(afineDomain.getCentering() == NodeCentered &&
  acoarseDomain.getCentering() == NodeCentered);
  fineDomain = afineDomain;
  coarseDomain = acoarseDomain;
}

template <int Dim> class Injection;

template <>
class Injection<2>: public Restrictor<2>{
public:
  using BaseClass = Restrictor<2>;

  Injection(const RectDomain<2>& afineDomain, const RectDomain<2>&
            acoarseDomain):Restrictor(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,2>& data,
             Tensor<Real,2>& res) const{
    auto& bx = coarseDomain;
    const Vec<int,2>& lo = bx.lo();
    const Vec<int,2>& hi = bx.hi();
    loop_box_2(bx,i,j){
        res(i,j) = data(2*i,2*j);
    }
  }
};

template <>
class Injection<1>: public Restrictor<1>{
public:
  using BaseClass = Restrictor<1>;

  Injection(const RectDomain<1>& afineDomain, const RectDomain<1>&
            acoarseDomain):Restrictor(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,1>& data,
             Tensor<Real,1>& res) const{
    auto& bx = coarseDomain;
    const Vec<int, 1>& lo = bx.lo();
    const Vec<int,1>& hi = bx.hi();
    loop_box_1(bx,i){
        res(i) = data(2*i);
    }
  }
};

template <int Dim> class FullWeighting;

template <>
class FullWeighting<2>: public Restrictor<2>{
public:
  using BaseClass = Restrictor<2>;

  FullWeighting(const RectDomain<2>& afineDomain, const RectDomain<2>&
                acoarseDomain):Restrictor(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,2>& data,
             Tensor<Real,2>& res) const{
    auto& bx = coarseDomain;
    const Vec<int,2>& lo = bx.lo();
    const Vec<int,2>& hi = bx.hi();
    loop_box_2(bx,i,j){
        res(i,j) = (4*data(2*i,2*j)+2*data(2*i+1,2*j)+2*data(2*i,2*j+1)+2*data(2*i-1,2*j)+2*data(2*i,2*j-1)+data(2*i+1,2*j+1)+data(2*i-1,2*j+1)+data(2*i+1,2*j-1)+data(2*i-1,2*j-1))/16.0;
    }
  }
  
};

template <>
class FullWeighting<1>: public Restrictor<1>{
public:
  using BaseClass = Restrictor<1>;

  FullWeighting(const RectDomain<1>& afineDomain, const RectDomain<1>&
                acoarseDomain):Restrictor(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,1>& data,
             Tensor<Real,1>& res) const{
    auto& bx = coarseDomain;
    const Vec<int,1>& lo = bx.lo();
    const Vec<int,1>& hi = bx.hi();
    loop_box_1(bx,i){
        res(i) = (2*data(2*i)+data(2*i-1)+data(2*i+1))/4.0;
    }
  }
  
};


  


#endif // _RESTRICTOR_H_
