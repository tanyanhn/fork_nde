#ifndef _WEIGHTEDJACOBI_H_
#define _WEIGHTEDJACOBI_H_

#include "Multigrid/Smoother.h"

template <int Dim> class WeightedJacobi;

template <>
class WeightedJacobi<1> : public Smoother<1>{
public:
  WeightedJacobi(const RectDomain<1>& adomain):Smoother(adomain){}
  
  void apply(const Tensor<Real,1>& phi, const Tensor<Real,1>&
             rhs, Tensor<Real,1>& res) const{
    const Vec<Real,1>& dx = domain.spacing();
    const Box<1>& bx = domain;
    assert(res.box().contain(bx));
    const Real b = 2.0/(dx[0]*dx[0]);
    loop_box_1(bx,i){
      Real tmp = (phi(i-1)+phi(i+1))/(dx[0]*dx[0]);
      Real r = (tmp - rhs(i))/b;
      res(i) = r * weight + phi(i) * (1 - weight);
    }
  }
protected:
  const Real weight = 2.0/3;
};

template <>
class WeightedJacobi<2> : public Smoother<2>{
public:
  WeightedJacobi(const RectDomain<2>& adomain):Smoother(adomain){}
  
  void apply(const Tensor<Real,2>& phi, const Tensor<Real,2>&
             rhs, Tensor<Real,2>& res) const{
    const Vec<Real,2>& dx = domain.spacing();
    const Box<2>& bx = domain;
    assert(res.box().contain(bx));
    const Real b = 2.0/(dx[0]*dx[0])+2.0/(dx[1]*dx[1]);
    loop_box_2(bx,i,j){
      Real tmp = (phi(i-1,j)+phi(i+1,j))/(dx[0]*dx[0])+(phi(i,j-1)+phi(i,j+1))/(dx[1]*dx[1]);
      Real r = (tmp - rhs(i,j))/b;
      res(i,j) = r * weight + phi(i,j) * (1 - weight);
    }
  }
protected:
  const Real weight = 4.0/5;
};



#endif // _WEIGHTEDJACOBI_H_
