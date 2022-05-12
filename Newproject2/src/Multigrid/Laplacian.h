#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_

#include "RegularGrid/RectDomain.h"
#include "Multigrid/WeightedJacobi.h"

template <int Dim> class Laplacian;

template <int Dim>
class Laplacian{
public:
  Laplacian(const RectDomain<Dim> adomain);

  void apply(const Tensor<Real,Dim>& phi, Tensor<Real,Dim>&
             res) const;

  void smooth(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
              rhs, Tensor<Real,Dim>& res) const;
  
  void computeResidul(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                      rhs, Tensor<Real,Dim>& res) const;

  Real computeNorm(const Tensor<Real,Dim>& phi, int p = 0) const;
protected:
  RectDomain<Dim> domain;
  Smoother<Dim>* psmoother;
};

template <int Dim>
Laplacian<Dim>::Laplacian(const RectDomain<Dim> adomain){
  assert(adomain.getCentering() == NodeCentered);
  assert(adomain.getNumGhost() == 1);
  domain = adomain;
  psmoother = new WeightedJacobi<Dim>{adomain};
}

template <>
void Laplacian<1>::apply(const Tensor<Real,1>& phi, Tensor<Real,1>&
                         res) const{
  const Box<1>& bx = domain;
  const Vec<Real,1>& dx = domain.spacing();
  assert(res.box().contain(bx));
  loop_box_1(bx,i){
    res(i) = (-2*phi(i)+phi(i-1)+phi(i+1))/(dx[0]*dx[0]);
  }
}

template <>
void Laplacian<2>::apply(const Tensor<Real,2>& phi, Tensor<Real,2>&
                         res) const{
  const Box<2>& bx = domain;
  const Vec<Real,2>& dx = domain.spacing();
  assert(res.box().contain(bx));
  loop_box_2(bx,i,j){
    res(i,j) = (-2*phi(i,j)+phi(i-1,j)+phi(i+1,j))/(dx[0]*dx[0]) +
      (-2*phi(i,j)+phi(i,j-1)+phi(i,j+1))/(dx[1]*dx[1]);
  }
}


template <int Dim>
void Laplacian<Dim>::smooth(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                            rhs, Tensor<Real,Dim>& res) const{
  psmoother->apply(phi,rhs,res);
}



template <int Dim>
void Laplacian<Dim>::computeResidul(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                                    rhs, Tensor<Real,Dim>& res) const{
   this->apply(phi,res);
   res = rhs - res;
}

template <int Dim>
Real Laplacian<Dim>::computeNorm(const Tensor<Real,Dim>& phi, int p) const{
  const Box<Dim>& bx = domain;
  if(p == 0)
    return norm(phi.slice(bx), 0);
  else if(p == 1)
    return norm(phi.slice(bx), 1) * prod(domain.spacing());
  else if(p == 2)
    return norm(phi.slice(bx), 2) * sqrt(prod(domain.spacing()));
  return 0.0;
}






#endif // _LAPLACIAN_H_
