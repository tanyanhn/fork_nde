#ifndef _IRPOSSIONOP_H_
#define _IRPOSSIONOP_H_

#include "IrregularGrid/IrRectDomain.h"
#include "Core/Tensor.h"
#include "Core/TensorSlice.h"

#include <vector>

template <int Dim> class IrPossionOP;

template <>
class IrPossionOP<2>{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;

  template <class T>
  using Vector = std::vector<T>;
  
  IrPossionOP(const IrRectDomain<2>& aDomain);

  void Laplacian(const Tensor<Real,2>& phi, Tensor<Real,2>& res);

  void relaxLaplacian(const Tensor<Real,2>& phi, const Tensor<Real,2>&
                      rhs, Tensor<Real,2>& res, Real weight, const
                      Func<2>& func);

  void relaxLaplacian(const Tensor<Real,2>& phi, const Tensor<Real,2>&
                      rhs, Tensor<Real,2>& res, Real weight, int);

  void computeResidul(const Tensor<Real,2>& phi, const Tensor<Real,2>&
                      rhs, Tensor<Real,2>& res);
  
  Real computeNorm(const Tensor<Real,2>& phi, int p = 0);
protected:
  IrRectDomain<2> Domain;
};


IrPossionOP<2>::IrPossionOP(const IrRectDomain<2>& aDomain):Domain(aDomain){}

void IrPossionOP<2>::Laplacian(const Tensor<Real,2>& phi,
                               Tensor<Real,2>& res){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const rVec& dx = rDomain.spacing();
  assert(res.box().contain(rDomain));
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (Domain.isInside(Node)){
      res(Node) = (-2*phi(i,j)+phi(i-1,j)+phi(i+1,j))/(dx[0]*dx[0]) +
        (-2*phi(i,j)+phi(i,j-1)+phi(i,j+1))/(dx[1]*dx[1]);
      }
  }
}


// void IrPossionOP<2>::relaxLaplacian(const Tensor<Real,2>& phi, const Tensor<Real,2>&
//                                     rhs, Tensor<Real,2>& res, Real
//                                     weight, const Func<2>& func){
//   const RectDomain<2>& rDomain = Domain.getDomain();
//   const rVec& dx = rDomain.spacing();
//   assert(res.box().contain(rDomain));
//   const Real b = 2.0/(dx[0]*dx[0])+2.0/(dx[1]*dx[1]);
//   loop_box_2(rDomain,i,j){
//     iVec Node{i,j};
//     if (Domain.isInside(Node)){
//       Real tmp = (phi(i-1,j)+phi(i+1,j))/(dx[0]*dx[0])+(phi(i,j-1)+phi(i,j+1))/(dx[1]*dx[1]);
//       Real r = (tmp - rhs(i,j))/b;
//       res(i,j) = r * weight + phi(i,j) * (1 - weight);
//       }
    
//   }
// }


void IrPossionOP<2>::relaxLaplacian(const Tensor<Real,2>& phi, const Tensor<Real,2>&
                                    rhs, Tensor<Real,2>& res, Real
                                    weight, const Func<2>& func){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const iVec& lo = rDomain.lo();
  const rVec& dx = rDomain.spacing();
  assert(res.box().contain(rDomain));
  //Tensor<Real,2> tmpphi = phi;
  const Real b = 2.0/(dx[0]*dx[0])+2.0/(dx[1]*dx[1]);
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (Domain.isInside(Node)){
      if (Domain.isRePoint(Node)){
        Real tmp = (phi(i-1,j)+phi(i+1,j))/(dx[0]*dx[0])+(phi(i,j-1)+phi(i,j+1))/(dx[1]*dx[1]);
        Real r = (tmp - rhs(i,j))/b;
        res(i,j) = r*weight + phi(i,j)*(1-weight);
      }
      else{
        Vec<Real,2> frac = Domain.getProjPoint(Node);
        rVec projpt = (Node+iVec{1,1}*frac-lo)*dx;
        Real projvalue;
        iVec Node1 = Node + iVec{0,-1};
        Real tmp1 = (phi(Node)-phi(Node1))*frac[1]+phi(Node);
        if (frac[0] >= 0){
          iVec Node2 = Node + iVec{1,0};
          iVec Node3 = Node + iVec{1,-1};
          Real tmp2 = (phi(Node2)-phi(Node3))*frac[1]+phi(Node2);
          projvalue = (tmp2-tmp1)*frac[0]+tmp1;
        }
        else{
          iVec Node2 = Node + iVec{-1,0};
          iVec Node3 = Node + iVec{-1,-1};
          Real tmp2 = (phi(Node2)-phi(Node3))*frac[1]+phi(Node2);
          projvalue = (tmp1-tmp2)*frac[0]+tmp1;
        }
        Real rvalue = func.F(projpt);
        res(i,j) = phi(i,j)+2*(rvalue - projvalue);
      }
    } 
  }
}

void IrPossionOP<2>::relaxLaplacian(const Tensor<Real,2>& phi, const Tensor<Real,2>&
                                    rhs, Tensor<Real,2>& res, Real weight, int){
  D2Func0 f;
  relaxLaplacian(phi, rhs, res, weight, f);
}

void IrPossionOP<2>::computeResidul(const Tensor<Real,2>& phi, const Tensor<Real,2>&
                                    rhs, Tensor<Real,2>& res){
  Laplacian(phi,res);
  res = rhs - res;
}


Real IrPossionOP<2>::computeNorm(const Tensor<Real,2>& phi, int p){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const Box<2>& bx = rDomain;
  if(p == 0)
    return norm(phi.slice(bx), 0);
  else if(p == 1)
    return norm(phi.slice(bx), 1) * prod(rDomain.spacing());
  else if(p == 2)
    return norm(phi.slice(bx), 2) * sqrt(prod(rDomain.spacing()));
  return 0.0;
}


#endif //_IRPOSSIONOP_H_
