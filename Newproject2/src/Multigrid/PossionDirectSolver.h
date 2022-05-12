#ifndef _POSSIONDIRECTSOLVER_H_
#define _POSSIONDIRECTSOLVER_H_

#include "Multigrid/Laplacian.h"
#include "RegularGrid/GhostFiller.h"
#include <cblas.h>
#include <lapacke.h>

template <int Dim> class PossionDirectSolver;

template <int Dim>
class PossionDirectSolver{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;

  PossionDirectSolver() = default;
  
  PossionDirectSolver(const RectDomain<Dim>& adomain);

  void solve(Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
             rhs, const char* const BCTypes, const int
             maxIter = 20) const;
protected:
  RectDomain<Dim> domain;
};

template <int Dim>
PossionDirectSolver<Dim>::PossionDirectSolver(const RectDomain<Dim>& adomain){
  assert(adomain.getCentering() == NodeCentered);
  //assert(adomain.getNumGhost() == 1);
  domain = adomain;
}

template <int Dim>
void PossionDirectSolver<Dim>::solve(Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                                     rhs, const char* const BCTypes, const int maxIter) const{
  Laplacian<Dim> LOP(domain);
  GhostFiller<Dim> GF(domain);
  Tensor<Real,Dim> tmpres(domain.getGhostedBox());
  for (int i = 0 ; i < maxIter ; i+=2){
    GF.fillAllSides(phi,BCTypes);
    LOP.smooth(phi,rhs,tmpres);
    GF.fillAllSides(tmpres,BCTypes);
    LOP.smooth(tmpres,rhs,phi);
  }
  GF.fillAllSides(phi,BCTypes);
}

// template <int Dim>
// void PossionDirectSolver<Dim>::Solve(Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
//                               rhs, const char* BCTypes, const Real
//                               weight, const int maxIter) const {
//   const Vec<Real,Dim>& dx = Domain.spacing();
//   Real h = prod(dx);
//   const Vec<int,Dim>& lo = Domain.lo();
//   const Vec<int,Dim>& hi = Domain.hi();
//   const int numl = 3;
//   Tensor<Real,2> Matrix(iVec{numl*numl,numl*numl});
//   Matrix = 0.0;
//   Tensor<Real,1> RHS(Vec<int,1>{numl*numl});
//   RHS = 0.0;
//   loop_box_2(Domain, i, j){
//     iVec node{i,j};
//     if (i != lo[1] && i != hi[1] && j != lo[0] && j != hi[0]){
//       Matrix(i+j*numl,i+j*numl) = -4/h;
//       Matrix(i+j*numl,i+1+j*numl) = 1/h;
//       Matrix(i+j*numl,i-1+j*numl) = 1/h;
//       Matrix(i+j*numl,i+(j+1)*numl) = 1/h;
//       Matrix(i+j*numl,i+(j-1)*numl) = 1/h;
//       RHS(i+j*numl) = rhs(node);
//     }
//     else{
//       Matrix(i+j*numl,i+j*numl) = 1.0;
//       RHS(i+j*numl) = 0.0;
//     }    
//   }
//   Tensor<int,1> ipiv(numl*numl);
//   auto info = LAPACKE_dgesv(LAPACK_COL_MAJOR, numl*numl, 1, Matrix.data(), numl*numl, ipiv.data(), RHS.data(), numl*numl);
//   if(info != 0)
//     throw std::runtime_error("Solve() - DGESV");
//   loop_box_2(Domain, i, j){
//     phi(i,j)=RHS(j*numl+i);
//   }
//   GhostFiller<Dim> GF(Domain);
//   GF.fillAllSides(phi,BCTypes,0);
// }

#endif // _POSSIONDIRECTSOLVER_H_
