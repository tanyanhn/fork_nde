#ifndef _IRBOTTOMSOLVER_H_
#define _IRBOTTOMSOLVER_H_

#include "IrregularGrid/IrPossionOP.h"
#include "IrregularGrid/IrGhostFiller.h"
#include <cblas.h>
#include <lapacke.h>

template <int Dim> class IrBottomSolver;

template <>
class IrBottomSolver<2>{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;

  IrBottomSolver() = default;
  
  IrBottomSolver(const IrRectDomain<2>& aDomain);

  void Solve(Tensor<Real,2>& phi, const Tensor<Real,2>&
             rhs, const char* BCTypes, const Real weight, const int maxIter = 20);
protected:
  IrRectDomain<2> Domain;
};


IrBottomSolver<2>::IrBottomSolver(const IrRectDomain<2>& aDomain):Domain(aDomain){}


// void IrBottomSolver<2>::Solve(Tensor<Real,2>& phi, const Tensor<Real,2>&
//                               rhs, const char* BCTypes, const Real weight, const int maxIter){
//   IrPossionOP<2> POP(Domain);
//   IrGhostFiller<2> GF(Domain);
//   Tensor<Real,2> tmpres(Domain.getDomain().getGhostedBox());
//   for (int i = 0 ; i < maxIter ; i+=2){
//     GF.fillAllSide(phi,BCTypes,0);
//     POP.relaxLaplacian(phi,rhs,tmpres,1.0);
//     GF.fillAllSide(tmpres,BCTypes,0);
//     POP.relaxLaplacian(tmpres,rhs,phi,1.0);
//   }
//   GF.fillAllSide(phi,BCTypes,0);
// }

void IrBottomSolver<2>::Solve(Tensor<Real,2>& phi, const Tensor<Real,2>&
                              rhs, const char* BCTypes, const Real
                              weight, const int maxIter){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const rVec& dx = rDomain.spacing();
  Real h = dx[0]*dx[1];
  const iVec& lo = rDomain.lo();
  const iVec& hi = rDomain.hi();
  const int numl = 5;
  Tensor<Real,2> Matrix(iVec{numl*numl,numl*numl});
  Matrix = 0.0;
  Tensor<Real,1> RHS(Vec<int,1>{numl*numl});
  RHS = 0.0;
  loop_box_2(rDomain, i, j){
    iVec node{i,j};
    if (i != lo[1] && i != hi[1] && j != lo[0] && j != hi[0]){
      Matrix(i+j*numl,i+j*numl) = -4/h;
      Matrix(i+j*numl,i+1+j*numl) = 1/h;
      Matrix(i+j*numl,i-1+j*numl) = 1/h;
      Matrix(i+j*numl,i+(j+1)*numl) = 1/h;
      Matrix(i+j*numl,i+(j-1)*numl) = 1/h;
      RHS(i+j*numl) = rhs(node);
    }
    else{
      Matrix(i+j*numl,i+j*numl) = 1.0;
      RHS(i+j*numl) = 0.0;
    }    
  }
  Tensor<int,1> ipiv(numl*numl);
  auto info = LAPACKE_dgesv(LAPACK_COL_MAJOR, numl*numl, 1, Matrix.data(), numl*numl, ipiv.data(), RHS.data(), numl*numl);
  if(info != 0)
    throw std::runtime_error("Solve() - DGESV");
  loop_box_2(rDomain, i, j){
    phi(i,j)=RHS(j*numl+i);
  }
  IrGhostFiller<2> GF(Domain);
  GF.fillAllSide(phi,BCTypes,0);
}

#endif // _IRBOTTOMSOLVER_H_
