#ifndef _IRMULTIGRID_H_
#define _IRMULTIGRID_H_

#include "IrregularGrid/IrPossionOP.h"
#include "IrregularGrid/IrBottomSolver.h"
#include "IrregularGrid/IrRestriction.h"
#include "IrregularGrid/IrInterpolation.h"
#include "IrregularGrid/IrFuncFiller.h"
#include <vector>
#include <memory>
#include <cstring>

template <int Dim> class IrMultigrid;

template <>
class IrMultigrid<2>{
public:
  template <class T>
  using Vector = std::vector<T>;

  using VPR = Vector<IrRestriction<2>*>;

  using VPI = Vector<IrInterpolation<2>*>;

  IrMultigrid() = default;

  IrMultigrid(const Vector<IrRectDomain<2> >& avDomain, const VPR&
              avpRestriction, const VPI& avpInterpolation);

  ~IrMultigrid();

  void setBCType(const char* aBCTypes);

  void setParam(const Real aweight, const int anumIter1, const int
                anumIter2, const int anumIter3, const Real areltol,
                const int amaxIter);

  const Vector<IrRectDomain<2> >& getvDomain() const;

  const char* getBCTypes() const;

  const Real getWeight() const;

  const Vec<int,3> getNumIter() const;

  const Real getReltol() const;

  const int getMaxIter() const;

  template <class TFunc>
  Real SolveVCycle(Tensor<Real,2>& phi, const Tensor<Real,2>& rhs,
  const TFunc& f);

  template <class TFunc>
  Real SolveFMCycle(Tensor<Real,2>& phi, const Tensor<Real,2>& rhs,
  const TFunc& f);

  protected:
  void VCycle(int depth, Tensor<Real,2>& phi, const Tensor<Real,2>&
              rhs);

  void FMCycle(int depth, Tensor<Real,2>& phi, const Tensor<Real,2>&
               rhs);
  
protected:
  // Operator listed here
  Vector<IrRectDomain<2> > vDomain;
  Vector<IrGhostFiller<2> > vGhostFiller;
  Vector<IrPossionOP<2> > vPossionOP;
  VPR vpRestriction;
  VPI vpInterpolation;
  IrBottomSolver<2> bottomSolver;
  // Boundry condition
  char BCTypes[4];
  // Other parameter
  Real weight;
  int numIter1;
  int numIter2;
  int numIter3;
  Real reltol;
  int maxIter;
};


IrMultigrid<2>::IrMultigrid(const Vector<IrRectDomain<2> >& avDomain, const VPR&
                            avpRestriction, const VPI&
                            avpInterpolation):vDomain(avDomain),
                                              vpRestriction(avpRestriction),vpInterpolation(avpInterpolation){
  const int size = avDomain.size();
  for (int i = 0 ; i < size ; i++){
    vGhostFiller.push_back(IrGhostFiller<2>{avDomain[i]});
    vPossionOP.push_back(IrPossionOP<2>{avDomain[i]});
  }  
  bottomSolver = IrBottomSolver<2>(avDomain[size-1]);
}


IrMultigrid<2>::~IrMultigrid(){}

void IrMultigrid<2>::setBCType(const char* aBCTypes){
  strcpy(BCTypes,aBCTypes);
}


void IrMultigrid<2>::setParam(const Real aweight, const int anumIter1,
                              const int anumIter2, const int anumIter3, const Real areltol,
                              const int amaxIter){
  weight = aweight;
  numIter1 = anumIter1;
  numIter2 = anumIter2;
  numIter3 = anumIter3;
  reltol = areltol;
  maxIter = amaxIter;
}


const std::vector<IrRectDomain<2> >& IrMultigrid<2>::getvDomain() const{
  return vDomain;
}

const char* IrMultigrid<2>::getBCTypes() const{
  return BCTypes;
}

const Real IrMultigrid<2>::getWeight() const{
  return weight;
}

const Real IrMultigrid<2>::getReltol() const{
  return reltol;
}

const Vec<int,3> IrMultigrid<2>::getNumIter() const{
  return Vec<int,3>{numIter1,numIter2,numIter3};
}

const int IrMultigrid<2>::getMaxIter() const{
  return maxIter;
}

void IrMultigrid<2>::VCycle(int depth, Tensor<Real,2>& phi, const Tensor<Real,2>& rhs){
  if (depth == (int)(vDomain.size()) - 1){
    bottomSolver.Solve(phi,rhs,BCTypes,weight,numIter3);
    return;
  }
  //std::cout << phi << std::endl;
  Tensor<Real,2> tmp(vDomain[depth].getDomain().getGhostedBox());
  for (int i = 0 ; i < numIter1 ; i+=2){
    vGhostFiller[depth].fillAllSide(phi,BCTypes,0);
    vPossionOP[depth].relaxLaplacian(phi,rhs,tmp,weight,0);
    vGhostFiller[depth].fillAllSide(tmp,BCTypes,0);
    vPossionOP[depth].relaxLaplacian(tmp,rhs,phi,weight,0);
  }
  vGhostFiller[depth].fillAllSide(phi,BCTypes,0);
  //std::cout << phi << std::endl;
  Tensor<Real,2> rsd(vDomain[depth].getDomain().getGhostedBox());
  vPossionOP[depth].computeResidul(phi,rhs,rsd);
  vGhostFiller[depth].fillAllSide(rsd,BCTypes,0);
  Tensor<Real,2> crsd(vDomain[depth+1].getDomain().getGhostedBox());
  vpRestriction[depth]->applyRestrict(rsd,crsd);
  Tensor<Real,2> cphi(vDomain[depth+1].getDomain().getGhostedBox());
  cphi = 0.0;
  VCycle(depth+1,cphi,crsd);
  vpInterpolation[depth]->applyInterpolate(cphi,phi);
  for (int i = 0 ; i < numIter2 ; i+=2){
    vGhostFiller[depth].fillAllSide(phi,BCTypes,0);
    vPossionOP[depth].relaxLaplacian(phi,rhs,tmp,weight,0);
    vGhostFiller[depth].fillAllSide(tmp,BCTypes,0);
    vPossionOP[depth].relaxLaplacian(tmp,rhs,phi,weight,0);
  }
  vGhostFiller[depth].fillAllSide(phi,BCTypes,0);
}


template <class TFunc>
Real IrMultigrid<2>::SolveVCycle(Tensor<Real,2>& phi, const Tensor<Real,2>& rhs,
                                 const TFunc& f){
  const IrRectDomain<2>& iDomain = vDomain[0];
  Tensor<Real,2> rsd(iDomain.getDomain().getGhostedBox());
  vGhostFiller[0].fillAllSide(phi,BCTypes,f);
  //std::cout << phi << std::endl;
  vPossionOP[0].computeResidul(phi,rhs,rsd);
  vGhostFiller[0].fillAllSide(rsd,BCTypes,0);
  //std::cout << rsd << std::endl;
  Real iRsd = vPossionOP[0].computeNorm(rsd);
  std::cout << iRsd << std::endl;
  Real Rsd[2];
  Rsd[0] = Rsd[1] = iRsd;
  int count = 1;
  Real relRsd,RsdRatio;
  Tensor<Real,2> res(iDomain.getDomain().getGhostedBox());
  while(1){
    Rsd[0] = Rsd[1];
    res = 0.0;
    VCycle(0,res,rsd);
    //std::cout << res << std::endl;
    phi = phi + res;
    vGhostFiller[0].fillAllSide(phi,BCTypes,f);
    vPossionOP[0].computeResidul(phi,rhs,rsd);
    vGhostFiller[0].fillAllSide(rsd,BCTypes,0);
    Rsd[1] = vPossionOP[0].computeNorm(rsd);
    std::cout << Rsd[1] << std::endl;
    relRsd = Rsd[1]/iRsd;
    RsdRatio = Rsd[1]/Rsd[0];
    std::cout << "Iter " << count << " , relrsd = " << relRsd <<\
      " , rsdratio = " << RsdRatio << std::endl;
    if (relRsd < reltol || count > maxIter)
      break;
    count++;
  }
  return relRsd;
}

void IrMultigrid<2>::FMCycle(int depth, Tensor<Real,2>& phi, const Tensor<Real,2>&
                             rhs){
  if (depth != (int)(vDomain.size()) - 1){
    Tensor<Real,2> crhs(vDomain[depth+1].getDomain().getGhostedBox());
    vpRestriction[depth]->applyRestrict(rhs,crhs);
    Tensor<Real,2> cphi(vDomain[depth+1].getDomain().getGhostedBox());
    cphi = 0.0;
    FMCycle(depth+1,cphi,crhs);
    vpInterpolation[depth]->applyInterpolate(cphi,phi);
  }
  VCycle(depth,phi,rhs);
}


template <class TFunc>
Real IrMultigrid<2>::SolveFMCycle(Tensor<Real,2>& phi, const Tensor<Real,2>& rhs,
                                  const TFunc& f){
   const IrRectDomain<2>& iDomain = vDomain[0];
   Tensor<Real,2> rsd(iDomain.getDomain().getGhostedBox());
   vGhostFiller[0].fillAllSide(phi,BCTypes,f);
   vPossionOP[0].computeResidul(phi,rhs,rsd);
   vGhostFiller[0].fillAllSide(rsd,BCTypes,0);
   Real iRsd = vPossionOP[0].computeNorm(rsd);
   Real Rsd[2];
   Rsd[0] = Rsd[1] = iRsd;
   int count = 1;
   Real relRsd,RsdRatio;
   Tensor<Real,2> res(iDomain.getDomain().getGhostedBox());
   while(1){
     Rsd[0] = Rsd[1];
     res = 0.0;
     FMCycle(0,res,rsd);
     phi = phi + res;
     vGhostFiller[0].fillAllSide(phi,BCTypes,f);
     vPossionOP[0].computeResidul(phi,rhs,rsd);
     vGhostFiller[0].fillAllSide(rsd,BCTypes,0);
     Rsd[1] = vPossionOP[0].computeNorm(rsd);
     relRsd = Rsd[1]/iRsd;
     RsdRatio = Rsd[1]/Rsd[0];
     std::cout << "Iter " << count << " , relrsd = " << relRsd <<       \
       " , rsdratio = " << RsdRatio << std::endl;
     if (relRsd < reltol || count > maxIter)
       break;
     count++;
   }
   return relRsd;
}
#endif // _IRMULTIGRID_H_
