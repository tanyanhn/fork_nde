#ifndef _MULTIGRIDSOLVER_H_
#define _MULTIGRIDSOLVER_H_

#include "Multigrid/Laplacian.h"
#include "Multigrid/PossionDirectSolver.h"
#include "Multigrid/Restrictor.h"
#include "Multigrid/Interpolator.h"
#include "Multigrid/MGFactory.h"
#include "RegularGrid/FuncFiller.h"
#include <vector>
#include <memory>
#include <cstring>

template <int Dim> class MultigridSolver;

struct MGParam{
  int numPreIter;
  int numPostIter;
  int numBottomIter;
  Real reltol;
  int maxIter;
};

template <int Dim>
class MultigridSolver{
public:

  template <class T>
  using Vector = std::vector<T>;

  using VPR = Vector<Restrictor<Dim>*>;

  using VPI = Vector<Interpolator<Dim>*>;

  MultigridSolver() = default;

  MultigridSolver(const Vector<RectDomain<Dim> >& avDomain, const VPR&
            avpRestrictor, const VPI& avpInterpolator);

  ~MultigridSolver();

  void setBCType(const char* aBCTypes);

  void setParam(const MGParam& aParam);

  const Vector<RectDomain<Dim> >& getvDomain() const;

  const char* getBCTypes() const;

  const MGParam getMGParam() const;

  const Vec<int,3> getNumIter() const;

  const Real getReltol() const;

  const int getMaxIter() const;
  
  Real solve(Tensor<Real,Dim>& phi, const Tensor<Real,Dim>& rhs,
             const std::array<ScalarFunction<Dim>*,3> pfuncs,
             bool useFMVCycle = 0) const;
protected:
  
  void VCycle(int depth, Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
              rhs) const;

  void FMVCycle(int depth, Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                            rhs) const;
  
protected:
  // Operator listed here
  Vector<RectDomain<Dim> > vDomain;
  Vector<GhostFiller<Dim> > vGhostFiller;
  Vector<Laplacian<Dim> > vLaplacian;
  VPR vpRestrictor;
  VPI vpInterpolator;
  PossionDirectSolver<Dim> bottomSolver;
  // Boundry condition
  char BCTypes[4];
  // Other parameter
  MGParam Param;
};

template <int Dim>
MultigridSolver<Dim>::MultigridSolver(const Vector<RectDomain<Dim> >& avDomain, const VPR&
                                      avpRestrictor, const VPI& avpInterpolator){
  vDomain = avDomain;
  vpRestrictor = avpRestrictor;
  vpInterpolator = avpInterpolator;
  const int size = avDomain.size();
  for (int i = 0 ; i < size ; i++){
    vGhostFiller.push_back(GhostFiller<Dim>{avDomain[i]});
    vLaplacian.push_back(Laplacian<Dim>{avDomain[i]});
  }  
  bottomSolver = PossionDirectSolver<Dim>(avDomain[size-1]);
}

template <int Dim>
MultigridSolver<Dim>::~MultigridSolver(){}

template <int Dim>
void MultigridSolver<Dim>::setBCType(const char* aBCTypes){
  strcpy(BCTypes,aBCTypes);
}

template <int Dim>
void MultigridSolver<Dim>::setParam(const MGParam& aParam){
  Param = aParam;
}

template <int Dim>
const std::vector<RectDomain<Dim> >& MultigridSolver<Dim>::getvDomain() const{
  return vDomain;
}

template <int Dim>
const char* MultigridSolver<Dim>::getBCTypes() const{
  return BCTypes;
}

template <int Dim>
const MGParam MultigridSolver<Dim>::getMGParam() const{
  return Param;
}

template <int Dim>
const Real MultigridSolver<Dim>::getReltol() const{
  return Param.reltol;
}

template <int Dim>
const Vec<int,3> MultigridSolver<Dim>::getNumIter() const{
  return Vec<int,3>{Param.numPreIter,Param.numPostIter,Param.numBottomIter};
}

template <int Dim>
const int MultigridSolver<Dim>::getMaxIter() const{
  return Param.maxIter;
}

template <int Dim>
void MultigridSolver<Dim>::VCycle(int depth, Tensor<Real,Dim>& phi,
                                  const Tensor<Real,Dim>& rhs) const{
  if (depth == (int)(vDomain.size()) - 1){
    bottomSolver.solve(phi,rhs,BCTypes,Param.numBottomIter);
    return;
  }
  Tensor<Real,Dim> tmp(vDomain[depth].getGhostedBox());
  for (int i = 0 ; i < Param.numPreIter ; i+=2){
    vGhostFiller[depth].fillAllSides(phi,BCTypes);
    vLaplacian[depth].smooth(phi,rhs,tmp);
    vGhostFiller[depth].fillAllSides(tmp,BCTypes);
    vLaplacian[depth].smooth(tmp,rhs,phi);
  }
  vGhostFiller[depth].fillAllSides(phi,BCTypes);
  Tensor<Real,Dim> rsd(vDomain[depth].getGhostedBox());
  vLaplacian[depth].computeResidul(phi,rhs,rsd);
  vGhostFiller[depth].fillAllSides(rsd,BCTypes);
  Tensor<Real,Dim> crsd(vDomain[depth+1].getGhostedBox());
  vpRestrictor[depth]->apply(rsd,crsd);
  Tensor<Real,Dim> cphi(vDomain[depth+1].getGhostedBox());
  cphi = 0.0;
  VCycle(depth+1,cphi,crsd);
  vpInterpolator[depth]->apply(cphi,phi);
  for (int i = 0 ; i < Param.numPostIter ; i+=2){
    vGhostFiller[depth].fillAllSides(phi,BCTypes);
    vLaplacian[depth].smooth(phi,rhs,tmp);
    vGhostFiller[depth].fillAllSides(tmp,BCTypes);
    vLaplacian[depth].smooth(tmp,rhs,phi);
  }
  vGhostFiller[depth].fillAllSides(phi,BCTypes);
}

template <int Dim>
void MultigridSolver<Dim>::FMVCycle(int depth, Tensor<Real,Dim>& phi, const Tensor<Real,Dim>&
                                    rhs) const{
  if (depth != (int)(vDomain.size()) - 1){
    Tensor<Real,Dim> crhs(vDomain[depth+1].getGhostedBox());
    vpRestrictor[depth]->apply(rhs,crhs);
    Tensor<Real,Dim> cphi(vDomain[depth+1].getGhostedBox());
    cphi = 0.0;
    FMVCycle(depth+1,cphi,crhs);
    vpInterpolator[depth]->apply(cphi,phi);
  }
  VCycle(depth,phi,rhs);
}


template <int Dim>
Real MultigridSolver<Dim>::solve(Tensor<Real,Dim>& phi, const Tensor<Real,Dim>& rhs,
                                 const
                                 std::array<ScalarFunction<Dim>*,3>
                                 pfuncs, bool useFMVCycle) const{
  const RectDomain<Dim>& Domain = vDomain[0];
  Tensor<Real,Dim> rsd(Domain.getGhostedBox());
  vGhostFiller[0].fillAllSides(phi,BCTypes,pfuncs);
  vLaplacian[0].computeResidul(phi,rhs,rsd);
  vGhostFiller[0].fillAllSides(rsd,BCTypes);
  Real iRsd = vLaplacian[0].computeNorm(rsd);
  Real Rsd[2];
  Rsd[0] = Rsd[1] = iRsd;
  int count = 1;
  Real relRsd,RsdRatio;
  Tensor<Real,Dim> res(Domain.getGhostedBox());
  while(1){
    Rsd[0] = Rsd[1];
    res = 0.0;
    if (useFMVCycle == 0)
      VCycle(0,res,rsd);
    else
      FMVCycle(0,res,rsd);
    phi = phi + res;
    vGhostFiller[0].fillAllSides(phi,BCTypes,pfuncs);
    vLaplacian[0].computeResidul(phi,rhs,rsd);
    vGhostFiller[0].fillAllSides(rsd,BCTypes);
    Rsd[1] = vLaplacian[0].computeNorm(rsd);
    relRsd = Rsd[1]/iRsd;
    RsdRatio = Rsd[1]/Rsd[0];
    std::cout << "Iter " << count << " , relrsd = " << relRsd <<\
      " , rsdratio = " << RsdRatio << std::endl;
    if (relRsd < Param.reltol || count > Param.maxIter)
      break;
    count++;
  }
  return relRsd;
}

template <int Dim>
std::unique_ptr<MultigridSolver<Dim> > createMultigridSolver1(const RectDomain<Dim>
                                                              aDomain, const char* aBCTypes){
  std::vector<RectDomain<Dim> > vD{aDomain};
  int an = (aDomain.size())[0]-1;
  for(int n = an ; n >= 4 ; n/=2){
    RectDomain<Dim> Domain = vD.back().coarsen();
    vD.push_back(Domain);
  }
  int size = vD.size();
  std::vector<Restrictor<Dim>*> vR;
  std::vector<Interpolator<Dim>*> vI;
  for (int i = 0 ; i < size - 1 ; i++){
    Restrictor<Dim>* Rst = new Injection<Dim>(vD[i],vD[i+1]);
    vR.push_back(Rst);
    Interpolator<Dim>* Itpl = new LinearInterpolator<Dim>(vD[i],vD[i+1]);
    vI.push_back(Itpl);
  }
  std::unique_ptr<MultigridSolver<Dim> > result = std::make_unique<MultigridSolver<Dim> >(vD,vR,vI);
  result->setBCType(aBCTypes);
  return result;
}

template <int Dim>
std::unique_ptr<MultigridSolver<Dim> > createMultigridSolver2(const RectDomain<Dim>
                                                              aDomain, const
                                                              char* aBCTypes){
  std::vector<RectDomain<Dim> > vD{aDomain};
  int an = (aDomain.size())[0]-1;
  for(int n = an ; n >= 4 ; n/=2){
    RectDomain<Dim> Domain = vD.back().coarsen();
    vD.push_back(Domain);
  }
  int size = vD.size();
  std::vector<Restrictor<Dim>*> vR;
  std::vector<Interpolator<Dim>*> vI;
  for (int i = 0 ; i < size - 1 ; i++){
    Restrictor<Dim>* Rst = new Injection<Dim>(vD[i],vD[i+1]);
    vR.push_back(Rst);
    Interpolator<Dim>* Itpl = new QuadraticInterpolator<Dim>(vD[i],vD[i+1]);
    vI.push_back(Itpl);
  }
  std::unique_ptr<MultigridSolver<Dim> > result = std::make_unique<MultigridSolver<Dim> >(vD,vR,vI);
  result->setBCType(aBCTypes);
  return result;
}

template <int Dim>
std::unique_ptr<MultigridSolver<Dim> > createMultigridSolver3(const RectDomain<Dim>
                                                              aDomain, const char* aBCTypes){
  std::vector<RectDomain<Dim> > vD{aDomain};
  int an = (aDomain.size())[0]-1;
  for(int n = an ; n >= 4 ; n/=2){
    RectDomain<Dim> Domain = vD.back().coarsen();
    vD.push_back(Domain);
  }
  int size = vD.size();
  std::vector<Restrictor<Dim>*> vR;
  std::vector<Interpolator<Dim>*> vI;
  for (int i = 0 ; i < size - 1 ; i++){
    Restrictor<Dim>* Rst = new FullWeighting<Dim>(vD[i],vD[i+1]);
    vR.push_back(Rst);
    Interpolator<Dim>* Itpl = new LinearInterpolator<Dim>(vD[i],vD[i+1]);
    vI.push_back(Itpl);
  }
  std::unique_ptr<MultigridSolver<Dim> > result = std::make_unique<MultigridSolver<Dim> >(vD,vR,vI);
  result->setBCType(aBCTypes);
  return result;
}

template <int Dim>
std::unique_ptr<MultigridSolver<Dim> > createMultigridSolver4(const RectDomain<Dim>
                                                              aDomain, const char* aBCTypes){
  std::vector<RectDomain<Dim> > vD{aDomain};
  int an = (aDomain.size())[0]-1;
  for(int n = an ; n >= 4 ; n/=2){
    RectDomain<Dim> Domain = vD.back().coarsen();
    vD.push_back(Domain);
  }
  int size = vD.size();
  std::vector<Restrictor<Dim>*> vR;
  std::vector<Interpolator<Dim>*> vI;
  for (int i = 0 ; i < size - 1 ; i++){
    Restrictor<Dim>* Rst = new FullWeighting<Dim>(vD[i],vD[i+1]);
    vR.push_back(Rst);
    Interpolator<Dim>* Itpl = new QuadraticInterpolator<Dim>(vD[i],vD[i+1]);
    vI.push_back(Itpl);
  }
  std::unique_ptr<MultigridSolver<Dim> > result = std::make_unique<MultigridSolver<Dim> >(vD,vR,vI);
  result->setBCType(aBCTypes);
  return result;
}

#endif // _MULTIGRID_H_
