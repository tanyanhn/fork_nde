#ifndef _IRMGFACTORY_H_
#define _IRMGFACTORY_H_

#include "IrregularGrid/IrMultigrid.h"
#include <string>
#include <fstream>

template <int Dim, class RST, class ITPL> class IrMGFactory;

template <class RST, class ITPL>
class IrMGFactory<2,RST,ITPL>{
public:

  using iVec = Vec<int,2>;

  using rVec = Vec<Real,2>;

  template <class T>
  using Vector = std::vector<T>;

  IrMGFactory(const IrRectDomain<2>& aDomain, const Func<1>* pfunc, const char* aBCTypes);

  IrMGFactory(const Vec<Real,2> &alo, const Vec<Real,2> &ahi, const
              int an, const Func<1>* pfunc, const char* aBCTypes);
  
  void setParam(const Real aweight, const int anumIter1, const int
                anumIter2, const int anumIter3, const Real areltol,
                const int amaxIter);

  const IrRectDomain<2>& getDomain() const;

  const char* getBCTypes() const;

  const Real getWeight() const;

  const Vec<int,3> getNumIter() const;

  const Real getReltol() const;

  const int getMaxIter() const;
  
  template <class TFunc>
  Tensor<Real,2> SolveVCycle(const TFunc& func);

  template <class TFunc>
  Tensor<Real,2> SolveFMCycle(const TFunc& func);

  template <class TFunc>
  Real computeError(const Tensor<Real,2>& res, const TFunc& func,
                    const int p = 0);

  void plot(const Tensor<Real,2>& res, const std::string &file);
protected:
  IrMultigrid<2> MG; 
};

template <class RST, class ITPL>
IrMGFactory<2,RST,ITPL>::IrMGFactory(const IrRectDomain<2>& aDomain,
                                     const Func<1>* pfunc, const char* aBCTypes){
  Vector<IrRectDomain<2> > vD{aDomain};
  int an = (aDomain.getDomain().size())[0]-1;
  for(int n = an ; n >= 8 ; n/=2){
    RectDomain<2> Domain = vD.back().getDomain().coarsen();
    IrRectDomain<2> iDomain(Domain,pfunc,1e-16);
    vD.push_back(iDomain);
  }
  int size = vD.size();
  Vector<IrRestriction<2>*> vR;
  Vector<IrInterpolation<2>*> vI;
  for (int i = 0 ; i < size - 1 ; i++){
    IrRestriction<2>* Rst = new RST(vD[i],vD[i+1]);
    vR.push_back(Rst);
    IrInterpolation<2>* Itpl = new ITPL(vD[i],vD[i+1]);
    vI.push_back(Itpl);
  }
  MG = IrMultigrid<2>(vD,vR,vI);
  MG.setBCType(aBCTypes);
}

template <class RST, class ITPL>
IrMGFactory<2,RST,ITPL>::IrMGFactory(const Vec<Real,2> &alo, const Vec<Real,2> &ahi, const
                                     int an, const Func<1>* pfunc, const char* aBCTypes){
  Vec<int,2> lo(0);
  Vec<int,2> hi(an);
  Box<2> bx(lo,hi);
  Vec<Real,2> dx = (ahi-alo)/an;
  RectDomain<2> aDomain(bx,dx,NodeCentered,1);
  IrRectDomain<2> aiDomain(aDomain,pfunc,1e-16);
  Vector<IrRectDomain<2> > vD{aiDomain};
  for(int n = an ; n >= 8 ; n/=2){
    RectDomain<2> Domain = vD.back().getDomain().coarsen();
    IrRectDomain<2> iDomain(Domain,pfunc,1e-16);
    vD.push_back(iDomain);
  }
  int size = vD.size();
  Vector<IrRestriction<2>*> vR;
  Vector<IrInterpolation<2>*> vI;
  for (int i = 0 ; i < size - 1 ; i++){
    IrRestriction<2>* Rst = new RST(vD[i],vD[i+1]);
    vR.push_back(Rst);
    IrInterpolation<2>* Itpl = new ITPL(vD[i],vD[i+1]);
    vI.push_back(Itpl);
  }
  MG = IrMultigrid<2>(vD,vR,vI);
  MG.setBCType(aBCTypes);
}

template <class RST, class ITPL>
void IrMGFactory<2,RST,ITPL>::setParam(const Real aweight, const int anumIter1, const int
                                     anumIter2, const int anumIter3, const Real areltol,
                                     const int amaxIter){
  MG.setParam(aweight,anumIter1,anumIter2,anumIter3,areltol,amaxIter);
}

template <class RST, class ITPL>
const IrRectDomain<2>& IrMGFactory<2,RST,ITPL>::getDomain() const{
  return (MG.getvDomain())[0];
}

template <class RST, class ITPL>
const char* IrMGFactory<2,RST,ITPL>::getBCTypes() const{
  return MG.getBCTypes();
}

template <class RST, class ITPL>
const Real IrMGFactory<2,RST,ITPL>::getWeight() const{
  return MG.getWeight();
}

template <class RST, class ITPL>
const Vec<int,3> IrMGFactory<2,RST,ITPL>::getNumIter() const{
  return MG.getNumIter();
}

template <class RST, class ITPL>
const Real IrMGFactory<2,RST,ITPL>::getReltol() const{
  return MG.getReltol();
}

template <class RST, class ITPL>
const int IrMGFactory<2,RST,ITPL>::getMaxIter() const{
  return MG.getMaxIter();
}

template <class RST, class ITPL>
template <class TFunc>
Tensor<Real,2> IrMGFactory<2,RST,ITPL>::SolveVCycle(const TFunc& func){
  const IrRectDomain<2>& Domain = (MG.getvDomain())[0];
  IrFuncFiller<2> FuncF(Domain);
  IrGhostFiller<2> GhostF(Domain);
  Tensor<Real,2> rhs(Domain.getDomain().getGhostedBox());
  FuncF.fillFunc(rhs,func);
  const char* BCTypes = MG.getBCTypes();
  GhostF.fillAllSide(rhs,BCTypes,func);
  //std::cout << rhs << std::endl;
  Tensor<Real,2> phi(Domain.getDomain().getGhostedBox());
  MG.SolveVCycle(phi,rhs,func);
  //std::cout << phi << std::endl;
  return phi;
}

template <class RST, class ITPL>
template <class TFunc>
Tensor<Real,2> IrMGFactory<2,RST,ITPL>::SolveFMCycle(const TFunc& func){
  const IrRectDomain<2>& Domain = (MG.getvDomain())[0];
  IrFuncFiller<2> FuncF(Domain);
  IrGhostFiller<2> GhostF(Domain);
  Tensor<Real,2> rhs(Domain.getDomain().getGhostedBox());
  FuncF.fillFunc(rhs,func);
  const char* BCTypes = MG.getBCTypes();
  GhostF.fillAllSide(rhs,BCTypes,func);
  Tensor<Real,2> phi(Domain.getDomain().getGhostedBox());
  MG.SolveFMCycle(phi,rhs,func);
  return phi;
}

template <class RST, class ITPL>
template <class TFunc>
Real IrMGFactory<2,RST,ITPL>::computeError(const Tensor<Real,2>& res,
                                           const TFunc& func, const int
                                           p){
  const IrRectDomain<2>& iDomain = (MG.getvDomain())[0];
  const RectDomain<2>& Domain = iDomain.getDomain();
  Tensor<Real,2> exactres(Domain);
  const Vec<int,2>& lo = Domain.lo();
  const Vec<Real,2>& dx = Domain.spacing();
  Tensor<Real,2> Rres = res.slice(Domain);
  loop_box_2(Domain,i,j){
    Vec<int,2> Node{i,j};
    if (iDomain.isInside(Node)){
      Vec<Real,2> rNode = (Node-lo)*dx;
      exactres(i,j) = func.F(rNode);
    }
    else
      Rres(i,j) = 0.0;
  }
  // std::cout << exactres << std::endl;
  // std::cout << Rres << std::endl;
  IrPossionOP<2> POP(iDomain);
  return POP.computeNorm(Rres-exactres,p);
}


template <class RST, class ITPL>
void IrMGFactory<2,RST,ITPL>::plot(const Tensor<Real,2>& res, const std::string &file){
  const IrRectDomain<2>& iDomain = (MG.getvDomain())[0];
  const RectDomain<2>& Domain = iDomain.getDomain();
  Tensor<Real,2> Rres = res.slice(Domain);
  const Vec<int,2> &lo = Domain.lo();
  const Vec<Real,2> &dx = Domain.spacing();
  std::ofstream os;
  os.open(file);
  os << "res=[\n" << std::endl;
  loop_box_2(Domain,i,j){
      Vec<int,2> Node{i,j};
      if (iDomain.isInside(Node)){
        Vec<Real,2> rNode = (Node - lo)*dx;
        os << rNode[0] << "," << rNode[1] << "," << Rres(i,j) << ";\n"  \
           << std::endl;
      }
  }
  os << "];\n" << std::endl;
  os << "x=res(:,1);y=res(:,2);z=res(:,3);\n" << std::endl;
  os << "T = delaunay(x,y);trimesh(T,x,y,z);\n" << std::endl;
  os.close();
}


#endif // _IRMGFACTORY_H_
