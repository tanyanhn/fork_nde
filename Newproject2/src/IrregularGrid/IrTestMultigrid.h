#ifndef _IRTESTMULTIGRID_H_
#define _IRTESTMULTIGRID_H_

#include "IrregularGrid/IrMGFactory.h"

template <int Dim, class RST, class ITPL> class IrTestMultigrid;

template <class RST, class ITPL>
class IrTestMultigrid<2,RST,ITPL>{
public:

  using iVec = Vec<int,2>;

  using rVec = Vec<Real,2>;
  
  template <class T>
  using Vector = std::vector<T>;

  IrTestMultigrid(const IrMGFactory<2,RST,ITPL>& aMGF);

  template <class TFunc>
  Vector<Real> encrypTestVCycle(const TFunc& f, const int numEncryp, const
                          int p = 0);

  template<class TFunc>
  void TestVCycle(const TFunc &f, const int numEncryp);
  
protected:
  IrMGFactory<2,RST,ITPL> MGF;
};

template <class RST, class ITPL>
IrTestMultigrid<2,RST,ITPL>::IrTestMultigrid(const IrMGFactory<2,RST,ITPL>& aMGF):MGF(aMGF){}

template <class RST, class ITPL>
template<class TFunc>
std::vector<Real> IrTestMultigrid<2,RST,ITPL>::encrypTestVCycle(const
                                                                TFunc &f, const int numEncryp, const int p){
  std::vector<Real> res(2*numEncryp+1);
  const char* BCTypes = MGF.getBCTypes();
  const IrRectDomain<2>& Domain = MGF.getDomain();
  int n = Domain.getDomain().size()[0]-1;
  const Real weight = MGF.getWeight();
  const Vec<int,3> numIter = MGF.getNumIter();
  const Real reltol = MGF.getReltol();
  const int maxIter = MGF.getMaxIter();
  const Func<1>* pfunc = Domain.getplbFunc();
  std::cout << "Grid " << n << "X" << n << ":" << std::endl;
  res[0] = MGF.computeError(MGF.SolveVCycle(f),f,p);
  std::vector<IrRectDomain<2> > vDomain{Domain};
  for(int i = 1 ; i <= numEncryp ; i++){
    RectDomain<2> nDomain = vDomain.back().getDomain().refine();
    IrRectDomain<2> niDomain(nDomain,pfunc,1e-16);
    vDomain.push_back(niDomain);
    IrMGFactory<2,RST,ITPL> nMGF(niDomain,pfunc,BCTypes);
    nMGF.setParam(weight,numIter[0],numIter[1],numIter[2],reltol,maxIter);
    n*=2;
    std::cout << "Grid " << n << "X" << n << ":" << std::endl;
    res[2*i] = nMGF.computeError(nMGF.SolveVCycle(f),f,p);
    res[2*i-1] = log(res[2*i-2]/res[2*i])/log(2);
  }
  return res;
}

template <class RST, class ITPL>
template<class TFunc>
void IrTestMultigrid<2,RST,ITPL>::TestVCycle(const TFunc &f, const int numEncryp){
  std::vector<Real> res1(2*numEncryp+1);
  std::vector<Real> res2(2*numEncryp+1);
  std::vector<Real> res3(2*numEncryp+1);
  const char* BCTypes = MGF.getBCTypes();
  const IrRectDomain<2>& Domain = MGF.getDomain();
  int n = Domain.getDomain().size()[0]-1;
  const Real weight = MGF.getWeight();
  const Vec<int,3> numIter = MGF.getNumIter();
  const Real reltol = MGF.getReltol();
  const int maxIter = MGF.getMaxIter();
  const Func<1>* pfunc = Domain.getplbFunc();
  std::cout << "Grid " << n << "X" << n << ":" << std::endl;
  Tensor<Real,2> phi = MGF.SolveVCycle(f);
  res1[0] = MGF.computeError(phi,f,0);
  res2[0] = MGF.computeError(phi,f,1);
  res3[0] = MGF.computeError(phi,f,2);
  std::vector<IrRectDomain<2> > vDomain{Domain};
  for(int i = 1 ; i <= numEncryp ; i++){
    RectDomain<2> nDomain = vDomain.back().getDomain().refine();
    IrRectDomain<2> niDomain(nDomain,pfunc,1e-16);
    vDomain.push_back(niDomain);
    IrMGFactory<2,RST,ITPL> nMGF(niDomain,pfunc,BCTypes);
    nMGF.setParam(weight,numIter[0],numIter[1],numIter[2],reltol,maxIter);
    n*=2;
    std::cout << "Grid " << n << "X" << n << ":" << std::endl;
    Tensor<Real,2> cphi = nMGF.SolveVCycle(f);
    res1[2*i] = nMGF.computeError(cphi,f,0);
    res2[2*i] = nMGF.computeError(cphi,f,1);
    res3[2*i] = nMGF.computeError(cphi,f,2);
    res1[2*i-1] = log(res1[2*i-2]/res1[2*i])/log(2);
    res2[2*i-1] = log(res2[2*i-2]/res2[2*i])/log(2);
    res3[2*i-1] = log(res3[2*i-2]/res3[2*i])/log(2);
  }
  std::cout << "----------------------------------------" << std::endl;
  const Vec<Real,2>& dx = Domain.getDomain().spacing();
  const Vec<int,2>& lo = Domain.getDomain().lo();
  const Vec<int,2>& hi = Domain.getDomain().hi();
  n = hi[0] - lo[0];
  std::cout << "\\noindent Domain: " << "$[" << lo[0] << "," << (hi[0]-lo[0])*dx[0] \
            << "]\\times[" << lo[1] << "," << (hi[1]-lo[1])*dx[1] <<    \
    "]" << "$\\\\" << std::endl;
  std::cout << "BCtype : " << BCTypes[0] << " , "<< BCTypes[1] << " , " \
            << BCTypes[2] << " , " << BCTypes[3] << std::endl;
  std::cout << "\\begin{table}[htbp]" << std::endl;
  std::cout << "\\centering\\begin{tabular}{c|";
  for (int i = 1 ; i <= 2*numEncryp+1 ; i++)
    std::cout << "c";
  std::cout << "}" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$n$&" << n;
  for (int i = 1 ; i <= numEncryp ; i++){
    n*=2;
    std::cout << "&ratio&" << n;
  }
  std::cout << "\\\\" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$\\|\\mathrm{E}\\|_1$&" << res2[0];
  for (int i = 1 ; i < 2*numEncryp+1 ; i++){
    std::cout << "&" << res2[i];
  }
  std::cout << "\\\\" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$\\|\\mathrm{E}\\|_2$&" << res3[0];
  for (int i = 1 ; i < 2*numEncryp+1 ; i++){
    std::cout << "&" << res3[i];
  }
  std::cout << "\\\\" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$\\|\\mathrm{E}\\|_{\\infty}$&" << res1[0];
  for (int i = 1 ; i < 2*numEncryp+1 ; i++){
    std::cout << "&" << res1[i];
  }
  std::cout << "\\\\" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\caption{}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}




#endif // _TESTMULTIGRID_H_
