#include "IrregularGrid/IrTestMultigrid.h"
#include <vector>
#include <cmath>

using namespace std;

using iVec = Vec<int,2>;
using rVec = Vec<Real,2>;
template <class T>
using Vector = std::vector<T>;

int main(int argc, char* argv[]){
  
  Vec<Real,2> lo(0);
  Vec<Real,2> hi(1);
  int n = 4;
  char BCTypes[] = {'D','D','D','D'};
  Func<1>* pBCFunc = new BCFunc;
  IrMGFactory<2,FullWeightingRestriction<2>,LinearInterpolation<2>> MGF(lo,hi,n,pBCFunc,BCTypes);
  const Real weight = 2.0/3;
  const int numIter1 = 10;
  const int numIter2 = 10;
  const int numIter3 = 10;
  const Real reltol = 1e-10;
  const int maxIter = 20;
  MGF.setParam(weight,numIter1,numIter2,numIter3,reltol,maxIter);

  D2Func2 f2;
  // Tensor<Real,2> phi = MGF.SolveVCycle(f1);
  // cout << MGF.computeError(phi,f1) << endl;
  IrTestMultigrid<2,FullWeightingRestriction<2>,LinearInterpolation<2>> TM(MGF);
  TM.TestVCycle(f2,3);
  


  
  // Func<1>* pBCFunc = new BCFunc;
  // iVec lo{0,0};
  // iVec hi{16,16};
  // Box<2> bx{lo,hi};
  // rVec dx{1.0/16,1.0/16};
  // RectDomain<2> Domain(bx,dx,NodeCentered,1);
  // IrRectDomain<2> iDomain(Domain,pBCFunc,1e-16);
  
  // cout << iDomain.getProjPoint(iVec{2,1}) << endl;
  
  // D2Func1 Func1;
  // IrFuncFiller<2> irFF(iDomain);
  // Tensor<Real,2> f(iDomain.getDomain().getGhostedBox());
  // irFF.fillFunc(f,Func1);
  // IrGhostFiller<2> irGF(iDomain);
  // char BCTypes[] = {'D','D','D','D'};
  // irGF.fillAllSide(f,BCTypes,Func1);

  // RectDomain<2> cDomain = Domain.coarsen();
  // IrRectDomain<2> ciDomain(cDomain,pBCFunc,1e-16);
  // Tensor<Real,2> cf(ciDomain.getDomain().getGhostedBox());
  // Injection<2>* IRst = new Injection<2>(iDomain,ciDomain);
  // IRst->applyRestrict(f,cf);
  // Tensor<Real,2> fcf(iDomain.getDomain().getGhostedBox());
  // LinearInterpolation<2>* IIpl = new
  // LinearInterpolation<2>(iDomain,ciDomain);
  // IIpl->applyInterpolate(cf,fcf);
  
  // IrPossionOP<2> irPP(iDomain);
  // Tensor<Real,2> res(iDomain.getDomain().getGhostedBox());
  // Tensor<Real,2> relaxres(iDomain.getDomain().getGhostedBox());
  // irPP.Laplacian(f,res,Func1,'D');
  // irPP.relaxLaplacian(res,f,relaxres,2.0/3);

  //IrBottomSolver<2> IBS(iDomain);
  
  //cout << f << endl;

  
  // cout << cf << endl;
  // cout << fcf << endl;
  // cout << res << endl;
  // cout << relaxres << endl;
  
  // const Vector<iVec>& vpoints = PS.getPoints(iVec{6,4});
  // for (auto p = vpoints.cbegin() ; p != vpoints.cend() ; ++p)
  //   cout << *p << endl;
  // const Vector<Real>& vweights = PS.getWeights(iVec{6,4});
  // for (auto p = vweights.cbegin() ; p != vweights.cend() ; ++p)
  //   cout << *p << endl;

  // PLG<2> plg(iDomain,iVec{4,2},'D');
  // const Vector<iVec>& vpoints = plg.getvpoints();
  // for (auto p = vpoints.cbegin() ; p != vpoints.cend() ; ++p)
  //   cout << *p << endl;
  // const Vector<Real>& vweights = plg.getvweights();
  // for (auto p = vweights.cbegin() ; p != vweights.cend() ; ++p)
  //   cout << *p << endl;
  // const Tensor<Real,2>& Matrix = plg.getMatrix();
  // cout << Matrix << endl;
  // Tensor<Real,1> rhs(Vec<int,1>{Matrix.size()[0]});
  // plg.computeRHS(f,rhs,Func1);
  // cout << rhs << endl;
  // plg.solvePLG(rhs);
  // const Vector<Real>& coes = plg.getcoes();
  //  for (auto p = coes.cbegin() ; p != coes.cend() ; ++p)
  //   cout << *p << endl;
}
