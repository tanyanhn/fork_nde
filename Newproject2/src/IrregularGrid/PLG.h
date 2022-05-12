#ifndef _PLG_H_
#define _PLG_H_

#include "IrregularGrid/IrRectDomain.h"
#include "IrregularGrid/QRHelper.h"
#include <utility>

template <int Dim> class PLG;

template <>
class PLG<2>{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;

  template <class T>
  using Vector = std::vector<T>;

  PLG(const IrRectDomain<2>& aDomain, const iVec& apt, const char aBCType);

  const Vector<iVec>& getvpoints() const;

  const Tensor<Real,2>& getMatrix() const;

  const Vector<Real>& getcoes() const;

  const Vector<Real>& getvweights() const;

  template <class TFunc>
  Real applyLaplacian(const Tensor<Real,2>& phi, const TFunc& func, const char BCType);

  //protected:
  void setPoints();
  void setWeights();
  void constructMatrix();
  
  template <class TFunc>
  void computeRHS(const Tensor<Real,2>& phi, Tensor<Real,1>& rhs,
                  const TFunc& func);

  void solvePLG(const Tensor<Real,1>& rhs);
protected:
  IrRectDomain<2> Domain;
  iVec targetpt;
  char BCType;
  Real frac[3];
  Vector<Real> coes;
  Vector<iVec> vpoints;
  Vector<Real> vweights;
  Tensor<Real,2> Matrix;
};

PLG<2>::PLG(const IrRectDomain<2>& aDomain, const iVec& apt, const
            char aBCType):Domain(aDomain),targetpt(apt),BCType(aBCType){
  assert(Domain.isRePoint(targetpt) == 0);
  setPoints();
  setWeights();
  constructMatrix();
}

const std::vector<Vec<int,2> >& PLG<2>::getvpoints() const{
  return vpoints;
}

const Tensor<Real,2>& PLG<2>::getMatrix() const{
  return Matrix;
}

const std::vector<Real>& PLG<2>::getcoes() const{
  return coes;
}

const std::vector<Real>& PLG<2>::getvweights() const{
  return vweights;
}

template <class TFunc>
Real PLG<2>::applyLaplacian(const Tensor<Real,2>& phi, const TFunc&
                            func, const char BCType){
  Tensor<Real,1> rhs(Vec<int,1>{Matrix.size()[0]});
  computeRHS(phi,rhs,func);
  solvePLG(rhs);
  return 2*coes[0]+2*coes[1];
}

void PLG<2>::setPoints(){
  vpoints.push_back(targetpt);
  iVec pt1 = targetpt + iVec{0,-1};
  iVec pt2 = targetpt + iVec{1,-1};
  iVec pt3 = targetpt + iVec{1,0};
  iVec pt4 = targetpt + iVec{1,1};
  iVec pt5 = targetpt + iVec{0,1};
  iVec pt6 = targetpt + iVec{-1,1};
  iVec pt7 = targetpt + iVec{-1,0};
  iVec pt8 = targetpt + iVec{-1,-1};
  if (Domain.isInside(pt1))
    vpoints.push_back(pt1);
  if (Domain.isInside(pt2))
    vpoints.push_back(pt2);
  if (Domain.isInside(pt3))
    vpoints.push_back(pt3);
  if (Domain.isInside(pt4))
    vpoints.push_back(pt4);
  if (Domain.isInside(pt5))
    vpoints.push_back(pt5);
  if (Domain.isInside(pt6))
    vpoints.push_back(pt6);
  if (Domain.isInside(pt7))
    vpoints.push_back(pt7);
  if (Domain.isInside(pt8))
    vpoints.push_back(pt8);
  int size = vpoints.size();
  if (size < 6){
    iVec pt9 = targetpt + iVec{0,2};
    if (Domain.isInside(pt9))
      vpoints.push_back(pt9);
    size++;
  }
}

void PLG<2>::setWeights(){
  const int size = vpoints.size();
  for (int i = 0 ; i < size ; i++){
    rVec tmp = (rVec)(vpoints[i] - targetpt);
    Real distance = norm(tmp,2);
    vweights.push_back(1.0/(1.0+distance));
  }
}

void PLG<2>::constructMatrix(){
  const iVec& lo = Domain.getDomain().lo();
  const rVec& dx = Domain.getDomain().spacing();
  rVec rpt = (targetpt - lo) * dx;
  iVec pt1 = targetpt + iVec{0,-1};
  iVec pt2 = targetpt + iVec{1,0};
  iVec pt3 = targetpt + iVec{-1,0};
  int size = vpoints.size();
  frac[0] = frac[1] = frac[2] = -1;
  if (!Domain.isInside(pt1)){
    rVec rpt1 = (pt1 - lo) * dx;
    Real y = (*Domain.getplbFunc())(rpt[0]);
    frac[0] = (rpt[1] - y)/(rpt[1] - rpt1[1]);
  }
  if (!Domain.isInside(pt2)){
    rVec rpt2 = (pt2 - lo) * dx;
    Real x = (*Domain.getplbFunc()).rev1(rpt[1]);
    frac[1] = (x - rpt[0])/(rpt2[0] - rpt[0]);
  }
  if (!Domain.isInside(pt3)){
    rVec rpt3 = (pt3 - lo) * dx;
    Real x = (*Domain.getplbFunc()).rev2(rpt[1]);
    frac[2] = (rpt[0] - x)/(rpt[0] - rpt3[0]);
  }
  int count = 0;
  for (int i = 0 ; i < 3 ; i++){
    if (frac[i] > 0)
      count++;
  }
  std::cout << frac[0] << "," << frac[1] << "," << frac[2] <<
  std::endl;
  Matrix.resize(iVec{size+count,6});
  for (int i = 0 ; i < size ; i++){
    int x = vpoints[i][0] - targetpt[0];
    int y = vpoints[i][1] - targetpt[1];
    Matrix(i,0) = x*x*vweights[i];Matrix(i,1) = y*y*vweights[i];
    Matrix(i,2) = x*y*vweights[i];Matrix(i,3) = x*vweights[i];
    Matrix(i,4) = y*vweights[i];Matrix(i,5) = vweights[i];
  }
  int k = size;
  if (frac[0] > 0){
    Real x = 0;
    Real y = -frac[0];
    if (BCType == 'D'){
      Matrix(k,0) = x*x;Matrix(k,1) = y*y;Matrix(k,2) = x*y;
      Matrix(k,3) = x;Matrix(k,4) = y;Matrix(k,5) = 1;
    } 
    k++;
  }
  if (frac[1] > 0){
    Real x = frac[1];
    Real y = 0;
    if (BCType == 'D'){
      Matrix(k,0) = x*x;Matrix(k,1) = y*y;Matrix(k,2) = x*y;
      Matrix(k,3) = x;Matrix(k,4) = y;Matrix(k,5) = 1;
    }
    k++;
  }
  if (frac[2] > 0){
    Real x = -frac[2];
    Real y = 0;
    if (BCType == 'D'){
      Matrix(k,0) = x*x;Matrix(k,1) = y*y;Matrix(k,2) = x*y;
      Matrix(k,3) = x;Matrix(k,4) = y;Matrix(k,5) = 1;
    }
    k++;
  }
 
}

template <class TFunc>
void PLG<2>::computeRHS(const Tensor<Real,2>& phi, Tensor<Real,1>&
                        rhs, const TFunc& func){
  const iVec& lo = Domain.getDomain().lo();
  const rVec& dx = Domain.getDomain().spacing();
  const int size = vpoints.size();
  int count = 0;
  for (int i = 0 ; i < 3 ; i++){
    if (frac[i] > 0)
      count++;
  }
  for (int i = 0 ; i < size ; i++){
    rhs(i) = phi(vpoints[i])*vweights[i];
  }
  int k = size;
  if (frac[0] > 0){
    Real x = 0;
    Real y = -frac[0];
    rVec rpt = (targetpt + rVec{x,y} - lo)*dx;
    if (BCType == 'D')
      rhs(k) = func.F(rpt);
    k++;
  }
  if (frac[1] > 0){
    Real x = frac[1];
    Real y = 0;
    rVec rpt = (targetpt + rVec{x,y} - lo)*dx;
    if (BCType == 'D')
      rhs(k) = func.F(rpt);
    k++;
  }
  if (frac[2] > 0){
    Real x = -frac[2];
    Real y = 0;
    rVec rpt = (targetpt + rVec{x,y} - lo)*dx;
    if (BCType == 'D')
      rhs(k) = func.F(rpt);
    k++;
  }
   if (frac[0] == 0 || frac[1] == 0 || frac[2] == 0){
     rVec rpt = (targetpt - lo) * dx;
     rhs(0) = func.F(rpt)*vweights[0];
   }
    
}

void PLG<2>::solvePLG(const Tensor<Real,1>& rhs){
  coes = QRSolve(Matrix.data(), rhs.data(),
  Matrix.size()[0], 6);
}

#endif // _PLG_H_
