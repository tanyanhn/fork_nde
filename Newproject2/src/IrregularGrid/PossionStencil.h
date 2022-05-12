#ifndef _POSSIONSTENCIL_H_
#define _POSSIONSTENCIL_H_

#include "IrregularGrid/IrRectDomain.h"
#include <vector>

template <int Dim> class PossionStencil;

template <>
class PossionStencil<2>{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;

  template <class T>
  using Vector = std::vector<T>;

  PossionStencil(const IrRectDomain<2>& aDomain);

  const IrRectDomain<2>& getIrDomain() const;

  const Vector<iVec>& getPoints(const iVec& pt) const;

  const Vector<Real>& getWeights(const iVec& pt) const;
public:
  struct stencil{
    Vector<iVec> points;
    Vector<Real> weights;
  };
    

protected:
  IrRectDomain<2> Domain;
  Tensor<stencil,2> vstencil;
};

PossionStencil<2>::PossionStencil(const IrRectDomain<2>& aDomain):Domain(aDomain){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const rVec& dx = rDomain.spacing();
  vstencil.resize(rDomain);
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (Domain.isInside(Node)){
      if (Domain.isRePoint(Node)){
        vstencil(Node).points.push_back(Node);
        vstencil(Node).weights.push_back(-2.0/(dx[0]*dx[0]) - 2.0/(dx[1]*dx[1]));
        vstencil(Node).points.push_back(Node-iVec{1,0});
        vstencil(Node).weights.push_back(1.0/(dx[0]*dx[0]));
        vstencil(Node).points.push_back(Node+iVec{1,0});
        vstencil(Node).weights.push_back(1.0/(dx[0]*dx[0]));
        vstencil(Node).points.push_back(Node-iVec{0,1});
        vstencil(Node).weights.push_back(1.0/(dx[1]*dx[1]));
        vstencil(Node).points.push_back(Node+iVec{0,1});
        vstencil(Node).weights.push_back(1.0/(dx[1]*dx[1]));
      }
      else {
        // PLG<2> plg(Domain,Node);
        // vstencil(Node).points = plg.Laplacian().first;
        // vstencil(Node).weights = plg.Laplacian().second;
      }
    }
  }
}

const IrRectDomain<2>& PossionStencil<2>::getIrDomain() const{
  return Domain;
}

const std::vector<Vec<int,2> >& PossionStencil<2>::getPoints(const iVec& pt) const{
  return vstencil(pt).points;
}

const std::vector<Real>& PossionStencil<2>::getWeights(const iVec& pt) const{
  return vstencil(pt).weights;
}




#endif // _POSSIONSTENCIL_H_
