#ifndef _GHOSTFILLER_H_
#define _GHOSTFILLER_H_

#include "RegularGrid/ScalarFunction.h"
#include "RegularGrid/RectDomain.h"
#include "Core/Tensor.h"
#include <array>


template <int Dim> class GhostFiller;

template <int Dim>
class GhostFiller{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;
  
  GhostFiller(const RectDomain<Dim>& adomain);
  
  void fillAllSides(Tensor<Real,Dim>& res, const char* const BCTypes, const
                    std::array<ScalarFunction<Dim>*,3> pfuncs =
                    std::array<ScalarFunction<Dim>*,3>{}) const;
public:
  enum side{
            low = -1,
            high = 1
  };
protected:
  void fillOneSide(Tensor<Real,Dim>& res, int dim, side s, const char
                   BCType, const ScalarFunction<Dim>* pfunc = nullptr)
    const;

protected:
  RectDomain<Dim> domain;
};

template <int Dim>
GhostFiller<Dim>::GhostFiller(const RectDomain<Dim>& adomain){
  assert(adomain.getCentering() == NodeCentered);
  assert(adomain.getNumGhost() == 1);
  domain = adomain;
}

// Case Dim = 2 -----------------------------------------------

template<>
void GhostFiller<2>::fillOneSide(Tensor<Real,2>& res, int dim, side
                                 s, const char BCType, const
                                 ScalarFunction<2>* pfunc) const{
  const Box<2>& bx = domain;
  const iVec& lo = bx.lo();
  const iVec& sz = bx.size();
  const rVec& dx = domain.spacing();
  if (pfunc == nullptr){
    ScalarFunction<2>* pD2func0 = new D2Func0;
    fillOneSide(res,dim,s,BCType,pD2func0);
    return;
  }
  if (dim == 1 && s == low){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[0] ; i++){
         iVec Node{i,0};
         rVec rNode = (Node - lo)*dx;
         res(0,i) = (*pfunc)(rNode);
         res(-1,i) = 6*res(1,i)-8*res(2,i)+3*res(3,i);
       }
     }
     else if (BCType == 'X'){
       for(int i = 0 ; i < sz[0] ; i++){
         res(-1,i) = 6*res(1,i)-8*res(2,i)+3*res(3,i);
       }
     }
     else if (BCType == 'N'){
       for(int i = 0 ; i < sz[0] ; i++){
         iVec Node{i,0};
         rVec rNode = (Node - lo)*dx;
         res(-1,i) = res(1,i)-2*dx[1]*(*pfunc)(rNode);
       }
     }
  }
  else if (dim == 1 && s == high){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[0] ; i++){
         iVec Node{i,sz[1]-1};
         rVec rNode = (Node - lo)*dx;
         res(sz[1]-1,i) = (*pfunc)(rNode);
         res(sz[1],i) = 6*res(sz[1]-2,i)-8*res(sz[1]-3,i)+3*res(sz[1]-4,i);
       }
     }
     else if (BCType == 'X'){
       for(int i = 0 ; i < sz[0] ; i++){
         res(sz[1],i) = 6*res(sz[1]-2,i)-8*res(sz[1]-3,i)+3*res(sz[1]-4,i);
       }
     }
     else if (BCType == 'N'){
       for(int i = 0 ; i < sz[0] ; i++){
         iVec Node{i,sz[1]-1};
         rVec rNode = (Node - lo)*dx;
         res(sz[1],i) = res(sz[1]-2,i)+2*dx[1]*(*pfunc)(rNode);
       }
     }
  }
  else if (dim == 0 && s == low){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[1] ; i++){
         iVec Node{0,i};
         rVec rNode = (Node - lo)*dx;
         res(i,0) = (*pfunc)(rNode);
         res(i,-1) = 6*res(i,1)-8*res(i,2)+3*res(i,3);
       }
     }
     else if (BCType == 'X'){
       for(int i = 0 ; i < sz[1] ; i++){
         res(i,-1) = 6*res(i,1)-8*res(i,2)+3*res(i,3);
       }
     }
     else if (BCType == 'N'){
       for(int i = 0 ; i < sz[1] ; i++){
         iVec Node{0,i};
         rVec rNode = (Node - lo)*dx;
         res(i,-1) = res(i,1)-2*dx[0]*(*pfunc)(rNode);
       }
     }
  }
  else if (dim == 0 && s == high){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[1] ; i++){
         iVec Node{sz[0]-1,i};
         rVec rNode = (Node - lo)*dx;
         res(i,sz[0]-1) = (*pfunc)(rNode);
         res(i,sz[0]) = 6*res(i,sz[0]-2)-8*res(i,sz[0]-3)+3*res(i,sz[0]-4);
       }
     }
     else if (BCType == 'X'){
       for(int i = 0 ; i < sz[1] ; i++){
         res(i,sz[0]) = 6*res(i,sz[0]-2)-8*res(i,sz[0]-3)+3*res(i,sz[0]-4);
       }
     }
     else if (BCType == 'N'){
       for(int i = 0 ; i < sz[1] ; i++){
         iVec Node{sz[0]-1,i};
         rVec rNode = (Node - lo)*dx;
         res(i,sz[0]) = res(i,sz[0]-2)+2*dx[0]*(*pfunc)(rNode);
       }
     }
  }
  else
    assert(0);
}


template<>
void GhostFiller<2>::fillAllSides(Tensor<Real,2>& res, const char* const
                                  BCTypes , const
                                  std::array<ScalarFunction<2>*,3>
                                  pfuncs) const{
  for(int i = 0 ; i < 4 ; i++){
    if (BCTypes[i] == 'N')
      fillOneSide(res,i/2,(-1+(i%2)*2 == 1)?high:low,BCTypes[i],pfuncs[i/2 + 1]);
  }
  for(int i = 0 ; i < 4 ; i++){
    if (BCTypes[i] != 'N')
      fillOneSide(res,i/2,(-1+(i%2)*2 == 1)?high:low,BCTypes[i],pfuncs[0]);
  }
}


// Case Dim = 1 -----------------------------------------------

template<>
void GhostFiller<1>::fillOneSide(Tensor<Real,1>& res, int dim, side
                                 s, const char BCType, const
                                 ScalarFunction<1>* pfunc) const{
  assert(dim == 0);
  const Box<1>& bx = domain;
  const Vec<int,1>& lo = bx.lo();
  const Vec<int,1>& sz = bx.size();
  const Vec<Real,1>& dx = domain.spacing();
  if (pfunc == nullptr){
    ScalarFunction<1>* pD1func0 = new D1Func0;
    fillOneSide(res,dim,s,BCType,pD1func0);
    return;
  }
  if (s == low){
    if (BCType == 'D'){
      Vec<int,1> Node{0};
      Vec<Real,1> rNode = (Node - lo)*dx;
      res(0) = (*pfunc)(rNode);
      res(-1) = 6*res(1)-8*res(2)+3*res(3);
    }
    else if (BCType == 'X'){
      res(-1) = 6*res(1)-8*res(2)+3*res(3);
    }
    else if (BCType == 'N'){
      for(int i = 0 ; i < sz[0] ; i++){
        Vec<int,1> Node{0};
        Vec<Real,1> rNode = (Node - lo)*dx;
        res(-1) = res(1)-2*dx[0]*(*pfunc)(rNode);
      }
    }
  }
  else if (s == high){
    if (BCType == 'D'){
      Vec<int,1> Node{sz[0]-1};
      Vec<Real,1> rNode = (Node - lo)*dx;
      res(sz[0]-1) = (*pfunc)(rNode);
      res(sz[0]) = 6*res(sz[0]-2)-8*res(sz[0]-3)+3*res(sz[0]-4);
    }
    else if (BCType == 'X'){
      res(sz[0]) = 6*res(sz[0]-2)-8*res(sz[0]-3)+3*res(sz[0]-4);
    }
    else if (BCType == 'N'){
      Vec<int,1> Node{sz[0]-1};
      Vec<Real,1> rNode = (Node - lo)*dx;
      res(sz[0]) = res(sz[0]-2)+2*dx[0]*(*pfunc)(rNode);
     }
  }
  else
    assert(0);
}

template<>
void GhostFiller<1>::fillAllSides(Tensor<Real,1>& res, const char* const
                                  BCTypes, const
                                  std::array<ScalarFunction<1>*,3>
                                  pfuncs) const{
  if (BCTypes[0] == 'N')
    fillOneSide(res,0,low,BCTypes[0],pfuncs[1]);
  else
    fillOneSide(res,0,low,BCTypes[0],pfuncs[0]);
  if (BCTypes[1] == 'N')
    fillOneSide(res,0,high,BCTypes[1],pfuncs[1]);
  else
    fillOneSide(res,0,high,BCTypes[1],pfuncs[0]);
}

#endif //_GHOSTFILLER_H_
