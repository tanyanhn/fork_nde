#ifndef _IRGHOSTFILLER_H_
#define _IRGHOSTFILLER_H_

#include "RegularGrid/RectDomain.h"
#include "RegularGrid/Func.h"
#include "Core/Tensor.h"


template <int Dim> class IrGhostFiller;

template <>
class IrGhostFiller<2>{
public:
  using iVec = Vec<int,2>;
  using rVec = Vec<Real,2>;
  
  IrGhostFiller(const IrRectDomain<2>& aDomain);

  void fillAllSide(Tensor<Real,2>& res, const char* BCType, const
                   Func<2>& func);

  void fillAllSide(Tensor<Real,2>& res, const char* BCType, int);
  
  // void fillIrSide(Tensor<Real,2>& res, const char BCType, const
  // Func<2>& func);

  //void fillIrSide(Tensor<Real,2>& res, const char BCType, int);
protected:
  void fillOneSide(Tensor<Real,2>& res, int dim, int side, const char
                   BCType, const Func<2>& func);

  void fillOneSide(Tensor<Real,2>& res, int dim, int side, const char
                   BCType, int);

  void fillIrSide1(Tensor<Real,2>& res, const char BCType, const
  Func<2>& func);

  void fillIrSide2(Tensor<Real,2>& res, const char BCType, const
  Func<2>& func);

  void fillIrSide1(Tensor<Real,2>& res, const char BCType, int);

  void fillIrSide2(Tensor<Real,2>& res, const char BCType, int);
  
protected:
  IrRectDomain<2> Domain;
};

IrGhostFiller<2>::IrGhostFiller(const IrRectDomain<2>& aDomain):Domain(aDomain){}


void IrGhostFiller<2>::fillOneSide(Tensor<Real,2>& res, int dim, int
                                 side, const char BCType, const
                                 Func<2>& func){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const iVec& lo = rDomain.lo();
  const iVec& sz = rDomain.size();
  const rVec& dx = rDomain.spacing(); 
  if (dim == 1 && side == -1){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[0] ; i++){
         iVec Node{i,0};
         rVec rNode = (Node - lo)*dx;
         res(0,i) = func.F(rNode);
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
         res(-1,i) = res(1,i)-2*dx[1]*func.Fy(rNode);
       }
     }
  }
  else if (dim == 1 && side == 1){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[0] ; i++){
         iVec Node{i,sz[1]-1};
         rVec rNode = (Node - lo)*dx;
         res(sz[1]-1,i) = func.F(rNode);
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
         res(sz[1],i) = res(sz[1]-2,i)+2*dx[1]*func.Fy(rNode);
       }
     }
  }
  else if (dim == 0 && side == -1){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[1] ; i++){
         iVec Node{0,i};
         rVec rNode = (Node - lo)*dx;
         res(i,0) = func.F(rNode);
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
         res(i,-1) = res(i,1)-2*dx[0]*func.Fx(rNode);
       }
     }
  }
  else if (dim == 0 && side == 1){
     if (BCType == 'D'){
       for(int i = 0 ; i < sz[1] ; i++){
         iVec Node{sz[0]-1,i};
         rVec rNode = (Node - lo)*dx;
         res(i,sz[0]-1) = func.F(rNode);
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
         res(i,sz[0]) = res(i,sz[0]-2)+2*dx[0]*func.Fx(rNode);
       }
     }
  }
  else
    assert(0);
}

void IrGhostFiller<2>::fillOneSide(Tensor<Real,2>& res, int dim, int side, const char
                                 BCType, int){
  D2Func0 f;
  fillOneSide(res,dim,side,BCType,f);
}

void IrGhostFiller<2>::fillIrSide1(Tensor<Real,2>& res, const char
                                   BCType, const Func<2>& func){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const iVec& lo = rDomain.lo();
  const iVec& sz = rDomain.size();
  const rVec& dx = rDomain.spacing();
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (!Domain.isInside(Node))
      res(Node) = 0.0;
  }
  for (int i = 0 ; i < sz[0] ; i++)
    res(i,-1) = 0.0;
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (Domain.isInside(Node) && (!Domain.isRePoint(Node))){
      rVec rpt = (Node - lo) * dx;
      iVec pt2 = Node + iVec{0,-1};
      rVec rpt2 = (pt2 - lo) * dx;
      if (!Domain.isInside(pt2)){
        Real y = (*Domain.getplbFunc())(rpt[0]);
        Real frac = (rpt[1] - y)/(rpt[1] - rpt2[1]);
        rVec bcpt{rpt[0],y};
        iVec tmpNode = Node + iVec{0,1};
        if (BCType == 'D'){
          if (frac < 1e-16){
            res(Node) = func.F(rpt);
            res(pt2) = 2*res(Node) - res(tmpNode);
          }
          else {
            //Real qfrac = frac*frac;
            //res(pt2) = (func.F(bcpt)-(1-qfrac)*res(Node)-0.5*(qfrac-frac)*res(tmpNode))/(0.5*(qfrac+frac));
            res(pt2) = (func.F(bcpt) - (1-frac)*res(Node))/frac;
          }
        }
      }  
    }
  }
}

void IrGhostFiller<2>::fillIrSide2(Tensor<Real,2>& res, const char
                                   BCType, const Func<2>& func){
  const RectDomain<2>& rDomain = Domain.getDomain();
  const iVec& lo = rDomain.lo();
  //const iVec& sz = rDomain.size();
  const rVec& dx = rDomain.spacing();
  loop_box_2(rDomain,i,j){
    iVec Node{i,j};
    if (Domain.isInside(Node) && (!Domain.isRePoint(Node))){
      rVec rpt = (Node - lo) * dx;
      iVec pt1 = Node + iVec{-1,0};
      rVec rpt1 = (pt1 - lo) * dx;
      if (!Domain.isInside(pt1) && fabs(res(pt1)) < 1e-16){
        Real x = (*Domain.getplbFunc()).rev2(rpt[1]);
        Real frac = (rpt[0] - x)/(rpt[0] - rpt1[0]);
        rVec bcpt{x,rpt[1]};
        iVec tmpNode = Node + iVec{1,0};
        if (BCType == 'D'){
          if (frac < 1e-16){
            res(Node) = func.F(rpt);
            res(pt1) = 2*res(Node) - res(tmpNode);
          }
          else {
            //Real qfrac = frac*frac;
            //res(pt1) = (func.F(bcpt)-(1-qfrac)*res(Node)-0.5*(qfrac-frac)*res(tmpNode))/(0.5*(qfrac+frac));
            res(pt1) = (func.F(bcpt) - (1-frac)*res(Node))/frac;
          }
        } 
      }
      iVec pt3 = Node + iVec{1,0};
      rVec rpt3 = (pt3 - lo) * dx;
      if (!Domain.isInside(pt3) && fabs(res(pt3)) < 1e-16){
        Real x = (*Domain.getplbFunc()).rev1(rpt[1]);
        Real frac = (x - rpt[0])/(rpt3[0] - rpt[0]);
        rVec bcpt{x,rpt[1]};
        iVec tmpNode = Node + iVec{-1,0};
        if (BCType == 'D'){
          if (frac < 1e-16){
            res(Node) = func.F(rpt);
            res(pt3) = 2*res(Node) - res(tmpNode);
          }
          else {
            //Real qfrac = frac*frac;
            //res(pt3) = (func.F(bcpt)-(1-qfrac)*res(Node)-0.5*(qfrac-frac)*res(tmpNode))/(0.5*(qfrac+frac));
            res(pt3) = (func.F(bcpt) - (1-frac)*res(Node))/frac;
          }
        } 
      }
    }
  }
}
// void IrGhostFiller<2>::fillIrSide(Tensor<Real,2>& res, const char
//                                   BCType, const Func<2>& func){
//   const RectDomain<2>& rDomain = Domain.getDomain();
//   const iVec& lo = rDomain.lo();
//   const iVec& sz = rDomain.size();
//   const rVec& dx = rDomain.spacing();
//   loop_box_2(rDomain,i,j){
//     iVec Node{i,j};
//     if (!Domain.isInside(Node))
//       res(Node) = 0.0;
//   }
//   for (int i = 0 ; i < sz[0] ; i++)
//     res(i,-1) = 0.0;
//   loop_box_2(rDomain,i,j){
//     iVec Node{i,j};
//     if (Domain.isInside(Node) && (!Domain.isRePoint(Node))){
//       rVec rpt = (Node - lo) * dx;
//       iVec pt1 = Node + iVec{-1,0};
//       rVec rpt1 = (pt1 - lo) * dx;
//       if (!Domain.isInside(pt1)){
//         Real x = (*Domain.getplbFunc()).rev2(rpt[1]);
//         Real frac = (rpt[0] - x)/(rpt[0] - rpt1[0]);
//         rVec bcpt{x,rpt[1]};
//         iVec tmpNode = Node + iVec{1,0};
//         if (BCType == 'D'){
//           if (frac < 1e-16){
//             res(Node) = func.F(rpt);
//             res(pt1) = 2*res(Node) - res(tmpNode);
//           }
//           else {
//             Real qfrac = frac*frac;
//             res(pt1) = (func.F(bcpt)-(1-qfrac)*res(Node)-0.5*(qfrac-frac)*res(tmpNode))/(0.5*(qfrac+frac));
//           }
//         } 
//       }
//       iVec pt2 = Node + iVec{0,-1};
//       rVec rpt2 = (pt2 - lo) * dx;
//       if (!Domain.isInside(pt2)){
//         Real y = (*Domain.getplbFunc())(rpt[0]);
//         Real frac = (rpt[1] - y)/(rpt[1] - rpt2[1]);
//         rVec bcpt{rpt[0],y};
//         iVec tmpNode = Node + iVec{0,1};
//         if (BCType == 'D'){
//           if (frac < 1e-16){
//             res(Node) = func.F(rpt);
//             res(pt2) = 2*res(Node) - res(tmpNode);
//           }
//           else {
//             Real qfrac = frac*frac;
//             res(pt2) = (func.F(bcpt)-(1-qfrac)*res(Node)-0.5*(qfrac-frac)*res(tmpNode))/(0.5*(qfrac+frac));
//           }
//         }
//       }
//       iVec pt3 = Node + iVec{1,0};
//       rVec rpt3 = (pt3 - lo) * dx;
//       if (!Domain.isInside(pt3)){
//         Real x = (*Domain.getplbFunc()).rev1(rpt[1]);
//         Real frac = (x - rpt[0])/(rpt3[0] - rpt[0]);
//         rVec bcpt{x,rpt[1]};
//         iVec tmpNode = Node + iVec{-1,0};
//         if (BCType == 'D'){
//           if (frac < 1e-16){
//             res(Node) = func.F(rpt);
//             res(pt3) = 2*res(Node) - res(tmpNode);
//           }
//           else {
//             Real qfrac = frac*frac;
//             res(pt3) = (func.F(bcpt)-(1-qfrac)*res(Node)-0.5*(qfrac-frac)*res(tmpNode))/(0.5*(qfrac+frac));
//           }
//         } 
//       }
//     }
//   }
// }

void IrGhostFiller<2>::fillIrSide1(Tensor<Real,2>& res, const char BCType, int){
  D2Func0 f;
  fillIrSide1(res,BCType,f);
}

void IrGhostFiller<2>::fillIrSide2(Tensor<Real,2>& res, const char BCType, int){
  D2Func0 f;
  fillIrSide2(res,BCType,f);
}

void IrGhostFiller<2>::fillAllSide(Tensor<Real,2>& res, const char*
                                 BCType, const Func<2>& func){
  fillIrSide1(res,BCType[0],func);
  for(int i = 1 ; i < 4 ; i++){
    if (BCType[i] == 'N')
      fillOneSide(res,i/2,-1+(i%2)*2,BCType[i],func);
  }
  for(int i = 1 ; i < 4 ; i++){
    if (BCType[i] != 'N')
      fillOneSide(res,i/2,-1+(i%2)*2,BCType[i],func);
  }
  fillIrSide2(res,BCType[0],func);
}

void IrGhostFiller<2>::fillAllSide(Tensor<Real,2>& res, const char*
                                 BCType, int){
  fillIrSide1(res,BCType[0],0);
  for(int i = 1 ; i < 4 ; i++){
    if (BCType[i] == 'N')
      fillOneSide(res,i/2,-1+(i%2)*2,BCType[i],0);
  }
  for(int i = 1 ; i < 4 ; i++){
    if (BCType[i] != 'N')
      fillOneSide(res,i/2,-1+(i%2)*2,BCType[i],0);
  }
  fillIrSide2(res,BCType[0],0);
}


#endif //_IRGHOSTFILLER_H_
