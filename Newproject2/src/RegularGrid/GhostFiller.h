#ifndef _GHOSTFILLER_H_
#define _GHOSTFILLER_H_

#include "RegularGrid/ScalarFunction.h"
#include "RegularGrid/RectDomain.h"
#include "Core/Tensor.h"
#include <mpi/mpi.h>
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
  void fillAllSidesMpi(Tensor<Real,Dim>& res, const char* const BCTypes, const
                    std::array<ScalarFunction<Dim>*,3> pfuncs =
                    std::array<ScalarFunction<Dim>*,3>{}) const;

  void transformInitMpi();
  void transformFreeMpi();

 public:
  enum side{
            low = -1,
            high = 1
  };
protected:
 void fillOneSide(Tensor<Real, Dim>& res,
                  int dim,
                  side s,
                  const char BCType,
                  const ScalarFunction<Dim>* pfunc = nullptr) const;

protected:
 RectDomain<Dim> domain;
 MPI_Datatype dataTypeMpi[Dim];
 int neighborCore[2 * Dim] = {-1};
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
  const iVec& hi = bx.hi();
  const rVec& dx = domain.spacing();
  if (pfunc == nullptr){
    ScalarFunction<2>* pD2func0 = new D2Func0;
    fillOneSide(res,dim,s,BCType,pD2func0);
    return;
  }
  if (dim == 1 && s == low){
     if (BCType == 'D'){
       for(int i = lo[1]; i <= hi[1] ; i++){
         iVec Node{lo[0], i};
         rVec rNode = (Node)*dx;
         res(lo[0], i) = (*pfunc)(rNode);  
         res(lo[0]-1,i) = 6*res(lo[0] + 1,i)-8*res(lo[0] + 2,i)+3*res(lo[0] + 3,i);
       }
     }
     else if (BCType == 'N'){
       for(int i = lo[1]; i <= hi[1] ; i++){
         iVec Node{lo[0], i};
         rVec rNode = (Node)*dx;
         res(lo[0]-1, i) = res(1,i)-2*dx[1]*(*pfunc)(rNode);
       }
     }
  }
  else if (dim == 1 && s == high){
     if (BCType == 'D'){
       for(int i = lo[1]; i <= hi[1] ; i++){
         iVec Node{hi[0], i};
         rVec rNode = (Node)*dx;
         res(hi[0], i) = (*pfunc)(rNode);
         res(hi[0]+1,i) = 6*res(hi[0] - 1,i)-8*res(hi[0] - 2,i)+3*res(hi[0] - 3,i);
       }
     }
     else if (BCType == 'N'){
       for(int i = lo[1] ; i < hi[1]; i++){
         iVec Node{hi[0], i};
         rVec rNode = (Node)*dx;
         res(hi[0] + 1,i) = res(hi[0] - 1,i)+2*dx[1]*(*pfunc)(rNode);
       }
     }
  }
  else if (dim == 0 && s == low){
     if (BCType == 'D'){
       for(int i = lo[0]; i <= hi[0] ; i++){
         iVec Node{i, lo[1]};
         rVec rNode = (Node)*dx;
         res(i, lo[1]) = (*pfunc)(rNode);
         res(i, lo[1]-1) = 6*res(i, lo[1] + 1)-8*res(i, lo[1] + 2)+3*res(i, lo[1] + 3);
       }
     }
     else if (BCType == 'N'){
       for(int i = lo[0]; i <= hi[0] ; i++){
         iVec Node{i, lo[1]};
         rVec rNode = (Node)*dx;
         res(i, lo[1]-1) = res(1,i)-2*dx[0]*(*pfunc)(rNode);
       }
     }
  }
  else if (dim == 0 && s == high){
     if (BCType == 'D'){
       for(int i = lo[0]; i <= hi[0] ; i++){
         iVec Node{i, hi[1]};
         rVec rNode = (Node)*dx;
         res(i, hi[1]) = (*pfunc)(rNode);
         res(i, hi[1]+1) = 6*res(i, hi[1] - 1)-8*res(i, hi[1] - 2)+3*res(i, hi[1] - 3);
       }
     }
     else if (BCType == 'N'){
       for(int i = lo[0] ; i < hi[0]; i++){
         iVec Node{hi[1], i};
         rVec rNode = (Node)*dx;
         res(i, hi[1] + 1) = res(hi[1] - 1,i)+2*dx[0]*(*pfunc)(rNode);
       }
     }
  }
  else
    assert(0);
}


template<>
void GhostFiller<2>::transformInitMpi() {
    int mpiSize;
    int mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    int nx = std::round(std::pow(mpiSize, 1.0 / 2));
    int i = mpiRank % nx;
    int j = mpiRank / nx;
    neighborCore[1] = mpiRank + nx;
    neighborCore[0] = mpiRank - nx;
    neighborCore[3] = mpiRank + 1;
    neighborCore[2] = mpiRank - 1;
    if (j == nx - 1)
      neighborCore[1] = MPI_PROC_NULL;
    else if (j == 0)
      neighborCore[0] = MPI_PROC_NULL;
    if (i == nx - 1)
      neighborCore[3] = MPI_PROC_NULL;
    else if (i == 0)
      neighborCore[2] = MPI_PROC_NULL;

    const auto& bx = domain.inflate(1);
    MPI_Type_vector(bx.size()[0] - 2 * (domain.getNumGhost()), 1, 1, MPI_DOUBLE, &dataTypeMpi[0]);
    MPI_Type_vector(bx.size()[1] - 2 * (domain.getNumGhost()), 1, bx.size()[0], MPI_DOUBLE, &dataTypeMpi[1]);
    MPI_Type_commit(&dataTypeMpi[0]);
    MPI_Type_commit(&dataTypeMpi[1]);
}


template <>
void GhostFiller<2>::fillAllSides(
    Tensor<Real, 2>& res,
    const char* const BCTypes,
    const std::array<ScalarFunction<2>*, 3> pfuncs) const {
  for (int i = 0; i < 4; i++) {
    if (BCTypes[i] == 'N')
      fillOneSide(res, i / 2, (-1 + (i % 2) * 2 == 1) ? high : low, BCTypes[i],
                  pfuncs[2 - i / 2]);
  }
  for (int i = 0; i < 4; i++) {
    if (BCTypes[i] == 'D')
      fillOneSide(res, i / 2, (-1 + (i % 2) * 2 == 1) ? high : low, BCTypes[i],
                  pfuncs[0]);
  }
}
template <>
void GhostFiller<2>::fillAllSidesMpi(
    Tensor<Real, 2>& res,
    const char* const BCTypes,
    const std::array<ScalarFunction<2>*, 3> pfuncs) const {
  for (int i = 0; i < 4; ++i) {
    if (BCTypes[i] == 'N') {
      fillOneSide(res, i / 2, (-1 + (i % 2) * 2 == 1) ? high : low, BCTypes[i],
                  pfuncs[2 - i / 2]);
    } else if (BCTypes[i] == 'D') {
      fillOneSide(res, i / 2, (-1 + (i % 2) * 2 == 1) ? high : low, BCTypes[i],
                  pfuncs[0]);
    }
  }
  const auto& bx = res.box();
  int nghost = domain.getNumGhost();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Sendrecv(&res(bx.lo()[0] + nghost, bx.lo()[1] + 2), 1, dataTypeMpi[0],
                 neighborCore[0], rank, &res(bx.lo()[0] + nghost, bx.hi()[1]),
                 1, dataTypeMpi[0], neighborCore[1], neighborCore[1],
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&res(bx.lo()[0] + nghost, bx.hi()[1] - 2), 1, dataTypeMpi[0],
                 neighborCore[1], rank, &res(bx.lo()[0] + nghost, bx.lo()[1]),
                 1, dataTypeMpi[0], neighborCore[0], neighborCore[0],
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&res(bx.lo()[0] + 2, bx.lo()[1] + nghost), 1, dataTypeMpi[1],
                 neighborCore[2], rank, &res(bx.hi()[0], bx.lo()[1] + nghost),
                 1, dataTypeMpi[1], neighborCore[3], neighborCore[3],
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&res(bx.hi()[0] - 2, bx.lo()[1] + nghost), 1, dataTypeMpi[1],
                 neighborCore[3], rank, &res(bx.lo()[0], bx.lo()[1] + nghost),
                 1, dataTypeMpi[1], neighborCore[2], neighborCore[2],
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

template<int Dim>
void GhostFiller<Dim>::transformFreeMpi(){
  MPI_Type_free(&dataTypeMpi[0]);
  if constexpr(Dim > 1)
    MPI_Type_free(&dataTypeMpi[1]);
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
      Vec<Real,1> rNode = (Node)*dx;
      res(0) = (*pfunc)(rNode);
      res(-1) = 6*res(1)-8*res(2)+3*res(3);
    }
    else if (BCType == 'X'){
      res(-1) = 6*res(1)-8*res(2)+3*res(3);
    }
    else if (BCType == 'N'){
      for(int i = 0 ; i < sz[0] ; i++){
        Vec<int,1> Node{0};
        Vec<Real,1> rNode = (Node)*dx;
        res(-1) = res(1)-2*dx[0]*(*pfunc)(rNode);
      }
    }
  }
  else if (s == high){
    if (BCType == 'D'){
      Vec<int,1> Node{sz[0]-1};
      Vec<Real,1> rNode = (Node)*dx;
      res(sz[0]-1) = (*pfunc)(rNode);
      res(sz[0]) = 6*res(sz[0]-2)-8*res(sz[0]-3)+3*res(sz[0]-4);
    }
    else if (BCType == 'X'){
      res(sz[0]) = 6*res(sz[0]-2)-8*res(sz[0]-3)+3*res(sz[0]-4);
    }
    else if (BCType == 'N'){
      Vec<int,1> Node{sz[0]-1};
      Vec<Real,1> rNode = (Node)*dx;
      res(sz[0]) = res(sz[0]-2)+2*dx[0]*(*pfunc)(rNode);
     }
  }
  else
    assert(0);
}

template<>
void GhostFiller<1>::transformInitMpi() {
    int mpiSize;
    int mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank == 0)
      neighborCore[0] = MPI_PROC_NULL;
    else
      neighborCore[0] = mpiRank - 1;
    if (mpiRank == mpiSize - 1)
      neighborCore[1] = MPI_PROC_NULL;
    else 
      neighborCore[1] = mpiRank + 1;
  MPI_Type_vector(1, 1, 1, MPI_DOUBLE, &dataTypeMpi[0]);
  MPI_Type_commit(&dataTypeMpi[0]);
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

template <>
void GhostFiller<1>::fillAllSidesMpi(
    Tensor<Real, 1>& res,
    const char* const BCTypes,
    const std::array<ScalarFunction<1>*, 3> pfuncs) const {
  for (int i = 0; i < 2; ++i) {
    if (BCTypes[i] == 'N') {
    fillOneSide(res,0,i == 0 ? low : high,BCTypes[i],pfuncs[1]);
    } else if (BCTypes[i] == 'D') {
    fillOneSide(res,0,i == 0 ? low : high,BCTypes[i],pfuncs[0]);
    } 
  }
  const auto& bx = res.box();
  MPI_Sendrecv(&res(bx.lo()[0] + 2), 1, dataTypeMpi[0], neighborCore[0], 0,
               &res(bx.hi()[0]), 1, dataTypeMpi[0], neighborCore[1], 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&res(bx.hi()[0] - 2), 1, dataTypeMpi[0], neighborCore[1], 0,
               &res(bx.lo()[0]), 1, dataTypeMpi[0], neighborCore[0], 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


#endif //_GHOSTFILLER_H_
