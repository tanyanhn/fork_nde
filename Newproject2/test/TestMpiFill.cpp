#include "Multigrid/MGFactory.h"
#include "Multigrid/TestMultigrid.h"
#include <iostream>
#include <mpi/mpi.h>

using namespace std;

using iVec = Vec<int,2>;
using rVec = Vec<Real,2>;
template <class T>
using Vector = std::vector<T>;

int printcore = 2;

using namespace std;
template <class T>
void printCuba(T* A, int m, int n, int a = -1, int b = -1) {
  cout.setf(std::ios::fixed);
  if (a == -1 && b == -1) {
    for (auto i = m-1; i > -1; --i) {
      for (auto j = 0; j < n; ++j)
        cout << A[i * n + j] << ' ';
      cout << '\n';
    }
  } else if (a == -1 && b != -1) {
    for (auto i = m-1; i > -1; --i) {
      int j = b;
      cout << A[i * n + j] << ' ';
      cout << '\n';
    }
  } else if (a != -1 && b == -1) {
    int i = a;
    for (auto j = 0; j < n; ++j)
      cout << A[i * n + j] << ' ';
    cout << '\n';
  } else {
    cout << A[a * n + b] << '\n';
  }
}


template<class T>
void fillCuba(T* A, int m, int n, T va) {
  for (auto i = 0; i < m; ++i) 
    for(auto j = 0; j < n; ++j)
      A[i * n + j] = va + i * 1000 + j * 100;
}

int main() {
  MPI_Init(nullptr, nullptr);
  int mpiSize;
  int mpiRank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  int l = 4 + 2;
  Box<2> bx(0, l - 2 -1);
  Tensor<double, 2> L(bx.inflate(1));
  fillCuba<double>(L.data(), l, l, mpiRank);

  RectDomain<2> domain(bx, 1.0 / (l - 2), NodeCentered, 1);
  GhostFiller<2> filler(domain);
  filler.transformInitMpi();
  ScalarFunction<2>* pD2Func1f = new D2Func1f;
  ScalarFunction<2>* pD2Func1F = new D2FuncTest;
  ScalarFunction<2>* pD2Func1Fx = new D2Func1Fx;
  ScalarFunction<2>* pD2Func1Fy = new D2Func1Fy;
  std::array<ScalarFunction<2>*,3> pfuncs{pD2Func1F,pD2Func1Fx,pD2Func1Fy};
  char BCType[5] = "DDDD";
  if (mpiRank == 0)
    BCType[1] = 'G', BCType[3] = 'G';
  else if (mpiRank == 1)
    BCType[1] = 'G', BCType[2] = 'G';
  else if (mpiRank == 2)
    BCType[0] = 'G', BCType[3] = 'G';
  else if (mpiRank == 1)
    BCType[0] = 'G', BCType[2] = 'G';

  if (mpiRank == printcore) {
    printCuba<double>(L.data(), l, l);
    cout << endl;
  }
  filler.fillAllSidesMpi(L, BCType, pfuncs);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpiRank == printcore) {
    printCuba<double>(L.data(), l, l);
    cout << endl;
  }
  filler.transformFreeMpi();
  MPI_Finalize();
}