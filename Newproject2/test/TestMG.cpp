#include "Multigrid/MGFactory.h"
#include "Multigrid/TestMultigrid.h"
#include "Multigrid/Laplacian.h"

using namespace std;
using iVec = Vec<int,2>;
using rVec = Vec<Real,2>;

template <class T>
using Vector = std::vector<T>;

int main(int argc, char* argv[]){
  registerFactory<2>();
  TestMultigrid<2> TM1("../../test/Input6.json");
  unregisterFactory<2>();
  //D1Func1 f1;
  // MGFactory<1,Injection<1>,LinearInterpolation<1>>
  // MGF("../../test/Input2.json");
  // TestMultigrid<1,Injection<1>,LinearInterpolation<1>> TM(MGF);
  // TM.TestVCycle(f1,3);
}
