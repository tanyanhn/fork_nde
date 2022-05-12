#include "Multigrid/MGFactory.h"
#include "Multigrid/TestMultigrid.h"

using namespace std;
using iVec = Vec<int,2>;
using rVec = Vec<Real,2>;

template <class T>
using Vector = std::vector<T>;

int main(int argc, char* argv[]){
  string name = argv[1];
  string file = "../../test/"+name+".json";
  const int num = 3;
  cout << "V-cycle test:" << endl;
  cout << "--------------------------- Test function F = \\exp(x*y) ---------------------------" << endl;
  ScalarFunction<2>* pD2Func1f = new D2Func1f;
  ScalarFunction<2>* pD2Func1F = new D2Func1F;
  ScalarFunction<2>* pD2Func1Fx = new D2Func1Fx;
  ScalarFunction<2>* pD2Func1Fy = new D2Func1Fy;
  std::array<ScalarFunction<2>*,3> pfuncs{pD2Func1F,pD2Func1Fx,pD2Func1Fy};
  TestMultigrid<2> TM(file);
  TM.test(pD2Func1f,pfuncs,num);
}
