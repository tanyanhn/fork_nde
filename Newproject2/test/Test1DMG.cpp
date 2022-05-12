#include "Multigrid/MGFactory.h"
#include "Multigrid/TestMultigrid.h"

using namespace std;
using iVec = Vec<int,1>;
using rVec = Vec<Real,1>;

template <class T>
using Vector = std::vector<T>;

int main(int argc, char* argv[]){
  string name = argv[1];
  string file = "../../test/"+name+".json";
  const int num = 3;
  cout << "V-cycle test:" << endl;
  cout << "--------------------------- Test function F = \\sin(\\pi x)---------------------------" << endl;
  ScalarFunction<1>* pD1Func1f = new D1Func1f;
  ScalarFunction<1>* pD1Func1F = new D1Func1F;
  ScalarFunction<1>* pD1Func1Fx = new D1Func1Fx;
  std::array<ScalarFunction<1>*,3> pfuncs{pD1Func1F,pD1Func1Fx};
  TestMultigrid<1> TM(file);
  TM.test(pD1Func1f,pfuncs,num);
}
