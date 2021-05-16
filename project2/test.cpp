#include "Multigrid.h"

int main(int argc, char* argv[]){
  std::pair<double,double> boud(1,exp(sin(1)));
  std::pair<int,double> crit(0,1e-8);
  int n = 64;
  double* init = new double[n-1];
  for (int i = 0 ; i < n-1 ; i++)
    init[i] = i;
  Multigrid<injection,quadratic> V(boud,init,crit);
  double* f = V.righthand(n);
  double* result = V.onestep_V_cycle(init,n,f,2000,2000);
  std::cout << "result:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << result[i] << std::endl;
  double* sol = V.ref_solution1(n);
  std::cout << "solution:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << sol[i] << std::endl;
  return 0;
}
