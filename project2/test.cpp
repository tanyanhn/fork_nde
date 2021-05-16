#include "Multigrid.h"

int main(int argc, char* argv[]){
  std::pair<double,double> boud(1,exp(sin(1)));
  std::pair<int,double> crit(1,0.00000001);
  int n = 1024;
  double* init = new double[n-1];
  for (int i = 0 ; i < n-1 ; i++)
    init[i] = i;
  Multigrid<full_weighting,linear> V(boud,init,crit);
  double* f = V.righthand(n);
  double* result = V.fm_cycle(n,f,2000,2000);
  std::cout << "result:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << result[i] << std::endl;
  double* sol = V.ref_solution(n);
  std::cout << "solution:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << sol[i] << std::endl;
  return 0;
}
