#include "Multigrid.h"

int main(int argc, char* argv[]){
  std::pair<double,double> boud(0,0);
  std::pair<int,double> crit(1,1e-8);
  int n = 256;
  double* init = new double[n-1];
  for (int i = 0 ; i < n-1 ; i++)
    init[i] = 0;
  Multigrid<full_weighting,linear,function2> V(boud,init,crit);
  double* f = V.righthand(n);
  double result = V.analysis_V_cycle(init,n,f,200,200);
  std::cout << "result:" << result << std::endl;
  /*
  std::cout << "result:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << result[i] << std::endl;
  double* sol = V.ref_solution(n);
  std::cout << "solution:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << sol[i] << std::endl;
  */
  return 0;
}
