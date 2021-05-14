#include "Multigrid.h"

int main(int argc, char* argv[]){
  std::pair<double,double> boud(-2.3,0.3);
  std::pair<int,double> crit(0,1e-8);
  int n = 8;
  double* init = new double[n-1];
  init[0]=1;init[1]=3;init[2]=6;
  init[3]=1.4;init[4]=0;init[5]=-2;init[6]=2.4;
  Multigrid<injection,quadratic> V(boud,init,crit);
  V.test_restriction(n);
  V.test_interpolation(n);
  double* p = V.get_initial();
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << p[i] << std::endl;
  return 0;
}
