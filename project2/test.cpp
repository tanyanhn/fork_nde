#include "Multigrid.h"

int main(int argc, char* argv[]){
  std::pair<double,double> boud(0,0);
  std::pair<int,double> crit(0,1e-8);
  int n = 8;
  double* init = new double[n-1];
  for (int i = 0 ; i < n-1 ; i++)
    init[i] = i;
  Multigrid<injection,quadratic> V(boud,init,crit);
  /*V.test_restriction(n);
  V.test_interpolation(n);
  double* p = V.get_initial();
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << p[i] << std::endl;*/
  double* A = V.lefthand(n);
  double* f = V.righthand(n);
  double weight = 1;
  double* result = V.weighted_Jacobi(A,f,init,n,weight,1);
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << result[i] << std::endl;
  return 0;
}
