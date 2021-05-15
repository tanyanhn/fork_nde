#include "Multigrid.h"

int main(int argc, char* argv[]){
  std::pair<double,double> boud(0,0);
  std::pair<int,double> crit(0,1e-8);
  int n = 8;
  double* init = new double[n-1];
  init[0]=1;init[1]=3;init[2]=6;
  init[3]=1.4;init[4]=0;init[5]=-2;init[6]=2.4;
  Multigrid<injection,quadratic> V(boud,init,crit);
  /*V.test_restriction(n);
  V.test_interpolation(n);
  double* p = V.get_initial();
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << p[i] << std::endl;*/
  double* A = V.lefthand(n);
  double* f = V.righthand(n);
  std::cout << "left:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++){
    for (int j = 0 ; j < n - 1 ; j++)
      std::cout << A[i*(n-1)+j] << " ";
    std::cout << std::endl;
  }
  std::cout << "right:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << f[i] << std::endl;
  return 0;
}
