#include "Factory.h"

int main(int argc, char* argv[]){
  std::string Restriction = "full_weighting";
  std::string Interpolation = "linear";
  std::string boundary = "nonhomogeneous";
  std::array<std::string,3> input = {Restriction,Interpolation,boundary};
  Grid *pGrid = MultigridFactory::Instance()->
    CreateMultigrid(input);
  std::pair<double,double> boud(1,exp(sin(1)));
  std::pair<int,double> crit(1,1e-8);
  int n = 256;
  double* init = new double[n-1];
  for (int i = 0 ; i < n-1 ; i++)
    init[i] = 0;
  pGrid->load_boundary(boud);
  pGrid->load_initial(init);
  pGrid->load_criteria(crit);
  double* f = pGrid->righthand(n);
  double result = pGrid->analysis_V_cycle(init,n,f,200,200);
  std::cout << "result:" << result << std::endl;
  /*
  std::cout << "result:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << result[i] << std::endl;
  double* sol = pGrid->ref_solution(n);
  std::cout << "solution:" << std::endl;
  for (int i = 0 ; i < n - 1 ; i++)
    std::cout << sol[i] << std::endl;
  */
  return 0;
}
