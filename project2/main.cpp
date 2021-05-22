#include "Factory.h"

int main(int argc, char* argv[]){
  std::string restriction,interpolation,boundary,cycles,garbage;
  std::pair<int,double> crit;
  int initial_type,criteria_type,analysis,index,n;
  double criteria;
  std::fstream Input(argv[1]);
  while (getline(Input,garbage,'\n')){
    Input >> index >> boundary >> restriction >> interpolation
	  >> cycles >> criteria_type >> criteria >> initial_type
	 >> n >> analysis;
    std::array<std::string,3> inputstring = {restriction,interpolation,boundary};
    Grid *pGrid = MultigridFactory::Instance()->
    CreateMultigrid(inputstring);
    std::pair<int,double> crit(criteria_type,criteria);
    std::pair<double,double> boud;
    double* init = new double[n-1];
    if (boundary == "nonhomogeneous"){
      boud.first = 1; boud.second = exp(sin(1));
    }
    else if (boundary == "homogeneous"){
      boud.first = 0; boud.second = 0;
    }
    else{
      std::cerr << "Wrong boundary condition!" << std::endl;
      exit(-1);
    }
    if (initial_type == 0){
      for (int i = 0 ; i < n-1 ; i++)
	init[i] = 0;
    }
    else{
      std::cerr << "Wrong initial condition!" << std::endl;
      exit(-1);
    }
    pGrid->load_boundary(boud);
    pGrid->load_initial(init);
    pGrid->load_criteria(crit);
    double* f = pGrid->righthand(n);
    if (cycles == "V_cycle"){
      if (analysis){
	std::cout << "Problem " << index << ": " << std::endl; 
	double result = pGrid->analysis_V_cycle(init,n,f,200,200);
	std::cout << "Max_norm of error vector is " << result << std::endl;
	std::cout << "Run V_cycle_" << n << ".m by matlab to get the reduction rate of the residuals." << std::endl;
      }
      else{
	std::cout << "Problem " << index << ": " << std::endl;
	double* result = pGrid->V_cycle(init,n,f,200,200);
	pGrid->show_result(result,n,cycles);
	std::cout << "Run V_cycle_show_" << n << ".m by matlab to get calculation effect of result." << std::endl;
      }
    }
    else if (cycles == "full_multigrid"){
      if (analysis){
	std::cout << "Problem " << index << ": " << std::endl;
	double result = pGrid->analysis_fm_cycle(n,f,2000,2000);
	std::cout << "Max_norm of error vector is " << result << std::endl;
      }
      else{
	std::cout << "Problem " << index << ": " << std::endl;
	double* result = pGrid->fm_cycle(n,f,2000,2000);
	pGrid->show_result(result,n,cycles);
	std::cout << "Run full_multigrid_show_" << n << ".m by matlab to get calculation effect of result." << std::endl;
      }
    }
  }
  Input.close();
  std::cout << "All done!" << std::endl;
  return 0;
}
