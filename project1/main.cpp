#include "Factory.h"
#include <iostream>

#define tol 1e-8

int main(int argc, char* argv[]){
  std::string garbage,Method,Initial,Inputdata;
  double dt,mu,T,time,result;
  std::pair<double,double> results;
  Vector4d u0,v;
  int Index,Order,N;
  int err_analysis,grid_refine_analysis;
  std::fstream Input(argv[1]);
  while (getline(Input,garbage,'\n')){
    Input >> Index >> Method >> Order
	 >> dt >> Initial  >> N
	 >>  err_analysis >>  grid_refine_analysis;
    TimeIntegrator *pIntegrator = TimeIntegratorFactory::Instance()->
      CreateIntegrator(std::make_pair(Method,Order));
    u0 = initial_load(Initial,mu,T);
    if ( Initial == "Initial1"){
      if ( grid_refine_analysis ){
	results = pIntegrator->Grid_Refine1(u0,dt,mu,T);
	std::cout << "Problem " << Index << ": "
	          << "Order = " << results.first
	          << "  Error constant = " << results.second << std::endl
		  << "Analysis can be saw by runing analysis.m file" << std::endl; 
      }
      else{
	if ( err_analysis ){
	  result = pIntegrator->err_Initial(time,u0,dt,mu,T);
	  std::cout << "Problem " << Index << ": " << Method << " , " << "order = " << Order << std::endl
		    << "err: " << result << ",CPU time: " << time << "(ms)" << std::endl
		    << "        orbit can be saw by runing .m file" << std::endl;
	}
	else{
	  v = pIntegrator->n_steps(time,u0,dt,mu,N);
	  std::cout << "Problem " << Index << ": "
		    << "relust: (" << v(0) << "," << v(1) << ",0)"
		    <<",CPU time: " << time << "(ms)" << std::endl
		    << "        orbit can be saw by runing .m file" << std::endl;
	}
      }
    }
    else if ( Initial == "Initial2"){
      if ( grid_refine_analysis ){
	result = pIntegrator->Grid_Refine2(tol,u0,dt,mu,N);
	std::cout << "Problem " << Index << ": "
		  << "Results and analysis can be saw by runing analysis.m file" << std::endl; 
      }
      else{
	if ( err_analysis ){
	  result = pIntegrator->err_Richardson(tol,time,u0,dt,mu,N);
	  std::cout << "Problem " << Index << ": " << Method << " , " << "order = " << Order << std::endl
		    << "err: " << result << ",CPU time: " << time << "(ms)" << std::endl
		    << "        orbit can be saw by runing .m file" << std::endl;
	}
	else{
	  v = pIntegrator->n_steps(time,u0,dt,mu,N);
	  std::cout << "Problem " << Index << ": "
		    << "relust: (" << v(0) << "," << v(1) << ",0)"
		    <<",CPU time: " << time << "(ms)" << std::endl
		    << "        orbit can be saw by runing .m file" << std::endl;
	}
      }
    }
  }
  Input.close();
  std::cout << "All done!" << std::endl;

  return 0;
}
