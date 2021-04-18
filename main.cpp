#include "Factory.cpp"
#include <iostream>

#define tol 1e-8

void getInput(std::fstream &_Input, int &_Index, std::string &_Method,
	      int &_Order, double &_dt, std::string &_Initial, int &_N,
	      int &_err_analysis, int &_grid_refine_analysis){
  _Input >> _Index >> _Method >> _Order
	 >> _dt >> _Initial  >> _N
	 >>  _err_analysis >>  _grid_refine_analysis;
}


int main(int argc, char* argv[]){
  std::string garbage,Method,Initial;
  double dt,mu,T,time,result;
  Vector4d u0,v;
  int Index,Order,N;
  int err_analysis,grid_refine_analysis;
  std::fstream Input("Input");
  while (getline(Input,garbage,'\n')){
    getInput(Input,Index,Method,Order,dt,Initial,N,err_analysis,grid_refine_analysis);
    TimeIntegrator *pIntegrator = TimeIntegratorFactory::Instance()->
      CreateIntegrator(std::make_pair(Method,Order));
    u0 = initial_load(Initial,mu,T);
    if ( Initial == "Initial1"){
      if ( grid_refine_analysis ){
	result = pIntegrator->Grid_Refine1(u0,dt,mu,T);
	std::cout << "Problem " << Index << ": "
	  //
		  << "Analysis can be saw by runing .m file" << std::endl; 
      }
      else{
	if ( err_analysis ){
	  result = pIntegrator->err_Initial(time,u0,dt,mu,T);
	  std::cout << "Problem " << Index << ": "
		    << "err: " << result << ",CPU time: " << time << "(ms)" << std::endl
		    << "orbit can be saw by runing .m file" << std::endl;
	}
	else{
	  v = pIntegrator->n_steps(time,u0,dt,mu,N);
	  std::cout << "Problem " << Index << ": "
		    << "relust: (" << v(0) << "," << v(1) << ",0)"
		    <<",CPU time: " << time << "(ms)" << std::endl
		    << "orbit can be saw by runing .m file" << std::endl;
	}
      }
    }
    else if ( Initial == "Initial2"){
      if ( grid_refine_analysis ){
	result = pIntegrator->Grid_Refine2(tol,u0,dt,mu,N);
	std::cout << "Problem " << Index << ": "
	  //
		  << "Analysis can be saw by runing .m file" << std::endl; 
      }
      else{
	if ( err_analysis ){
	  result = pIntegrator->err_Richardson(tol,time,u0,dt,mu,N);
	  std::cout << "Problem " << Index << ": "
		    << "err: " << result << ",CPU time: " << time << "(ms)" << std::endl
		    << "orbit can be saw by runing .m file" << std::endl;
	}
	else{
	  v = pIntegrator->n_steps(time,u0,dt,mu,N);
	  std::cout << "Problem " << Index << ": "
		    << "relust: (" << v(0) << "," << v(1) << ",0)"
		    <<",CPU time: " << time << "(ms)" << std::endl
		    << "orbit can be saw by runing .m file" << std::endl;
	}
      }
    }
  }
  Input.close();
  std::cout << "All done!" << std::endl;

  return 0;
}
