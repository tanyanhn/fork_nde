#include "Factory.cpp"


int main(){
  std::cout <<"Test:"<<std::endl;
  Info_Table AB;
  AB.load_data("Adams_Bashforth",1);
  Vector4d u0;
  double mu;
  double T;
  u0 = initial_load("Initial1",mu,T);
  double dt = 0.001;
  double time = 0;
  double err = err_richardson(10e-8,time,u0,dt,mu,4,10,AB,1);
  std::cout << err << std::endl;
  /*typedef TimeIntegratorFactory<TimeIntegrator<BDFs>> BDF_C;
  BDF_C A;
  std::cout << (A.get_pointer()->get_table())->get_type() << std::endl;
  Vector4d u0;
  double mu,T;
  u0 = initial_load("Initial1",mu,T);
  double dt = 0.00025;
  double time;
  double v = A.get_pointer()->err_Initial(time,u0,dt,mu,4,T);
  std::cout << time << std::endl;
  std::cout << v << std::endl;*/
  return 0;
}
