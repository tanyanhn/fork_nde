#include "Factory.cpp"

int main(){
  std::cout <<"Test"<<std::endl;
  TimeIntegrator<BDFs> BDF;
  std::cout << (BDF.get_table())->get_type() << std::endl;
  Vector4d u0;
  double mu,T;
  u0 = initial_load("Initial1",mu,T);
  double dt = 0.0001;
  double time;
  double v = BDF.err_Initial(time,u0,dt,mu,4,T);
  std::cout << time << std::endl;
  std::cout << v << std::endl;
  return 0;
}
