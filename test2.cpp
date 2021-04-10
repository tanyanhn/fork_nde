#include "Factory.cpp"

int main(){
  std::cout <<"Test"<<std::endl;
  TimeIntegrator<Runge_Kutta> RK;
  std::cout << (RK.get_table())->get_type() << std::endl;
  Vector4d u0;
  double mu;
  u0 = initial_load("Initial1",mu);
  double dt = 0.001;
  double time;
  Vector4d v = RK.n_steps(time,u0,dt,mu,4,20000);
  std::cout << v << std::endl;
  return 0;
}
