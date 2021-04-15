#include "Factory.cpp"


int main(){
  std::cout <<"Test:"<<std::endl;
  Info_Table BDF;
  BDF.load_data("BDFs",3);
  //Info_Table AM;
  //AM.load_data("Adams_Moulton",2);
  Vector4d u0;
  double mu;
  double T;
  u0 = initial_load("Initial1",mu,T);
  double dt = 0.0002;
  double time = 0;
  //double err = err_initial(time,u0,dt,mu,3,T,BDF,3);
  double cr = grid_refine_err1(u0,dt,mu,T,4);
  std::cout << cr << std::endl;
  
  /* typedef TimeIntegratorFactory<TimeIntegrator<Runge_Kutta>> RK_C;
  RK_C A;
  std::cout << (A.get_pointer()->get_table())->get_type() << std::endl;
  Vector4d u0;
  double mu,T;
  u0 = initial_load("Initial2",mu,T);
  double dt = 0.01;
  double time;
  double v = A.get_pointer()->err_Richardson(10e-8,time,u0,dt,mu,0,2000);
  std::cout << time << std::endl;
  std::cout << v << std::endl;*/
  return 0;
}
