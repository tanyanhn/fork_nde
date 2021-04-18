#include "Factory.cpp"


int main(){
  std::cout <<"Test:"<<std::endl;
  /*Info_Table BDF;
  BDF.load_data("BDFs",3);
  //Info_Table AM;
  //AM.load_data("Adams_Moulton",2);
  Vector4d u0;
  double mu;
  double T;
  u0 = initial_load("Initial2",mu,T);
  double dt = 0.001;
  double time = 0;
  //double err = err_initial(time,u0,dt,mu,3,T,BDF,3);
  double err = grid_refine_err2(10e-12,u0,dt,mu,10000,4);
  std::cout << err << std::endl;*/


  std::string method = "Adams_Bashforth";
  int order = 4;
  std::string initial = "Initial2";
  //typedef std::pair<const std::string,int> fpair;
  TimeIntegrator *pIntegrator = TimeIntegratorFactory::Instance()->
    CreateIntegrator(std::make_pair(method,order));
  
  std::cout << pIntegrator->get_table()->get_type() << std::endl;
  Vector4d u0;
  double mu,T;
  u0 = initial_load(initial,mu,T);
  double dt = 0.001;
  double time;
  Vector4d v = pIntegrator->n_steps(time,u0,dt,mu,20000);
  std::cout << time << std::endl;
  std::cout << v << std::endl;
  return 0;
};


