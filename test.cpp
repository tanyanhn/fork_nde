#include "Info.h"
#include "cal.h"

int main( int argc, char* agrv[]){
  Info_Table AB;
  AB.load_data("Adams_Bashforth",1);
  /*std::cout << AB.get_n_table() << " " << AB.get_type() << std::endl;
  const Info& AB_info1 = AB.get_info(1);
  const Info& AB_info2 = AB.get_info(2);
  const Info& AB_info3 = AB.get_info(3);
  const Info& AB_info4 = AB.get_info(4);
  for (int i = 0; i < 3 ; i++)
    std::cout << AB_info.get_coe(i) << ",";
  std::cout << "\b " << std::endl;
  Info_Table AM;
  AM.load_data("Adams_Moulton",2);
  std::cout << AM.get_n_table() << " " << AM.get_type() << std::endl;
  const Info& AM_info = AM.get_info(7);
  for (int i = 0; i < 5 ; i++)
    std::cout << AM_info.get_coe(i) << ",";
  std::cout << "\b " << std::endl;
  Info_Table BDF;
  BDF.load_data("BDFs",3);
  std::cout << BDF.get_n_table() << " " << BDF.get_type() << std::endl;
  const Info& BDF_info = BDF.get_info(3);
  for (int i = 0; i < 4 ; i++)
    std::cout << BDF_info.get_coe(i) << ",";
  std::cout << "\b " << std::endl;*/
  Vector4d u0;
  double mu;
  u0 = initial_load("Initial2",&mu);
  double dt = 0.001;
  /*Vector4d v = AB_one_step(u0,dt,mu,1,AB);
  std::cout << v << std::endl;
  u[0] = v;
  v = AB_one_step(u,dt,mu,1,AB);
  std::cout << v << std::endl;
  u[0] = v;
  v = AB_one_step(u,dt,mu,1,AB);
  std::cout << v << std::endl;
  u[0] = v;
  v = AB_one_step(u,dt,mu,1,AB);
  std::cout << v << std::endl;*/
  Vector4d v = AB_method(u0,dt,mu,1,1,AB);
  std::cout <<  v << std::endl;
  Vector4d w = Newton(u0,mu,-dt,u0);
  std::cout << w << std::endl;
  Vector4d z = Newton(u0,mu,-0.5*dt,u0+0.5*dt*f(u0,mu));
  std::cout << z << std::endl;
  
  return 0;
}
