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
  u0 = initial_load("Initial1",&mu);
  double dt = 0.01;
  /*Vector4d v1 = AB_one_step(u,dt,mu,1,AB);
  std::cout << v1 << std::endl;
  u[1] = v1;
  Vector4d v2 = AB_one_step(u,dt,mu,2,AB);
  std::cout << v2 << std::endl;
  u[2] = v2;
  Vector4d v3 = AB_one_step(u,dt,mu,3,AB);
  u[3] = v3;
  Vector4d v4 = AB_one_step(u,dt,mu,4,AB);
  std::cout << v4 << std::endl;*/
  Vector4d v = AB_method(u0,dt,mu,1,3,AB);
  std::cout << v << std::endl;
  return 0;
}
