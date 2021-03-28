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
  std::cout << "\b " << std::endl;*/
  Info_Table AM;
  AM.load_data("Adams_Moulton",2);
  /*std::cout << AM.get_n_table() << " " << AM.get_type() << std::endl;
  const Info& AM_info = AM.get_info(7);
  for (int i = 0; i < 5 ; i++)
    std::cout << AM_info.get_coe(i) << ",";
  std::cout << "\b " << std::endl;*/
  Info_Table BDF;
  BDF.load_data("BDFs",3);
  /*std::cout << BDF.get_n_table() << " " << BDF.get_type() << std::endl;
  const Info& BDF_info = BDF.get_info(3);
  for (int i = 0; i < 4 ; i++)
    std::cout << BDF_info.get_coe(i) << ",";
  std::cout << "\b " << std::endl;*/
  Vector4d u0;
  double mu;
  double T;
  u0 = initial_load("Initial1",&mu,&T);
  // std::cout << u0 << std::endl << mu << "," << T << std::endl;
  double dt = 0.00005;
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
  double time = 0;
  /*
  Vector4d v = AB_method(&time,u0,dt,mu,4,20000,AB);
  std::cout << "AB: " << v << std::endl;
  std::cout << "AB CPU time:" << time << std::endl;
  Vector4d w = AM_method(&time,u0,dt,mu,5,20000,AM);
  std::cout << "AM: " << w << std::endl;
  std::cout << "AM CPU time:" << time << std::endl;
  Vector4d z = BDF_method(&time,u0,dt,mu,4,20000,BDF);
  std::cout << "BDF: " << z << std::endl;
  std::cout << "BDF CPU time:" << time << std::endl;
  Vector4d r = RK_method(&time,u0,dt,mu,20000);
  std::cout << "RK: "<< r << std::endl;
  std::cout << "RK CPU time:" << time << std::endl;
  */
  //double err1 = err_initial(&time,u0,dt,mu,5,T,AM,2);
  //std::cout << "AM CPU time : " << time << std::endl;
  //std::cout << "AM err : " << err1 << std::endl;
  //double err2 = err_initial(&time,u0,dt,mu,T,4);
  //std::cout << "RK CPU time : " << time << std::endl;
  //std::cout << "RK err : " << err2 << std::endl;
  
  return 0;
}
