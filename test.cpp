#include "Info.h"

int main( int argc, char* agrv[]){
  Info_Table AB;
  AB.load_data("Adams_Bashforth",1);
  std::cout << AB.get_n_table() << " " << AB.get_type() << std::endl;
  const Info& AB_info = AB.get_info(3);
  for (int i = 0; i < 3 ; i++)
    std::cout << AB_info.get_coe(i) << ",";
  std::cout << "\b " << std::endl;
  /*Info_Table AM;
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
  return 0;
}
