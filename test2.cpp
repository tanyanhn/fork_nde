#include "Factory.cpp"

int main(){
  std::cout <<"Test"<<std::endl;
  TimeIntegrator<Adams_Bashforth> AB;
  std::cout << (AB.get_table())->get_n_table() << std::endl;
  return 0;
}
