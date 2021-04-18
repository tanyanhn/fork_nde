#include "Factory.cpp"

void getInput(std::fstream &_Input, std::string &_Method, int &_Order, double &_dt, std::string &_Initial, int &_N, bool &_err_analysis, bool &_grid_refine_analysis){
  _Input >> _Method >> _Order >> _dt
	 >> _Initial  >> _N
	 >> _err_analysis >> _grid_refine_analysis;
}


int main(int argc, char* argv[]){
  std::string garbage,Method,Initial;
  double dt,T;
  Vector4d u0;
  int Order,N;
  bool err_analysis,grid_refine_analysis;
  std::fstream Input("Input");
  getline(Input,garbage);
  while (getline(Input,garbage,'\n')){
    getInput(Input,Method,Order,dt,Initial,N,err_analysis,grid_refine_analysis);
  }
  std::cout << Method << std::endl << Order << std::endl << dt <<std::endl << Initial <<std::endl << N <<std::endl << err_analysis <<std::endl << grid_refine_analysis << std::endl;

  return 0;
}
