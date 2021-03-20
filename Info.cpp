#include "Info.h"

Info::Info(){
  type = 0;
}

Info::~Info(){
  delete [] coe_info;
}

void Info::set(int _acc, int _type, double* _coe_info){
  accuracy = _acc;
  type = _type;
  if (type  == 1 || type == 2)
    n_coe = _acc;
  else if (type == 3)
    n_coe = _acc + 1;
  coe_info = new double [n_coe];
  for (int i = 0 ; i < n_coe ; i++ )
    coe_info[i] = _coe_info[i];
}

int Info::get_acc() const{
  return accuracy;
}

int Info::get_n_coe() const{
  return n_coe;
}

double Info::get_coe(int _idx) const{
  if (_idx >= n_coe)
    std::cerr << "error! coefficient no exist." << std::endl;
  else
    return coe_info[_idx];
}
