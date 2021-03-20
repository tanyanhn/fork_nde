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

int Info::get_type() const{
  return type;
}

int Info::get_n_coe() const{
  return n_coe;
}

double Info::get_coe(int _idx) const{
  if (_idx >= n_coe){
    std::cerr << "error! coefficient no exist." << std::endl;
    return 0;
  }
  else
    return coe_info[_idx];
}


Info_Table::Info_Table(){
  table = NULL;
}

Info_Table::~Info_Table(){
  delete [] table;
}

void Info_Table::load_data(const char _file[], int _type){
  if (table != NULL){
    std::cerr << "Table is not empty!" << std::endl;
    return;
  }
  type = _type;
  std::fstream data(_file);
  data >> n_table;
  table = new Info [n_table];
  int k = 0;
  while (!data.eof()){
    int acc,n_coe;
    double* coe;
    data >> acc >> n_coe;
    coe = new double [n_coe];
    for (int i = 0; i < n_coe; i++)
      data >> coe[i];
    table[k].set(acc,type,coe);
    k++;
    delete [] coe;
  }
  data.close();
}

int Info_Table::get_n_table() const{
  return n_table;
}

int Info_Table::get_type() const{
  return type;
}

const Info& Info_Table::get_info(int _acc) const{
  int i = 0;
  for (; i < n_table ; i++){
    if (table[i].get_acc() == _acc)
      break;
  }
  if (table[i].get_acc() == _acc)
    return table[i];
  else{
    std::cout << "Not Found!" << std::endl;
    return table[n_table-1];
  }
}


