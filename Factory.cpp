#include "Factory.h"

template<class T>
T& TimeIntegratorFactory<T>::Instance(){
  if (!pInstance_)
    pInstance_ = new T;
  return pInstance_;
}

class Adams_Bashforth{
  static Info_Table* create_table(){
    Info_Table* AB = new Info_Table;
    AB->load_data("Adams_Bashforth",1);
    return AB;
  }
};


template<class CalPolicy>
TimeIntegrator<CalPolicy>::TimeIntegrator(){
  ITable = CalPolicy().create_table();
}










