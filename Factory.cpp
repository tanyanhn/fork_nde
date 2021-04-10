#include "Factory.h"

template<class T>
T& TimeIntegratorFactory<T>::Instance(){
  if (!pInstance_)
    pInstance_ = new T;
  return pInstance_;
}

class Adams_Bashforth{
 public:
  static Info_Table* create_table(){
    Info_Table* AB = new Info_Table;
    AB->load_data("Adams_Bashforth",1);
    return AB;
  }
};

class Adams_Moulton{
 public:
   static Info_Table* create_table(){
     Info_Table* AM = new Info_Table;
     AM->load_data("Adams_Moulton",2);
     return AM;
  } 
};

class BDFs{
 public:
   static Info_Table* create_table(){
     Info_Table* BDF = new Info_Table;
     BDF->load_data("BDFs",3);
     return BDF;
  } 
};

class Runge_Kutta{
 public:
  static Info_Table* create_table(){
    return NULL;
  }
};

template<class CalPolicy>
TimeIntegrator<CalPolicy>::TimeIntegrator(){
  ITable = CalPolicy().create_table();
}


template<class CalPolicy>
TimeIntegrator<CalPolicy>::~TimeIntegrator(){}
 
template<class CalPolicy>
Info_Table* TimeIntegrator<CalPolicy>::get_table(){
  return ITable;
}




