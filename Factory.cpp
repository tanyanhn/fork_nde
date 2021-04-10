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
    Info_Table* RK = new Info_Table;
    RK->set_type(4);
    return RK;
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

template<class CalPolicy>
Vector4d TimeIntegrator<CalPolicy>::one_step(Vector4d* u, const double dt, const double mu, int _acc){
  switch(ITable->get_type()){
  case 1: return AB_one_step(u,dt,mu,_acc,*ITable);
  case 2: return AM_one_step(u,dt,mu,_acc,*ITable);
  case 3: return BDF_one_step(u,dt,mu,_acc,*ITable);
  case 4: return RK_one_step(u[0],dt,mu);
  default: std::cerr << "No matching!"  << std::endl;
  }
}

template<class CalPolicy>
Vector4d TimeIntegrator<CalPolicy>::n_steps(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N){
  switch(ITable->get_type()){
  case 1: return AB_method(time,u0,dt,mu,_acc,N,*ITable);
  case 2: return AM_method(time,u0,dt,mu,_acc,N,*ITable);
  case 3: return BDF_method(time,u0,dt,mu,_acc,N,*ITable);
  case 4: return RK_method(time,u0,dt,mu,N);
  default: std::cerr << "No matching!"  << std::endl;
  }
}


