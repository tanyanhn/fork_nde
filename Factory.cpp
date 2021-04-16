#include "Factory.h"

template<class T>
TimeIntegratorFactory<T>::TimeIntegratorFactory(){
  pInstance_ = new T;
}

template<class T>
TimeIntegratorFactory<T>::~TimeIntegratorFactory(){};

template<class T>
T* TimeIntegratorFactory<T>::get_pointer(){
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

template<class CalPolicy,int acc>
Method<CalPolicy,acc>::Method(){
  ITable = CalPolicy().create_table();
}

template<class CalPolicy,int acc>
Method<CalPolicy,acc>::~Method(){}
 
template<class CalPolicy,int acc>
Info_Table* Method<CalPolicy,acc>::get_table(){
  return ITable;
}

template<class CalPolicy,int acc>
Vector4d Method<CalPolicy,acc>::one_step(Vector4d* u, const double dt, const double mu){
  switch(ITable->get_type()){
  case 1: return AB_one_step(u,dt,mu,acc,*ITable);
  case 2: return AM_one_step(u,dt,mu,acc,*ITable);
  case 3: return BDF_one_step(u,dt,mu,acc,*ITable);
  case 4: return RK_one_step(u[0],dt,mu);
  default: std::cerr << "No matching!"  << std::endl;
  }
}

template<class CalPolicy,int acc>
Vector4d Method<CalPolicy,acc>::n_steps(double& time, Vector4d u0, const double dt, const double mu, int N){
  switch(ITable->get_type()){
  case 1: return AB_method(time,u0,dt,mu,acc,N,*ITable);
  case 2: return AM_method(time,u0,dt,mu,acc,N,*ITable);
  case 3: return BDF_method(time,u0,dt,mu,acc,N,*ITable);
  case 4: return RK_method(time,u0,dt,mu,N);
  default: std::cerr << "No matching!"  << std::endl;
  }
}

template<class CalPolicy,int acc>
double Method<CalPolicy,acc>::err_Initial(double& time, Vector4d u0, const double dt, const double mu, const double T){
  switch(ITable->get_type()){
  case 1: case 2: case 3:
    return err_initial(time,u0,dt,mu,acc,T,*ITable,ITable->get_type());
  case 4:
    return err_initial(time,u0,dt,mu,T,4);
  default:
    std::cerr<< "No matching!" << std::endl;
  }
}

template<class CalPolicy,int acc>
double Method<CalPolicy,acc>::err_Richardson(double tol, double& time, Vector4d u0, double dt, const double mu, int N){
  switch(ITable->get_type()){
  case 1: case 2: case 3:
    return err_richardson(tol,time,u0,dt,mu,acc,N,*ITable,ITable->get_type());
  case 4:
    return err_richardson(tol,time,u0,dt,mu,N,4);
  default:
    std::cerr<< "No matching!" << std::endl;
  }
}

template<class CalPolicy,int acc>
double Method<CalPolicy,acc>::Grid_Refine1(Vector4d u0, double dt, const double mu, const double T){
  switch(ITable->get_type()){
  case 1: case 2: case 3:
    return grid_refine_err1(u0,dt,mu,acc,T,*ITable,ITable->get_type());
  case 4:
    return grid_refine_err1(u0,dt,mu,T,4);
  default:
    std::cerr<< "No matching!" << std::endl;
  }
}





template<class CalPolicy,int acc>
double Method<CalPolicy,acc>::Grid_Refine2(double tol,Vector4d u0, double dt, const double mu, int N){
  switch(ITable->get_type()){
  case 1: case 2: case 3:
    return grid_refine_err2(tol,u0,dt,mu,acc,N,*ITable,ITable->get_type());
  case 4:
    return grid_refine_err2(tol,u0,dt,mu,N,4);
  default:
    std::cerr<< "No matching!" << std::endl;
  }

      
}
