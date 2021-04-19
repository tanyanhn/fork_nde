#include "Factory.h"

TimeIntegrator::TimeIntegrator(){};
TimeIntegrator::~TimeIntegrator(){};


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
std::pair<double,double> Method<CalPolicy,acc>::Grid_Refine1(Vector4d u0, double dt, const double mu, const double T){
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

TimeIntegratorFactory* TimeIntegratorFactory::pInstance_ = 0;

TimeIntegratorFactory* TimeIntegratorFactory::Instance(){
  if (!pInstance_)
    pInstance_ = new TimeIntegratorFactory;
  return pInstance_;
}

TimeIntegratorFactory::TimeIntegratorFactory(){};
TimeIntegratorFactory::TimeIntegratorFactory(const TimeIntegratorFactory&){};

namespace{
  template <int acc>
  TimeIntegrator* CreateAB(){
    return new Method<Adams_Bashforth,acc>;
  }

  template <int acc>
  TimeIntegrator* CreateAM(){
    return new Method<Adams_Moulton,acc>;
  }

  template <int acc>
  TimeIntegrator* CreateBDF(){
    return new Method<BDFs,acc>;
  }

  TimeIntegrator* CreateRK(){
    return new Method<Runge_Kutta,0>;
  }

  const bool _AB1_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Bashforth",1),CreateAB<1>);
  const bool _AB2_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Bashforth",2),CreateAB<2>);
  const bool _AB3_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Bashforth",3),CreateAB<3>);
  const bool _AB4_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Bashforth",4),CreateAB<4>);

  const bool _AM2_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Moulton",2),CreateAM<2>);
  const bool _AM3_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Moulton",3),CreateAM<3>);
  const bool _AM4_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Moulton",4),CreateAM<4>);
  const bool _AM5_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Adams_Moulton",5),CreateAM<5>);

  const bool _BDF1_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("BDFs",1),CreateBDF<1>);
  const bool _BDF2_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("BDFs",2),CreateBDF<2>);
  const bool _BDF3_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("BDFs",3),CreateBDF<3>);
  const bool _BDF4_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("BDFs",4),CreateBDF<4>);

  const bool _RK_registered =
    TimeIntegratorFactory::Instance()->RegisterIntegrator(std::pair<const std::string,int>("Runge_Kutta",0),CreateRK);  
}

bool TimeIntegratorFactory::RegisterIntegrator(const std::pair<const std::string,int> &MethodId, IntegratorCreator createFn){
  return callbacks_.insert(CallbackMap::value_type(MethodId,createFn)).second;
}

bool TimeIntegratorFactory::UnregisterIntegrator(const std::pair<const std::string,int> &MethodId){
  return callbacks_.erase(MethodId) == 1;
}


TimeIntegrator* TimeIntegratorFactory::CreateIntegrator(const std::pair<const std::string,int> &MethodId) const{
  CallbackMap::const_iterator i = callbacks_.find(MethodId);
  if (i == callbacks_.cend()){
    throw std::runtime_error("Unknown Method ID");
  }
  return (i->second());
}
