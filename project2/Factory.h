#ifndef __MGGS_FACTORY__
#define __MGGS_FACTORY__

#include "Multigrid.h"
#include <map>
#include <array>

class MultigridFactory{
 public :
  typedef Grid* (*MultigridCreator)();
  static MultigridFactory* Instance();
  
 public:
  bool RegisterMultigrid(const std::array<std::string,3> GridId, MultigridCreator createGd);
  bool UnregisterMultigrid(const std::array<std::string,3> GridId);
  Grid* CreateMultigrid(const std::array<std::string,3> GridId) const;

 private:
  MultigridFactory();
  MultigridFactory(const MultigridFactory&);
  static MultigridFactory* pInstance_;

 private:
  typedef std::map<const std::array<std::string,3>, MultigridCreator> CallbackMap;
  CallbackMap callbacks_;
};


MultigridFactory* MultigridFactory::pInstance_ = 0;

MultigridFactory* MultigridFactory::Instance(){
  if (!pInstance_)
    pInstance_ = new MultigridFactory;
  return pInstance_;
}

MultigridFactory::MultigridFactory(){};
MultigridFactory::MultigridFactory(const MultigridFactory&){};


namespace{
  Grid* CreateGrid_hfl(){
    return new Multigrid<full_weighting,linear,function2>;
  }

  Grid* CreateGrid_hfq(){
    return new Multigrid<full_weighting,quadratic,function2>;
  }

  Grid* CreateGrid_hil(){
    return new Multigrid<injection,linear,function2>;
  }

  Grid* CreateGrid_hiq(){
    return new Multigrid<injection,quadratic,function2>;
  }

  Grid* CreateGrid_nfl(){
    return new Multigrid<full_weighting,linear,function1>;
  }
   
  Grid* CreateGrid_nfq(){
    return new Multigrid<full_weighting,quadratic,function1>;
  }

  Grid* CreateGrid_nil(){
    return new Multigrid<injection,linear,function1>;
  }

  Grid* CreateGrid_niq(){
    return new Multigrid<injection,quadratic,function1>;
  }


  
  std::array<std::string,3> hfl={"full_weighting","linear","homogeneous"};
  const bool _hfl_registered =
    MultigridFactory::Instance()->RegisterMultigrid(hfl,CreateGrid_hfl);

  std::array<std::string,3> hfq={"full_weighting","quadratic","homogeneous"};
  const bool _hfq_registered =
    MultigridFactory::Instance()->RegisterMultigrid(hfq,CreateGrid_hfq);

  std::array<std::string,3> hil={"injection","linear","homogeneous"};
  const bool _hil_registered =
    MultigridFactory::Instance()->RegisterMultigrid(hil,CreateGrid_hil);

  std::array<std::string,3> hiq={"injection","quadratic","homogeneous"};
  const bool _hiq_registered =
    MultigridFactory::Instance()->RegisterMultigrid(hiq,CreateGrid_hiq);

  std::array<std::string,3> nfl={"full_weighting","linear","nonhomogeneous"};
  const bool _nfl_registered =
    MultigridFactory::Instance()->RegisterMultigrid(nfl,CreateGrid_nfl);

  std::array<std::string,3> nfq={"full_weighting","quadratic","nonhomogeneous"};
  const bool _nfq_registered =
    MultigridFactory::Instance()->RegisterMultigrid(nfq,CreateGrid_nfq);

  std::array<std::string,3> nil={"injection","linear","nonhomogeneous"};
  const bool _nil_registered =
    MultigridFactory::Instance()->RegisterMultigrid(nil,CreateGrid_nil);

  std::array<std::string,3> niq={"injection","quadratic","nonhomogeneous"};
  const bool _niq_registered =
    MultigridFactory::Instance()->RegisterMultigrid(niq,CreateGrid_niq);
}


bool MultigridFactory::RegisterMultigrid(const std::array<std::string,3> GridId, MultigridCreator createGd){
  return callbacks_.insert(CallbackMap::value_type(GridId,createGd)).second;
}

bool MultigridFactory::UnregisterMultigrid(const std::array<std::string,3> GridId){
  return callbacks_.erase(GridId) == 1;
}


Grid* MultigridFactory::CreateMultigrid(const std::array<std::string,3> GridId) const{
  CallbackMap::const_iterator i = callbacks_.find(GridId);
  if (i == callbacks_.cend()){
    throw std::runtime_error("Unknown Grid ID");
  }
  return (i->second());
}




















#else
//do nothing
#endif
