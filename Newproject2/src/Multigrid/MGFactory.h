#ifndef _MGFACTORY_H_
#define _MGFACTORY_H_

#include "Multigrid/MultigridSolver.h"
#include <memory>
#include <map>


template <int Dim> class MGFactory;

template <int Dim>
class MGFactory{
public:
  using CreateMultigridSolverCallback = std::unique_ptr<MultigridSolver<Dim> >
    (*)(RectDomain<Dim>, const char*);
private:
  using CallbackMap = std::map<std::string, CreateMultigridSolverCallback>;
public:
  static MGFactory& createFactory(){
    static MGFactory object;
    return object;
  }

  bool registerMultigridSolver(std::string multigridId, CreateMultigridSolverCallback createFn){
    return callbacks_.insert(typename CallbackMap::value_type(multigridId, createFn)).
      second;
  }

  bool unregisterMultigridSolver(std::string multigridId){
    return callbacks_.erase(multigridId) == 1;
  }

  template<class ...TS>
  std::unique_ptr<MultigridSolver<Dim>> createMultigridSolver(std::string
                                                  multigridId, TS && ...args){
    auto it = callbacks_.find(multigridId);
    if (it == callbacks_.end())
      {
        throw std::runtime_error("Unknown MultigridSolver ID. ");
      }
    return (it->second)(std::forward<TS>(args)...);
  }

private:
  MGFactory() = default;
  MGFactory(const MGFactory&) = default;
  MGFactory& operator=(const MGFactory&) = default;
  ~MGFactory() = default;
  
private:
  CallbackMap callbacks_;
  
};

template <int Dim>
void registerFactory(){
  MGFactory<Dim>& object = MGFactory<Dim>::createFactory();
  object.registerMultigridSolver("Injection+LinearInterpolation",
                           createMultigridSolver1);
  object.registerMultigridSolver("Injection+QuadraticInterpolation",
                           createMultigridSolver2);
  object.registerMultigridSolver("FullWeightingRestriction+LinearInterpolation",
                           createMultigridSolver3);
  object.registerMultigridSolver("FullWeightingRestriction+QuadraticInterpolation",
                           createMultigridSolver4);
  return;
}

template <int Dim>
void unregisterFactory(){
  MGFactory<Dim>& object = MGFactory<Dim>::createFactory();
  object.unregisterMultigridSolver("Injection+LinearInterpolation");
  object.unregisterMultigridSolver("Injection+QuadraticInterpolation");
  object.unregisterMultigridSolver("FullWeightingRestriction+LinearInterpolation");
  object.unregisterMultigridSolver("FullWeightingRestriction+QuadraticInterpolation");
  return;
}

#endif //  _MGFACTORY_H_
