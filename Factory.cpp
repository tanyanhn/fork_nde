#include "Factory.h"

template<class T>
T& TimeIntegratorFactory<T>::Instance(){
  if (!pInstance_)
    pInstance_ = new T;
  return pInstance_;
}












