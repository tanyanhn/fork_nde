#include "Info.h"
#include <eigen/Eigen/Dense>
#include <cmath>

using namespace Eigen;

Vector4d f(Vector4d u, double mu);
Vector4d initial_load(const char _file[],double* _mu);
Vector4d AB_one_step(Vector4d* u, const double dt, double mu, int _acc, const Info& info);
Vector4d AB_one_step(Vector4d* u, const double dt, double mu, int _acc, const Info_Table& table);
Vector4d AB_method(Vector4d u0, const double dt, double mu, int _acc, int N, const Info_Table& table);


Vector4d f(Vector4d u, double mu){
  double v0,v1,v2,v3;
  double tmp1 = u(0) + mu - 1;
  double tmp2 = u(0) + mu;
  double den1 = pow(u(1) * u(1) + tmp1 * tmp1, 3/2);
  double den2 = pow(u(1) * u(1) + tmp2 * tmp2, 3/2);
  v0 = u(2);
  v1 = u(3);
  v2 = 2 * u(3) + u(0) - mu * tmp1 / den1 - (1 - mu) * tmp2 / den2;
  v3 = -2 * u(2) + u(1) - mu * u(1) / den1 - (1 - mu) * u(1) / den2;
  Vector4d v(v0,v1,v2,v3);
  return v;
}

Vector4d initial_load(const char _file[],double* _mu){
  std::fstream data(_file);
  data >>*_mu;
  Vector4d u;
  for (int i = 0 ; i < 4 ; i++)
    data >> u(i);
  data.close();
  return u;
}

Vector4d AB_one_step(Vector4d* u, const double dt, double mu, int _acc, const Info& info){
  Vector4d v=u[_acc-1];
  for (int i = 1; i <= _acc ; i++)
    v += f(u[_acc-i],mu)*info.get_coe(i-1)*dt;
  return v;
}

Vector4d AB_one_step(Vector4d* u, const double dt, double mu, int _acc, const Info_Table& table){
  const Info& AB_info=table.get_info(_acc);
  return AB_one_step(u,dt,mu,_acc,AB_info);
}

Vector4d AB_method(Vector4d u0, const double dt, double mu, int _acc, int N, const Info_Table& table){
  Vector4d* u;
  u = new Vector4d [_acc];
  u[0] = u0;
  int k = 1;
  while (k < _acc){
    Vector4d v = AB_one_step(u,dt,mu,k,table);
    u[k] = v;
    k++;
  }
  std::ofstream os;
  os.open("AB_")
  const Info& AB_info = table.get_info(_acc);
  for (int i = 0; i < N - _acc + 1; i++){
    Vector4d v = AB_one_step(u,dt,mu,_acc,AB_info);
    for (int j = 0; j < _acc - 1 ; j++)
      u[j] = u[j+1];
    u[_acc-1] = v;
  }
  return u[_acc-1];
}
