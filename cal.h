#include "Info.h"
#include <eigen/Eigen/Dense>
#include <cmath>

using namespace Eigen;

Vector4d f(Vector4d u, double mu);
Vector4d initial_load(const char _file[],double* _mu);

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
