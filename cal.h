#ifndef __MGGS_CAL__
#define __MGGS_CAL__
#include "Info.h"
#include <eigen/Eigen/Dense>
#include <cmath>
#include <ctime>
#include <string>

using namespace Eigen;

Vector4d f(Vector4d u, const double mu);
Vector4d F(Vector4d u, const double mu, const double k, Vector4d b);
Matrix4d DF(Vector4d u, const double mu, const double k);
Vector4d Newton(Vector4d u0, const double mu, const double k, Vector4d b);

Vector4d initial_load(const char _file[], double& _mu);
Vector4d initial_load(const char _file[], double& _mu, double& _T);

Vector4d AB_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info& info);
Vector4d AB_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info_Table& table);
Vector4d AB_method(Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table);
Vector4d AB_method(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table);

Vector4d AM_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info& info);
Vector4d AM_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info_Table& table);
Vector4d AM_method(Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table);
Vector4d AM_method(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table);

Vector4d BDF_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info& info);
Vector4d BDF_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info_Table& table);
Vector4d BDF_method(Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table);
Vector4d BDF_method(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table);

Vector4d RK_one_step(Vector4d u, const double dt, const double mu);
Vector4d RK_method(Vector4d u0, const double dt, const double mu, int N);
Vector4d RK_method(double& time, Vector4d u0, const double dt, const double mu, int N);
  
double err_initial(double& time, Vector4d u0, const double dt, const double mu, const int _acc, double T, const Info_Table& table, int type);
double err_initial(double& time, Vector4d u0, const double dt, const double mu, double T, int type);

Vector4d f(Vector4d u, const double mu){
  double v0,v1,v2,v3;
  double tmp1 = u(0) + mu - 1;
  double tmp2 = u(0) + mu;
  double den1 = pow(u(1) * u(1) + tmp1 * tmp1, 1.5);
  double den2 = pow(u(1) * u(1) + tmp2 * tmp2, 1.5);
  v0 = u(2);
  v1 = u(3);
  v2 = 2 * u(3) + u(0) - mu * tmp1 / den1 - (1 - mu) * tmp2 / den2;
  v3 = -2 * u(2) + u(1) - mu * u(1) / den1 - (1 - mu) * u(1) / den2;
  Vector4d v(v0,v1,v2,v3);
  return v;
}

Vector4d F(Vector4d u, const double mu, const double k, Vector4d b){
  Vector4d v = u + k * f(u,mu) - b;
  return v;
}

Matrix4d DF(Vector4d u, const double mu, const double k){
  Matrix4d D;
  double tmp1 = u(0) + mu - 1;
  double tmp2 = u(0) + mu;
  double den1 = pow(u(1) * u(1) + tmp1 * tmp1, 2.5);
  double den2 = pow(u(1) * u(1) + tmp2 * tmp2, 2.5);
  double result1 = 1 - mu * (u(1) * u(1) - 2 * tmp1 * tmp1) / den1 -
    (1 - mu) * (u(1) * u(1) - 2 * tmp2 * tmp2) / den2;
  double result2 = 3 * mu * u(1) * tmp1 / den1 + 3 * (1 - mu) * u(1) * tmp2 / den2;
  double result3 = 1 - mu * (-2 * u(1) * u(1) + tmp1 * tmp1) / den1 -
    (1 - mu) * (-2 * u(1) * u(1) + tmp2 * tmp2) /den2;
  D << 1,0,k,0,
    0,1,0,k,
    k*result1,k*result2,1,2*k,
    k*result2,k*result3,-2*k,1;
  return D;
}

Vector4d Newton(Vector4d u0, const double mu, const double k, Vector4d b){
  Vector4d v = u0;
  Matrix4d M = DF(v,mu,k);
  Vector4d r = -F(v,mu,k,b);
  Vector4d dv;
  dv = M.fullPivHouseholderQr().solve(r);
  v += dv;
  r = -F(v,mu,k,b);
  while( r.norm() > 10e-16){
    M = DF(v,mu,k);
    dv = M.fullPivHouseholderQr().solve(r);
    v += dv;
    r = -F(v,mu,k,b);
  }
  return v;
}


Vector4d initial_load(const char _file[],double& _mu){
  std::fstream data(_file);
  data >>_mu;
  Vector4d u;
  for (int i = 0 ; i < 4 ; i++)
    data >> u(i);
  data.close();
  return u;
}

Vector4d initial_load(const char _file[], double& _mu, double& _T){
  std::fstream data(_file);
  data >>_mu;
  Vector4d u;
  for (int i = 0 ; i < 4 ; i++)
    data >> u(i);
  data >>_T;
  data.close();
  return u;
}


Vector4d AB_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info& info){
  Vector4d v=u[_acc-1];
  for (int i = 1; i <= _acc ; i++)
    v += f(u[_acc-i],mu)*info.get_coe(i-1)*dt;
  return v;
}

Vector4d AB_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info_Table& table){
  const Info& AB_info=table.get_info(_acc);
  return AB_one_step(u,dt,mu,_acc,AB_info);
}

Vector4d AB_method(Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table){
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
  std::string ss = "AB_" + std::to_string(_acc) + "_" + std::to_string(N) + ".m";
  const char *s = ss.c_str();
  os.open(s);
  os << "a=[\n";
  for (int i = 0; i < _acc ; i++)
    os << u[i](0) << "," << u[i](1) << ";\n";
  const Info& AB_info = table.get_info(_acc);
  for (int i = 0; i < N - _acc + 1; i++){
    Vector4d v = AB_one_step(u,dt,mu,_acc,AB_info);
    for (int j = 0; j < _acc - 1 ; j++)
      u[j] = u[j+1];
    u[_acc-1] = v;
    os << v(0) << "," << v(1) << ";\n";
  }
  os << "];\n";
  os << "x = a(:,1);";
  os << "y = a(:,2);";
  os << "plot(x,y)";
  os.close();
  return u[_acc-1];
}

Vector4d AB_method(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table){
  clock_t t1 = clock();
  Vector4d v = AB_method(u0,dt,mu,_acc,N,table);
  clock_t t2 = clock();
  time = (double)(t2 - t1) / CLOCKS_PER_SEC * 1000;
  return v;
}

Vector4d AM_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info& info){
  Vector4d b = u[_acc-2];
  for ( int i = 2; i <= _acc; i++)
    b += f(u[_acc-i],mu)*info.get_coe(i-1)*dt;
  Vector4d v = Newton(u[_acc-2],mu,-dt*info.get_coe(0),b);
  return v;
}

Vector4d AM_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info_Table& table){
  const Info& AM_info=table.get_info(_acc);
  return AM_one_step(u,dt,mu,_acc,AM_info);
}

Vector4d AM_method(Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table){
  Vector4d* u;
  u = new Vector4d [_acc-1];
  u[0] = u0;
  int k = 1;
  while (k < _acc - 1){
    Vector4d v = AM_one_step(u,dt,mu,k+1,table);
    u[k] = v;
    k++;
  }
  std::ofstream os;
  std::string ss = "AM_" + std::to_string(_acc) + "_" + std::to_string(N) + ".m";
  const char *s = ss.c_str();
  os.open(s);
  os << "a=[\n";
  for (int i = 0; i < _acc - 1 ; i++)
    os << u[i](0) << "," << u[i](1) << ";\n";
  const Info& AM_info = table.get_info(_acc);
  for (int i = 0; i < N - _acc + 2; i++){
    Vector4d v = AM_one_step(u,dt,mu,_acc,AM_info);
    for (int j = 0; j < _acc - 2 ; j++)
      u[j] = u[j+1];
    u[_acc-2] = v;
    os << v(0) << "," << v(1) << ";\n";
  }
  os << "];\n";
  os << "x = a(:,1);";
  os << "y = a(:,2);";
  os << "plot(x,y)";
  os.close();
  return u[_acc-2];
}

Vector4d AM_method(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table){
  clock_t t1 = clock();
  Vector4d v = AM_method(u0,dt,mu,_acc,N,table);
  clock_t t2 = clock();
  time = (double)(t2 - t1) / CLOCKS_PER_SEC * 1000;
  return v;
}


Vector4d BDF_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info& info){
  Vector4d b;
  b << 0,0,0,0;
  for ( int i = 1; i <= _acc; i++)
    b += u[_acc-i]*info.get_coe(i-1);
  Vector4d v = Newton(u[_acc-1],mu,-dt*info.get_coe(_acc),-b);
  return v;
}

Vector4d BDF_one_step(Vector4d* u, const double dt, const double mu, int _acc, const Info_Table& table){
  const Info& BDF_info=table.get_info(_acc);
  return BDF_one_step(u,dt,mu,_acc,BDF_info);
}

Vector4d BDF_method(Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table){
  Vector4d* u;
  u = new Vector4d [_acc];
  u[0] = u0;
  int k = 1;
  while (k < _acc){
    Vector4d v = BDF_one_step(u,dt,mu,k,table);
    u[k] = v;
    k++;
  }
  std::ofstream os;
  std::string ss = "BDF_" + std::to_string(_acc) + "_" + std::to_string(N) + ".m";
  const char *s = ss.c_str();
  os.open(s);
  os << "a=[\n";
  for (int i = 0; i < _acc ; i++)
    os << u[i](0) << "," << u[i](1) << ";\n";
  const Info& BDF_info = table.get_info(_acc);
  for (int i = 0; i < N - _acc + 1; i++){
    Vector4d v = BDF_one_step(u,dt,mu,_acc,BDF_info);
    for (int j = 0; j < _acc - 1 ; j++)
      u[j] = u[j+1];
    u[_acc-1] = v;
    os << v(0) << "," << v(1) << ";\n";
  }
  os << "];\n";
  os << "x = a(:,1);";
  os << "y = a(:,2);";
  os << "plot(x,y)";
  os.close();
  return u[_acc-1];
}

Vector4d BDF_method(double& time, Vector4d u0, const double dt, const double mu, int _acc, int N, const Info_Table& table){
  clock_t t1 = clock();
  Vector4d v = BDF_method(u0,dt,mu,_acc,N,table);
  clock_t t2 = clock();
  time = (double)(t2 - t1) / CLOCKS_PER_SEC * 1000;
  return v;
}


Vector4d RK_one_step(Vector4d u, const double dt, const double mu){
  Vector4d k1 = dt * f(u,mu);
  Vector4d k2 = dt * f(u + 0.5 * k1,mu);
  Vector4d k3 = dt * f(u + 0.5 * k2,mu);
  Vector4d k4 = dt * f(u + k3,mu);
  return u + (k1 + 2 * k2 + 2 * k3 + k4)/6;
}


Vector4d RK_method(Vector4d u0, const double dt, const double mu, int N){
  Vector4d u = u0;
  std::ofstream os;
  std::string ss = "RK_" + std::to_string(N) + ".m";
  const char *s = ss.c_str();
  os.open(s);
  os << "a=[\n";
  for (int i = 0; i < N ; i++){
    Vector4d v = RK_one_step(u,dt,mu);
    u = v;
    os << u(0) << "," << u(1) << ";\n";
  }
  os << "];\n";
  os << "x = a(:,1);";
  os << "y = a(:,2);";
  os << "plot(x,y)";
  os.close();
  return u;
}

Vector4d RK_method(double& time, Vector4d u0, const double dt, const double mu, int N){
  clock_t t1 = clock();
  Vector4d v = RK_method(u0,dt,mu,N);
  clock_t t2 = clock();
  time = (double)(t2 - t1) / CLOCKS_PER_SEC * 1000;
  return v;
}

double err_initial(double& time, Vector4d u0, const double dt, const double mu, const int _acc, double T, const Info_Table& table, int type){
  int N = round(T/dt);
  Vector4d v;
  double err = -1;
  switch(type){
  case 1:
    v = AB_method(time,u0,dt,mu,_acc,N,table);
    err = (u0.head(2)-v.head(2)).norm();
    return err;
  case 2:
    v = AM_method(time,u0,dt,mu,_acc,N,table);
    err = (u0.head(2)-v.head(2)).norm();
    return err;
  case 3:
    v = BDF_method(time,u0,dt,mu,_acc,N,table);
    err = (u0.head(2)-v.head(2)).norm();
    return err;
  default:
    std::cerr << "No matching type!" << std::endl;
    exit(-1);    
  }
}

double err_initial(double& time, Vector4d u0, const double dt, const double mu, double T, int type){
  if ( type == 4){
    int N = round(T/dt);
    Vector4d v = RK_method(time,u0,dt,mu,N);
    double err = (u0.head(2)-v.head(2)).norm();
    return err;
  }
  else{
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
}


#else
//do nothing
#endif
