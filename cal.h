#ifndef __MGGS_CAL__
#define __MGGS_CAL__
#include "Info.h"
#include <eigen/Eigen/Dense>
#include <cmath>
#include <ctime>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace Eigen;

Vector4d f(Vector4d u, const double mu);
Vector4d F(Vector4d u, const double mu, const double k, Vector4d b);
Matrix4d DF(Vector4d u, const double mu, const double k);
Vector4d Newton(Vector4d u0, const double mu, const double k, Vector4d b);

Vector4d initial_load(const std::string &_file, double& _mu);
Vector4d initial_load(const std::string &_file, double& _mu, double& _T);

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

double max_norm(const Vector4d u1, const Vector4d u2);
double max_norm(const Vector2d u1, const Vector2d u2);

double err_initial(double& time, Vector4d u0, const double dt, const double mu, const int _acc, double T, const Info_Table& table, int type);
double err_initial(double& time, Vector4d u0, const double dt, const double mu, double T, int type);

Vector2d extrapolate(const Vector2d u1, const Vector2d u2, int j);

double err_richardson(double tol, double& time, Vector4d u0, double dt, const double mu, const int _acc, int N, const Info_Table& table, int type);
double err_richardson(double tol, double& time, Vector4d u0, double dt, const double mu, int N, int type);


std::pair<double,double> grid_refine_err1(Vector4d u0, double dt, const double mu, const int _acc, const double T, const Info_Table& table, int type);
std::pair<double,double> grid_refine_err1(Vector4d u0, double dt, const double mu, const double T, int type);
double grid_refine_err2(double tol, Vector4d u0, double dt, const double mu, const int _acc, int N,const Info_Table& table, int type);
double grid_refine_err2(double tol, Vector4d u0, double dt, const double mu, int N, int type);


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

double max_norm(const Vector4d u1, const Vector4d u2){
  double y1 = fabs(u1(0)-u2(0));
  double y2 = fabs(u1(1)-u2(1));
  if (y1 >= y2)
    return y1;
  else
    return y2;
}

double max_norm(const Vector2d u1, const Vector2d u2){
  double y1 = fabs(u1(0)-u2(0));
  double y2 = fabs(u1(1)-u2(1));
  if (y1 >= y2)
    return y1;
  else
    return y2;
}

Vector4d initial_load(const std::string &_file, double& _mu){
  std::fstream data(_file);
  data >>_mu;
  Vector4d u;
  for (int i = 0 ; i < 4 ; i++)
    data >> u(i);
  data.close();
  return u;
}

Vector4d initial_load(const std::string &_file, double& _mu, double& _T){
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
    Vector4d v = RK_one_step(u[k-1],dt,mu);
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
  os << "x = a(:,1);\n";
  os << "y = a(:,2);\n";
  os << "plot(x,y);";
  os << "title('Adams Bashforth,p = " << _acc << ",dt = " << dt << ",N = " << N << "');\n";
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
    Vector4d v = RK_one_step(u[k-1],dt,mu);
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
  os << "x = a(:,1);\n";
  os << "y = a(:,2);\n";
  os << "plot(x,y);";
  os << "title('Adams Moulton,p = " << _acc << ",dt = " << dt << ",N = " << N << "');\n";
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
    Vector4d v = RK_one_step(u[k-1],dt,mu);
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
  os << "x = a(:,1);\n";
  os << "y = a(:,2);\n";
  os << "plot(x,y);";
  os << "title('BDFs,p = " << _acc << ",dt = " << dt << ",N = " << N << "');\n";
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
  os << "x = a(:,1);\n";
  os << "y = a(:,2);\n";
  os << "plot(x,y);";
  os << "title('Runge Kutta,dt = " << dt << ",N = " << N << "');\n";
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
  int N = floor(T/dt);
  double rt = T - N*dt;
  Vector4d v,v1;
  double err = -1;
  switch(type){
  case 1:
    v = AB_method(time,u0,dt,mu,_acc,N,table);
    v1 = RK_one_step(v,rt,mu);
    err = max_norm(u0,v1);
    return err;
  case 2:
    v = AM_method(time,u0,dt,mu,_acc,N,table);
    v1 = RK_one_step(v,rt,mu);
    err = max_norm(u0,v1);
    return err;
  case 3:
    v = BDF_method(time,u0,dt,mu,_acc,N,table);
    v1 = RK_one_step(v,rt,mu);
    err = max_norm(u0,v1);
    return err;
  default:
    std::cerr << "No matching type!" << std::endl;
    exit(-1);    
  }
}

double err_initial(double& time, Vector4d u0, const double dt, const double mu, double T, int type){
  if ( type == 4){
    int N = floor(T/dt);
    double rt = T - N*dt;
    Vector4d v = RK_method(time,u0,dt,mu,N);
    Vector4d v1;
    v1 = RK_one_step(v,rt,mu);
    double err = max_norm(u0,v1);
    return err;
  }
  else{
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
}

Vector2d extrapolate(const Vector2d u1, const Vector2d u2, int j){
  double tmp = pow(4,j);
  Vector2d v;
  v(0) = (tmp*u1(0) - u2(0))/(tmp - 1);
  v(1) = (tmp*u1(1) - u2(1))/(tmp - 1);
  return v;
}

double err_richardson(double tol, double& time, Vector4d u0,double dt, const double mu, const int _acc, int N, const Info_Table& table, int type){
  Vector4d (*pf1)(Vector4d,const double,const double,int,int,const Info_Table&);
  Vector4d (*pf2)(double&,Vector4d,const double,const double,int,int,const Info_Table&);
  switch(type){
  case 1:
    pf1 = AB_method;
    pf2 = AB_method;
    break;
  case 2:
    pf1 = AM_method;
    pf2 = AM_method;
    break;
  case 3:
    pf1 = BDF_method;
    pf2 = BDF_method;
    break;
  default:
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
  double err;
  int i;
  Vector2d uu;
  typedef Matrix<Vector2d,Dynamic,1> VectorXv;
  VectorXv v1(1),v2(1);
  Vector4d u = pf2(time,u0,dt,mu,_acc,N,table);
  v1(0) << u(0),u(1);
  uu = v1(0);
  dt = 0.5*dt;
  N = 2*N;
  for ( i = 2; i < 7 ; i++){
    v2.resize(i,1);
    u = pf1(u0,dt,mu,_acc,N,table);
    v2(0) << u(0),u(1);
    for (int j = 1 ; j < i ; j++)
      v2(j) = extrapolate(v2(j-1),v1(j-1),j);
    err = max_norm(v2(i-1),v1(i-2));
    if (err < tol)
      break;
    v1.resize(i,1);
    v1 = v2;
    dt = 0.5*dt;
    N = 2*N;
  }
  return max_norm(uu,v2(i-1));
}

double err_richardson(double tol, double& time, Vector4d u0, double dt, const double mu, int N, int type){
  if (type == 4){
    double err;
    typedef Matrix<Vector2d,Dynamic,1> VectorXv;
    int i;
    Vector2d uu;
    VectorXv v1(1),v2(1);
    Vector4d u = RK_method(time,u0,dt,mu,N);
    v1(0) << u(0),u(1);
    uu = v1(0);
    dt = 0.5*dt;
    N = 2*N;
    for (i = 2; i < 7 ; i++){
      v2.resize(i,1);
      u = RK_method(u0,dt,mu,N);
      v2(0) << u(0),u(1);
      for (int j = 1 ; j < i ; j++)
	v2(j) = extrapolate(v2(j-1),v1(j-1),j);
      err = max_norm(v2(i-1),v1(i-2));
      if (err < tol)
	break;
      v1.resize(i,1);
      v1 = v2;
      dt = 0.5*dt;
      N = 2*N;
    }
    return max_norm(uu,v2(i-1));
  }
  else{
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
}

std::pair<double,double> grid_refine_err1(Vector4d u0, double dt, const double mu, const int _acc, const double T, const Info_Table& table, int type){
  int N = 4;
  std::string ss,Method;
  switch (type){
  case 1:
    ss = "AB_" + std::to_string(_acc) + "_Init1_analysis.m";
    Method = "Adams Bashforth";
    break;
  case 2:
    ss = "AM_" + std::to_string(_acc) + "_Init1_analysis.m";
    Method = "Adams Moulton";
    N--;
    break;
  case 3:
    ss = "BDF_" + std::to_string(_acc) + "_Init1_analysis.m";
    Method = "BDFs";
    N--;
    break;
  default:
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
  std::ofstream os;
  const char *s = ss.c_str();
  os.open(s);
  os << "x=[\n";
  double time;
  double* err = new double[N];
  for (int i = 0; i < N ; i++){
    err[i] = err_initial(time,u0,dt,mu,_acc,T,table,type);
    int step = round(T/dt);
    os << dt << "," << step << "," << time << "," << std::setprecision(16) << err[i] << ";\n";
    dt = 0.5 * dt;
  }
  os << "];\n";
  os << "dt = x(:,1);\n";
  os << "steps = x(:,2);\n";
  os << "CPU_time = x(:,3);\n";
  os << "err = x(:,4);\n";
  double cr = 1,constant = 1;
  const double l2 = log(2);
  for (int i = 0; i < N - 1 ; i++){
    double mult = std::max(1.0,log(err[i]/err[i+1])/l2);
    cr *= mult;
  }
  cr = pow(cr,1.0/(N - 1));
  for (int i = N - 1; i >= 0 ; i--){
    dt = 2 * dt;
    double mult = err[i]/pow(dt,cr);
    constant *= mult;
  }
  constant = pow(constant,1.0/N);
  os << "cr=" << cr << ";\n";
  os << "subplot(1,3,1);\n";
  os << "plot(dt,CPU_time,'r-*');xlabel('dt');ylabel('CPU time(ms)');"
     << "title('" << Method << ",p = " << _acc << "');\n";
  os << "subplot(1,3,2);\n";
  os << "plot(dt,err,'m-*');xlabel('dt');ylabel('err');"
     << "title('" << Method << ",p = " << _acc << "');\n";
  os << "subplot(1,3,3);\n";
  os << "plot(dt.^cr,err,'b-*');xlabel(strcat('dt^p,p=',num2str(cr)));ylabel('err');"
     << "title('" << Method << ",p = " << _acc << "');\n";
  os.close();
  std::pair<double,double> res(cr,constant);
  return res;
}

std::pair<double,double> grid_refine_err1(Vector4d u0, double dt, const double mu, const double T, int type){
  const int N = 4;
  std::string ss = "RK_Init1_analysis.m";
  std::ofstream os;
  const char *s = ss.c_str();
  os.open(s);
  os << "x=[\n";
  double time;
  double* err = new double[N];
  for (int i = 0; i < N ; i++){
    err[i] = err_initial(time,u0,dt,mu,T,4);
    int step = round(T/dt);
    os << dt << "," << step << "," << time << "," << std::setprecision(16) << err[i] << ";\n";
    dt = 0.5 * dt;
  }
  os << "];\n";
  os << "dt = x(:,1);\n";
  os << "steps = x(:,2);\n";
  os << "CPU_time = x(:,3);\n";
  os << "err = x(:,4);\n";
  double cr = 1,constant = 1;
  const double l2 = log(2);
  for (int i = 0; i < N - 1 ; i++){
    double mult = std::max(1.0,log(err[i]/err[i+1])/l2);
    cr *= mult;
  }
  cr = pow(cr,1.0/(N - 1));
  for (int i = N - 1; i >= 0 ; i--){
    dt = 2 * dt;
    double mult = err[i]/pow(dt,cr);
    constant *= mult;
  }
  constant = pow(constant,1.0/N);
  os << "cr=" << cr << ";\n";
  os << "subplot(1,3,1);\n";
  os << "plot(dt,CPU_time,'r-*');xlabel('dt');ylabel('CPU time(ms)');"
     << "title('Runge Kutta');\n";
  os << "subplot(1,3,2);\n";
  os << "plot(dt,err,'m-*');xlabel('dt');ylabel('err');"
     << "title('Runge Kutta');\n";
  os << "subplot(1,3,3);\n";
  os << "plot(dt.^cr,err,'b-*');xlabel(strcat('dt^p,p=',num2str(cr)));ylabel('err');"
     << "title('Runge Kutta');\n";
  os.close();
  std::pair<double,double> res(cr,constant);
  return res;
}

double grid_refine_err2(double tol, Vector4d u0, double dt, const double mu, const int _acc, int N,const Info_Table& table, int type){
  Vector4d (*pf)(double&,Vector4d,const double,const double,int,int,const Info_Table&);
  std::string ss,Method;
  switch(type){
  case 1:
    pf = AB_method;
    Method = "Adams Bashforth";
    ss = "AB_" + std::to_string(_acc) + "_Init2_analysis.m";
    break;
  case 2:
    pf = AM_method;
    Method = "Adams Moulton";
    ss = "AM_" + std::to_string(_acc) + "_Init2_analysis.m";
    break;
  case 3:
    pf = BDF_method;
    Method = "BDFs";
    ss = "BDF_" + std::to_string(_acc) + "_Init2_analysis.m";
    break;
  default:
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
  std::ofstream os;
  const char *s = ss.c_str();
  os.open(s);
  os << "x=[\n";
  double err,time;
  int i;
  Vector2d uu;
  typedef Matrix<Vector2d,Dynamic,1> VectorXv;
  VectorXv v1(1),v2(1);
  Vector4d u = pf(time,u0,dt,mu,_acc,N,table);
  os << dt << "," << N << "," << time << "," << std::setprecision(16) << u[0] << "," << std::setprecision(16) << u[1] << ";\n";
  v1(0) << u(0),u(1);
  uu = v1(0);
  dt = 0.5*dt;
  N = 2*N;
  for ( i = 2; i < 7 ; i++){
    v2.resize(i,1);
    u = pf(time,u0,dt,mu,_acc,N,table);
    if ( i < 5 ){
    os << dt << "," << N << "," << time << "," << std::setprecision(16) << u[0] << "," << std::setprecision(16) << u[1] << ";\n";
    }
    v2(0) << u(0),u(1);
    for (int j = 1 ; j < i ; j++)
      v2(j) = extrapolate(v2(j-1),v1(j-1),j);
    err = max_norm(v2(i-1),v1(i-2));
    v1.resize(i,1);
    v1 = v2;
    dt = 0.5*dt;
    N = 2*N;
    if (err < tol)
      break;
  }
  os << "];\n";
  os << "dt = x(:,1);\n";
  os << "steps = x(:,2);\n";
  os << "CPU_time = x(:,3);\n";
  os << "u = [" << v2(i-2)(0) << "," << v2(i-2)(1) << "];\n";
  os << "[err,index]=max(abs(x(:,4:5)-u),[],2);\n";
  os << "len=length(err);\n";
  os << "temp=err(2:end);temp(len)=0;\n";
  os << "result=max(1,log2(err(1:len-1)./temp(1:len-1)));\n";
  os << "cr=geomean(result);\n";
  os << "result2=err./(dt.^cr);\n";
  os << "constant=geomean(result2);\n";
  os << "subplot(1,3,1);\n";
  os << "plot(dt,CPU_time,'r-*');xlabel('dt');ylabel('CPU time(ms)');"
     << "title('" << Method << ",p = " << _acc << "');\n";
  os << "subplot(1,3,2);\n";
  os << "plot(dt,err,'m-*');xlabel('dt');ylabel('err');"
     << "title('" << Method << ",p = " << _acc << "');\n";
  os << "subplot(1,3,3);\n";
  os << "plot(dt.^cr,err,'b-*');xlabel(strcat('dt^p,p=',num2str(cr)));ylabel('err');"
     << "title('" << Method << ",p = " << _acc << "');\n";
  os.close();
  return max_norm(uu,v2(i-2));
}

double grid_refine_err2(double tol, Vector4d u0, double dt, const double mu, int N, int type){
   if (type == 4){
    double err,time;
    std::string ss = "RK_Init2_analysis.m";
    std::ofstream os;
    const char *s = ss.c_str();
    os.open(s);
    os << "x=[\n";
    typedef Matrix<Vector2d,Dynamic,1> VectorXv;
    int i;
    Vector2d uu;
    VectorXv v1(1),v2(1);
    Vector4d u = RK_method(time,u0,dt,mu,N);
    os << dt << "," << N << "," << time << "," << std::setprecision(16) << u[0] << "," << std::setprecision(16) << u[1] << ";\n";
    v1(0) << u(0),u(1);
    uu = v1(0);
    dt = 0.5*dt;
    N = 2*N;
    for (i = 2; i < 7 ; i++){
      v2.resize(i,1);
      u = RK_method(time,u0,dt,mu,N);
      if (i < 5){
      os << dt << "," << N << "," << time << "," << std::setprecision(16) << u[0] << "," << std::setprecision(16) << u[1] << ";\n";
      }
      v2(0) << u(0),u(1);
      for (int j = 1 ; j < i ; j++)
	v2(j) = extrapolate(v2(j-1),v1(j-1),j);
      err = max_norm(v2(i-1),v1(i-2));
      if (err < tol)
	break;
      v1.resize(i,1);
      v1 = v2;
      dt = 0.5*dt;
      N = 2*N;
    }
    os << "];\n";
    os << "dt = x(:,1);\n";
    os << "steps = x(:,2);\n";
    os << "CPU_time = x(:,3);\n";
    os << "u = [" << v2(i-1)(0) << "," << v2(i-1)(1) << "];\n";
    os << "[err,index]=max(abs(x(:,4:5)-u),[],2);\n";
    os << "len=length(err);\n";
    os << "temp=err(2:end);temp(len)=0;\n";
    os << "result=max(1,log2(err(1:len-1)./temp(1:len-1)));\n";
    os << "cr=geomean(result);\n";
    os << "result2=err./(dt.^cr);\n";
    os << "constant=geomean(result2);\n";
    os << "subplot(1,3,1);\n";
    os << "plot(dt,CPU_time,'r-*');xlabel('dt');ylabel('CPU time(ms)');"
       << "title('Runge Kutta');\n";
    os << "subplot(1,3,2);\n";
    os << "plot(dt,err,'m-*');xlabel('dt');ylabel('err');"
       << "title('Runge Kutta');\n";
    os << "subplot(1,3,3);\n";
    os << "plot(dt.^cr,err,'b-*');xlabel(strcat('dt^p,p=',num2str(cr)));ylabel('err');"
       << "title('Runge Kutta');\n";
    os.close();
    return max_norm(uu,v2(i-1));
  }
  else{
    std::cerr << "No matching type!" << std::endl;
    exit(-1);
  }
}

#else
//do nothing
#endif
