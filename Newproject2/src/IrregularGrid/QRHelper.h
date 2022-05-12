#ifndef QRHELPER_H_
#define QRHELPER_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <Core/Vec.h>

using Point = Vec<double,2>;

std::vector<double> QRSolve(double* A, const double* b, int m, int n);

std::vector<Point> QRSolve(double* A, double* b1, double* b2, int m,
                           int n);

void QR(double* A, double* Q, int m, int n);

std::vector<double> SolveQRUpper(double* R, double* c, int m, int n);

std::vector<double> QRSolve(double* A, const double* b, int m, int n){
  double* Q = new double[m*m];
  QR(A,Q,m,n);
  double * c = new double[n];
  for (int k = 0 ; k < n ; k++){
    c[k] = 0.0;
    for (int i = 0 ; i < m ; i++)
      c[k] += Q[k+i*m]*b[i];
  }
  std::vector<double> res = SolveQRUpper(A,c,m,n);
  delete [] c;
  delete [] Q;
  return res;
}

std::vector<Point> QRSolve(double* A, double* b1, double* b2, int m,
                           int n){
  double* Q = new double[m*m];
  QR(A,Q,m,n);
  double * c1 = new double[n];
  double * c2 = new double[n];
  for (int k = 0 ; k < n ; k++){
    c1[k] = c2[k] = 0.0;
    for (int i = 0 ; i < m ; i++){
      c1[k] += Q[k+i*m]*b1[i];
      c2[k] += Q[k+i*m]*b2[i];
    }
  }
  std::vector<double> res1 = SolveQRUpper(A,c1,m,n);
  std::vector<double> res2 = SolveQRUpper(A,c2,m,n);
  std::vector<Point> res;
  for (int i = 0 ; i < n ; i++)
    res.push_back(Point{res1[i],res2[i]});
  delete [] c1;
  delete [] c2;
  delete [] Q;
  return res;
}

// Get Q^T by Q and R by A
void QR(double* A, double* Q, int m, int n){
  int i, j, k, nn, jj;
  double u, alpha, w, t;
  
  if (m < n){
    std::cout << "Wrong! m < n." << std::endl;
    exit(1);
  }
  for (i = 0; i <= m - 1; i++)
    for (j = 0; j <= m - 1; j++){
      Q[i+j*m] = 0.0;
      if (i == j) Q[i+j*m] = 1.0;
    }
  nn = n;
  if (m == n) nn = m - 1;
  for (k = 0; k <= nn - 1; k++){
    u = 0.0;
    for (i = k; i <= m - 1; i++){
      w = fabs(A[i+k*m]);
      if (w > u) u = w;
    }
    alpha = 0.0;
    for (i = k; i <= m - 1; i++){
      t = A[i+k*m] / u; alpha = alpha + t * t;
    }
    if (A[k+k*m] > 0.0) u = -u;
    alpha = u * sqrt(alpha);
    if (fabs(alpha) + 1.0 == 1.0){
      std::cout << "Wrong!" << std::endl;
      exit(1);
    }
    u = sqrt(2.0*alpha*(alpha - A[k+k*m]));
    if ((u + 1.0) != 1.0){
      A[k+k*m] = (A[k+k*m] - alpha) / u;
      for (i = k + 1; i <= m - 1; i++) A[i+k*m] = A[i+k*m] / u;
      for (j = 0; j <= m - 1; j++){
        t = 0.0;
        for (jj = k; jj <= m - 1; jj++)
          t = t + A[jj+k*m] * Q[jj+j*m];
        for (i = k; i <= m - 1; i++)
          Q[i+j*m] = Q[i+j*m] - 2.0*t*A[i+k*m];
      }
      for (j = k + 1; j <= n - 1; j++){
        t = 0.0;
        for (jj = k; jj <= m - 1; jj++)
          t = t + A[jj+k*m] * A[jj+j*m];
        for (i = k; i <= m - 1; i++)
          A[i+j*m] = A[i+j*m] - 2.0*t*A[i+k*m];
      }
      A[k+k*m] = alpha;
      for (i = k + 1; i <= m - 1; i++)  A[i+k*m] = 0.0;
    }
  }
  // for (i = 0; i <= m - 2; i++)
  //   for (j = i + 1; j <= m - 1; j++){
  //     t = Q[i+j*m]; Q[i+j*m] = Q[j+i*m]; Q[j+i*m] = t;
  //   }
}

std::vector<double> SolveQRUpper(double* R, double* c, int m, int n){
  std::vector<double> res(n);
  res[n-1] = c[n-1] / R[n-1+(n-1)*m];
  for (int i = n - 2 ; i >= 0 ; i--){
    double tmp = c[i];
    for (int j = i + 1 ; j < n ; j++)
      tmp -= R[i+j*m]*res[j];
    res[i] = tmp / R[i+i*m];
  }
  return res;
}


#endif
