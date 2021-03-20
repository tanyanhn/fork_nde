#include "Info.h"
#include <eigen/Eigen/Dense>

using namespace Eigen;

Vector4d f(Vector4d u);

Vector4d f(Vector4d u){
  u(0) = u(0) * 2;
  u(2) += u(1);
  return u;
}
