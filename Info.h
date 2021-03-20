#ifndef __MGGS_INFO__
#define __MGGS_INFO__
#include <iostream>
#include <fstream>

class Info{
 private:
  int accuracy;
  int type;
  int n_coe;
  double* coe_info;
 public:
  Info();
  ~Info();
  void set(int _acc, int _type, double* _coe_info);
  int get_acc() const;
  int get_n_coe() const;
  double get_coe(int _idx) const;
};









#else
//do nothing
#endif

