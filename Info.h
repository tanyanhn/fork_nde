/**
 * @file   Info.h
 * @author Lizhiqi <3180103041@zju.edu.cn>
 * @date   Tue Apr 27 13:53:48 2021
 * 
 * @brief  Data structure for formula of
 * linear multistep methods
 * 
 */

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
  int get_type() const;
  int get_n_coe() const;
  double get_coe(int _idx) const;
};

class Info_Table{
 private:
  int n_table;
  int type;
  Info* table;
 public:
  Info_Table();
  ~Info_Table();
  void load_data(const char _file[],const int _type);
  void set_type(const int _type);
  int get_n_table() const;
  int get_type() const;
  const Info& get_info(int _acc) const;
};







#else
//do nothing
#endif

