#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "RegularGrid/RectDomain.h"
#include "Core/Tensor.h"

template <int Dim> class Interpolator;

template <int Dim>
class Interpolator{
public:
  Interpolator(const RectDomain<Dim>& afineDomain, const RectDomain<Dim>&
                acoarseDomain);
  virtual ~Interpolator() = default;
  virtual void apply(const Tensor<Real,Dim>& data,
                     Tensor<Real,Dim>& res) const = 0;
protected:
  RectDomain<Dim> fineDomain;
  RectDomain<Dim> coarseDomain;
};

template <int Dim>
Interpolator<Dim>::Interpolator(const RectDomain<Dim>& afineDomain, const RectDomain<Dim>&
                                acoarseDomain){
  assert(afineDomain.getCentering() == NodeCentered &&
  acoarseDomain.getCentering() == NodeCentered);
  fineDomain = afineDomain;
  coarseDomain = acoarseDomain;
}

template <int Dim> class LinearInterpolator;

template <>
class LinearInterpolator<2>: public Interpolator<2>{
public:
  using BaseClass = Interpolator<2>;

  LinearInterpolator(const RectDomain<2>& afineDomain, const RectDomain<2>&
                     acoarseDomain):Interpolator(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,2>& data,
             Tensor<Real,2>& res) const{
    const auto& bx = fineDomain;
    const Vec<int, 2>& lo = bx.lo();
    const Vec<int, 2>& hi = bx.hi();
    loop_box_2(bx, i, j) {
      if (i % 2 == 0 && j % 2 == 0)
        res(i, j) = data(i / 2, j / 2);
      else if (i % 2 == 0 && j % 2 != 0)
        res(i, j) = (data(i / 2, j / 2) + data(i / 2, j / 2 + 1)) / 2.0;
      else if (i % 2 != 0 && j % 2 == 0)
        res(i, j) = (data(i / 2, j / 2) + data(i / 2 + 1, j / 2)) / 2.0;
      else
        res(i, j) = (data(i / 2, j / 2) + data(i / 2 + 1, j / 2) +
                     data(i / 2, j / 2 + 1) + data(i / 2 + 1, j / 2 + 1)) /
                    4.0;
    }
  }
};

template <>
class LinearInterpolator<1>: public Interpolator<1>{
public:
  using BaseClass = Interpolator<1>;

  LinearInterpolator(const RectDomain<1>& afineDomain, const RectDomain<1>&
                     acoarseDomain):Interpolator(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,1>& data,
             Tensor<Real,1>& res) const{
    const auto& bx = fineDomain;
    const Vec<int,1>& lo = bx.lo();
    const Vec<int, 1>& hi = bx.hi();
    loop_box_1(bx, i) {
      if (i % 2 == 0)
        res(i) = data(i / 2);
      else
        res(i) = (data(i / 2) + data(i / 2 + 1)) / 2.0;
    }
  }
};

template <int Dim> class QuadraticInterpolator;

template <>
class QuadraticInterpolator<2>: public Interpolator<2>{
public:
  using BaseClass = Interpolator<2>;

  QuadraticInterpolator(const RectDomain<2>& afineDomain, const RectDomain<2>&
                        acoarseDomain):Interpolator(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,2>& data,
             Tensor<Real,2>& res) const{
    const auto& bx = fineDomain;
    const Vec<int,2>& sz = bx.size();
    const Vec<int,2>& lo = bx.lo();
    const Vec<int, 2>& hi = bx.hi();
    loop_box_2(bx, i, j) {
      if (i % 2 == 0 && j % 2 == 0)
        res(i, j) = data(i / 2, j / 2);
      else if (i % 2 == 0 && j % 2 != 0) {
        if (j / 2 + 2 < sz[1] - 4)
          res(i, j) = (3 * data(i / 2, j / 2) + 6 * data(i / 2, j / 2 + 1) -
                       data(i / 2, j / 2 + 2)) /
                      8.0;
        else
          res(i, j) = (3 * data(i / 2, j / 2 + 1) + 6 * data(i / 2, j / 2) -
                       data(i / 2, j / 2 - 1)) /
                      8.0;
      } else if (i % 2 != 0 && j % 2 == 0) {
        if (i / 2 + 2 < sz[0] - 4)
          res(i, j) = (3 * data(i / 2, j / 2) + 6 * data(i / 2 + 1, j / 2) -
                       data(i / 2 + 2, j / 2)) /
                      8.0;
        else
          res(i, j) = (3 * data(i / 2 + 1, j / 2) + 6 * data(i / 2, j / 2) -
                       data(i / 2 - 1, j / 2)) /
                      8.0;
      } else {
        Real tmp1, tmp2;
        if (j / 2 + 2 < sz[1] - 4)
          tmp1 = (3 * data(i / 2, j / 2) + 6 * data(i / 2, j / 2 + 1) -
                  data(i / 2, j / 2 + 2)) /
                 8.0;
        else
          tmp1 = (3 * data(i / 2, j / 2 + 1) + 6 * data(i / 2, j / 2) -
                  data(i / 2, j / 2 - 1)) /
                 8.0;
        if (i / 2 + 2 < sz[0] - 4)
          tmp2 = (3 * data(i / 2, j / 2) + 6 * data(i / 2 + 1, j / 2) -
                  data(i / 2 + 2, j / 2)) /
                 8.0;
        else
          tmp2 = (3 * data(i / 2 + 1, j / 2) + 6 * data(i / 2, j / 2) -
                  data(i / 2 - 1, j / 2)) /
                 8.0;
        res(i, j) = (tmp1 + tmp2) / 2.0;
      }
    }
  }
};

template <>
class QuadraticInterpolator<1>: public Interpolator<1>{
public:
  using BaseClass = Interpolator<1>;

  QuadraticInterpolator(const RectDomain<1>& afineDomain, const RectDomain<1>&
                        acoarseDomain):Interpolator(afineDomain,acoarseDomain){}
  
  void apply(const Tensor<Real,1>& data,
             Tensor<Real,1>& res) const{
    const auto& bx = fineDomain;
    const Vec<int,1>& sz = bx.size();
    const Vec<int,1>& lo = bx.lo();
    const Vec<int, 1>& hi = bx.hi();
    loop_box_1(bx, i) {
      if (i % 2 == 0)
        res(i) = data(i / 2);
      else {
        if (i / 2 + 2 < sz[0] - 4)
          res(i) =
              (3 * data(i / 2) + 6 * data(i / 2 + 1) - data(i / 2 + 2)) / 8.0;
        else
          res(i) =
              (3 * data(i / 2 + 1) + 6 * data(i / 2) - data(i / 2 - 1)) / 8.0;
      }
    }
  }
};


  


#endif // _INTERPOLATOR_H_
