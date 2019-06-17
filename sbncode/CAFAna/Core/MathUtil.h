#ifndef MATHUTIL_H
#define MATHUTIL_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file MathUtil.h                                                     //
//                                                                      //
// Simple mathematical functions, initially to prevent the abuse of     //
// pow()                                                                //
// <bckhouse@caltech.edu>						//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>

namespace util{

  /// \name Simple mathematical functions
  //@{

  /// More efficient square function than pow(x,2)
  template<class T> inline T sqr(T x){return x*x;}

  /// More efficient cube function than pow(x,3)
  template<class T> inline T cube(T x){return x*x*x;}

  /// More efficient exponentiation function than pow(x,n) for small n
  template<class T> inline T ipow(T x, unsigned int n)
  {
    T ret = 1; 
    if (n == 0) return ret;
    for(unsigned int i = 1; i <= n; ++i) ret *= x; 
    return ret;
  }


  /// 2D Euclidean distance
  inline double pythag(double x, double y)
  {
    return sqrt(sqr(x)+sqr(y));
  }

  /// 3D Euclidean distance
  inline double pythag(double x, double y, double z)
  {
    return sqrt(sqr(x)+sqr(y)+sqr(z));
  }

} // end namespace

#endif // MATHUTIL_H
