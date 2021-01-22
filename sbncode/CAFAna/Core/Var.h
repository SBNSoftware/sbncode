#pragma once

// This file defines the basic Var object. For specific variables, and examples
// of how to implement your own, see Vars.h

#include <functional>
#include <set>
#include <string>

#include "CAFAna/Core/Binning.h"

#include "StandardRecord/Proxy/FwdDeclare.h"

namespace caf{class SRSpill; class SRSpillTruthBranch; class SRSlice;}

namespace ana
{

  /// Most useful for combining weights.
  // need to declare these functions first
  // (this needed because of reasons detailed in https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Making_New_Friends;
  //  we want a 'one-to-one' relationship.)
  template<class T> class _Var;
  template<class T> _Var<T> operator*(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Var<T> operator/(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Var<T> operator+(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Var<T> operator-(const _Var<T>& a, const _Var<T>& b);

  /// Template for Var and SpillVar
  template<class T> class _Var
  {
  public:
    /// The type of the function part of a var
    typedef double (VarFunc_t)(const T* sr);

    /// std::function can wrap a real function, function object, or lambda
    _Var(const std::function<VarFunc_t>& fun);

    /// Allows a variable to be called with double value = myVar(sr) syntax
    double operator()(const T* sr) const
    {
      return fFunc(sr);
    }

    /// Vars with the same definition will have the same ID
    int ID() const {return fID;}

    static int MaxID();
  protected:
    // Give this guy access to the constructor that sets ID
    friend _Var<T> operator*<>(const _Var<T>& a, const _Var<T>& b);
    friend _Var<T> operator/<>(const _Var<T>& a, const _Var<T>& b);
    friend _Var<T> operator+<>(const _Var<T>& a, const _Var<T>& b);
    friend _Var<T> operator-<>(const _Var<T>& a, const _Var<T>& b);

    _Var(const std::function<VarFunc_t>& fun, int id)
      : fFunc(fun), fID(id)
    {
    }

    std::function<VarFunc_t> fFunc;

    int fID;
    /// The next ID that hasn't yet been assigned
    static int fgNextID;
  };

  /// \brief Representation of a variable to be retrieved from a \ref
  /// caf::StandardRecord object
  ///
  /// A Var consists of a function, taking a StandardRecord and returning the
  /// value of the variable (which may be some complicated function).
  typedef _Var<caf::SRSliceProxy> Var;

  /// \brief Equivalent of \ref Var acting on \ref caf::SRSpill
  typedef _Var<caf::SRSpillProxy> SpillVar;

  /// \brief For Vars where literally all you need is a single CAF variable
  ///
  /// eg Var myVar = SIMPLEVAR(my.var.str);
  /// NB lack of quotes quotes around my.var.str
#define SIMPLEVAR(CAFNAME) Var([](const caf::SRSliceProxy* sr){return sr->CAFNAME;})

#define SIMPLESPILLVAR(CAFNAME) SpillVar([](const caf::SRSpillProxy* sr){return sr->CAFNAME;})

  /// The simplest possible Var, always 1. Used as a default weight.
  const Var kUnweighted([](const caf::SRSliceProxy*){return 1;});

  const SpillVar kSpillUnweighted([](const caf::SRSpillProxy*){return 1;});

  /// \brief Variable formed from two input variables
  ///
  /// The binning of each variable has to be given to allow conversion into a
  /// 1D binned representation.
  template<class T> _Var<T>
  Var2D(const _Var<T>& a, const Binning& binsa,
        const _Var<T>& b, const Binning& binsb);

  /// \brief Variable formed from two input variables
  ///
  /// The binning of each variable has to be given to allow conversion into a
  /// 1D binned representation.
  template<class T> _Var<T>
  Var2D(const _Var<T>& a, int na, double a0, double a1,
        const _Var<T>& b, int nb, double b0, double b1);

  /// This is just like a Var2D, but useful for 3D Spectra
  template<class T> _Var<T>
  Var3D(const _Var<T>& a, const Binning& binsa,
        const _Var<T>& b, const Binning& binsb,
        const _Var<T>& c, const Binning& binsc);

  /// This is just like a Var2D, but useful for 3D Spectra
  template<class T> _Var<T>
  Var3D(const _Var<T>& a, int na, double a0, double a1,
        const _Var<T>& b, int nb, double b0, double b1,
        const _Var<T>& c, int nc, double c0, double c1);

  // TODO - remove once no more legacy callers exist
  #define SpillTruthVar2D Var2D
  #define SpillTruthVar3D Var3D

  /// Use to rescale another variable.
  Var Scaled(const Var& v, double s);

  /// Use to weight events up and down by some factor
  Var Constant(double c);

  /// Use to take sqrt of a var
  Var Sqrt(const Var& v);
} // namespace
