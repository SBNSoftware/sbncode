#pragma once

// This file defines the basic Cut object. For specific cuts, and examples of
// how to implement your own, see Cuts.h

#include <functional>
#include <set>
#include <string>

#include "CAFAna/Core/Var.h"

namespace caf{class StandardRecord; class SRSpill; class SRSpillTruthBranch; class SRSlice;}

namespace ana
{
  template<class T> class _Cut;

  template<class T> _Cut<T> operator&&(const _Cut<T>& a,
                                             const _Cut<T>& b);
  template<class T> _Cut<T> operator||(const _Cut<T>& a,
                                             const _Cut<T>& b);
  template<class T> _Cut<T> operator!(const _Cut<T>& a);

  typedef double (ExposureFunc_t)(const caf::SRSpill* spill);

  /// Template for Cut and SpillCut
  template<class T> class _Cut
  {
  public:
    /// The type of the function part of a cut
    typedef bool (CutFunc_t)(const T* sr);

    /// std::function can wrap a real function, function object, or lambda
    _Cut(const std::function<CutFunc_t>& func,
               const std::function<ExposureFunc_t>& liveFunc = 0,
               const std::function<ExposureFunc_t>& potFunc = 0);

    /// Allows a cut to be called with bool result = myCut(sr) syntax
    bool operator()(const T* sr) const
    {
      return fFunc(sr);
    }

    /// Provide a Livetime function if your cut is a timing cut etc
    double Livetime(const caf::SRSpill* spill) const
    {
      return fLiveFunc ? fLiveFunc(spill) : -1;
    }

    /// Could be useful for cuts on specific batches?
    double POT(const caf::SRSpill* spill) const
    {
      return fPOTFunc ? fPOTFunc(spill) : -1;
    }

    /// Cuts with the same definition will have the same ID
    int ID() const {return fID;}

    static int MaxID();

  protected:
    friend std::function<ExposureFunc_t> CombineExposures(const std::function<ExposureFunc_t>& a, const std::function<ExposureFunc_t>& b);

    // Give these guys access to the constructor that sets fID.
    friend _Cut<T> operator&&<>(const _Cut<T>& a,
				      const _Cut<T>& b);
    friend _Cut<T> operator||<>(const _Cut<T>& a,
				      const _Cut<T>& b);
    friend _Cut<T> operator!<>(const _Cut<T>& a);
    _Cut(const std::function<CutFunc_t>& fun,
               const std::function<ExposureFunc_t>& liveFunc,
               const std::function<ExposureFunc_t>& potFunc,
               int id)
      : fFunc(fun), fLiveFunc(liveFunc), fPOTFunc(potFunc), fID(id)
    {
    }

    std::function<CutFunc_t> fFunc;
    std::function<ExposureFunc_t> fLiveFunc, fPOTFunc;

    int fID;
    /// The next ID that hasn't yet been assigned
    static int fgNextID;
  };

  /// \brief Representation of a cut (selection) to be applied to a \ref
  /// caf::StandardRecord object
  ///
  /// A Cut consists of a function, taking a StandardRecord and returning a
  /// boolean indicating if that event passes the cut.
  ///
  /// Cut objects may be combined with the standard boolean operations && ||
  /// and !
  typedef _Cut<caf::SRSliceProxy> SliceCut;
  typedef _Cut<caf::SRSliceProxy> Cut;

  /// \brief Equivalent of \ref Cut acting on \ref caf::SRSpill. For use in
  /// spill-by-spill data quality cuts
  typedef _Cut<caf::SRSpillProxy> SpillCut;

  /// \brief Cut designed to be used over the nuTree, ie all neutrinos, not
  /// just those that got slices.
  typedef _Cut<caf::SRSpillTruthBranch> SpillTruthCut;

  template<class T> _Cut<T> operator>(const _Var<T>& v, double c);
  template<class T> _Cut<T> operator<(const _Var<T>& v, double c);
  template<class T> _Cut<T> operator>=(const _Var<T>& v, double c);
  template<class T> _Cut<T> operator<=(const _Var<T>& v, double c);
  template<class T> _Cut<T> operator==(const _Var<T>& v, double c);
  template<class T> _Cut<T> operator!=(const _Var<T>& v, double c);

  template<class T> _Cut<T> operator>(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Cut<T> operator<(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Cut<T> operator>=(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Cut<T> operator<=(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Cut<T> operator==(const _Var<T>& a, const _Var<T>& b);
  template<class T> _Cut<T> operator!=(const _Var<T>& a, const _Var<T>& b);

  template<class T> _Cut<T> operator>(double c, const _Var<T>& v);
  template<class T> _Cut<T> operator<(double c, const _Var<T>& v);
  template<class T> _Cut<T> operator>=(double c, const _Var<T>& v);
  template<class T> _Cut<T> operator<=(double c, const _Var<T>& v);
  template<class T> _Cut<T> operator!=(double c, const _Var<T>& v);

  /// The simplest possible cut: pass everything, used as a default
  const Cut kNoCut([](const caf::SRSliceProxy*){return true;});

  /// The simplest possible cut: pass everything, used as a default
  const SpillCut kNoSpillCut([](const caf::SRSpillProxy*){return true;});
} // namespace
