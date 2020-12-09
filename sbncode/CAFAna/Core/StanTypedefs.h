///  StanTypedefs.h:
///    Typedefs of various types templated over stan::math::var,
///    centralized here for convenience.
///    They don't get put in the header files associated with the types themselves
///    because there's an annoying chain of declarations needed for the typedef'ing,
///    and it's much easier to maintain if it's in a single place.
#pragma once

namespace stan
{
  namespace math
  {
    class var;
  }
}

namespace osc
{
  template <typename T> class _IOscCalcAdjustable;
  typedef _IOscCalcAdjustable<stan::math::var> IOscCalcAdjustableStan;

  template <typename T> class _IOscCalc;
  typedef _IOscCalc<stan::math::var> IOscCalcStan;

  template <typename T> class _OscCalcDMP;
  typedef _OscCalcDMP<stan::math::var> OscCalcDMPStan;

  template <typename T> class _OscCalcPMNS;
  typedef _OscCalcPMNS<stan::math::var> OscCalcPMNSStan;

  template <typename T> class _OscCalcPMNSOpt;
  typedef _OscCalcPMNSOpt<stan::math::var> OscCalcPMNSOptStan;

  namespace analytic{template <typename T> class _OscCalc;}
  typedef analytic::_OscCalc<stan::math::var> OscCalcAnalyticStan;
}

namespace ana
{
  // note: typedefs over forward-declared types are fragile.
  // if either the underlying type changes, its forward declaration
  // needs to be updated here.

  // ---------------------
  // vars
  template <typename T> class _IFitVar;
  typedef _IFitVar<stan::math::var> IFitVarStan;

  template <typename T> class _IConstrainedFitVar;
  typedef _IConstrainedFitVar<stan::math::var> IConstrainedFitVarStan;

}
