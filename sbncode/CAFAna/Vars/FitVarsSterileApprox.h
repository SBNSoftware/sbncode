#pragma once

#include "CAFAna/Core/IFitVar.h"

#pragma once

namespace ana
{
  /// \f$ \Delta m^2 \f$
  class FitDmSqSterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq";}
    virtual std::string LatexName() const {return "#Deltam^{2} (eV^{2})";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1e6;}
  };

  /// \Delta m^2 \f$
  const FitDmSqSterile kFitDmSqSterile = FitDmSqSterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{\mu\mu} \f$
  class FitSinSq2ThetaMuMu: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2thmm";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{#mu#mu}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{\mu\mu} \f$
  const FitSinSq2ThetaMuMu kFitSinSq2ThetaMuMu = FitSinSq2ThetaMuMu();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{\mu e} \f$
  class FitSinSq2ThetaMuE: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2thme";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{#mue}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{\mu e} \f$
  const FitSinSq2ThetaMuE kFitSinSq2ThetaMuE = FitSinSq2ThetaMuE();

} // namespace
