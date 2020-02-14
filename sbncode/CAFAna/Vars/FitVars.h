#pragma once

#include "CAFAna/Core/IFitVar.h"

#include <limits>

namespace ana
{
  /// \f$ \theta_{13} \f$
  class FitTheta13: public IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th13";}
    virtual std::string LatexName() const {return "#theta_{13}";}
  };

  /// \f$ \theta_{13} \f$
  const FitTheta13 kFitTheta13 = FitTheta13();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{13} \f$
  class FitSinSq2Theta13: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2th13";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{13}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{13} \f$
  const FitSinSq2Theta13 kFitSinSq2Theta13 = FitSinSq2Theta13();

  //----------------------------------------------------------------------

  /// \f$ \delta_{CP}/\pi \f$
  class FitDeltaInPiUnits: public IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "delta(pi)";}
    virtual std::string LatexName() const {return "#delta / #pi";}
  };

  /// \f$ \delta_{CP}/\pi \f$
  const FitDeltaInPiUnits kFitDeltaInPiUnits = FitDeltaInPiUnits();

  //----------------------------------------------------------------------
  /// \f$ \theta_{13} \f$
  class FitTheta23: public IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th23";}
    virtual std::string LatexName() const {return "#theta_{23}";}
  };

  /// \f$ \theta_{13} \f$
  const FitTheta23 kFitTheta23 = FitTheta23();
  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{23} \f$
  class FitSinSqTheta23: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ssth23";}
    virtual std::string LatexName() const {return "sin^{2}#theta_{23}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{23} \f$
  const FitSinSqTheta23 kFitSinSqTheta23 = FitSinSqTheta23();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{23} \f$
  class FitSinSq2Theta23: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2th23";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{23}";}


    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{23} \f$
  const FitSinSq2Theta23 kFitSinSq2Theta23 = FitSinSq2Theta23();

  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{32} \f$
  class FitDmSq32: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq32";}
    virtual std::string LatexName() const {return "#Deltam^{2}_{32}";}

    // "1eV^2 splitting should be enough for anyone"
    // OscCalculatorPMNS freaks out at large splittings
    virtual double LowLimit() const {return -1;}
    virtual double HighLimit() const {return +1;}
  };

  /// \f$ \Delta m^2_{32} \f$
  const FitDmSq32 kFitDmSq32 = FitDmSq32();

  //-------------------------------------------------------------------------

  /// \f$ \Delta m^2_{32}\times10^3{\rm eV}^2 \f$
  class FitDmSq32Scaled: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq32scaled";}
    virtual std::string LatexName() const {return "#Deltam^{2}_{32} (10^{-3} eV^{2})";}

    // "1eV^2 splitting should be enough for anyone"
    // OscCalculatorPMNS freaks out at large splittings
    virtual double LowLimit() const {return -1000;}
    virtual double HighLimit() const {return +1000;}
  };

  /// \f$ \Delta m^2_{32}\times10^3{\rm eV}^2 \f$
  const FitDmSq32Scaled kFitDmSq32Scaled = FitDmSq32Scaled();

  //----------------------------------------------------------------------

  /// \f$ \tan^2\theta_{12} \f$
  class FitTanSqTheta12: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "tsth12";}
    virtual std::string LatexName() const {return "tan^{2}#theta_{12}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return std::numeric_limits<double>::max();}
  };

  /// \f$ \tan^2\theta_{12} \f$
  const FitTanSqTheta12 kFitTanSqTheta12 = FitTanSqTheta12();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{12} \f$
  class FitSinSq2Theta12: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2th12";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{12}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{12} \f$
  const FitSinSq2Theta12 kFitSinSq2Theta12 = FitSinSq2Theta12();

  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{21} \f$
  class FitDmSq21: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq21";}
    virtual std::string LatexName() const {return "#Deltam^{2}_{21}";}

    // "1eV^2 splitting should be enough for anyone"
    // OscCalculatorPMNS freaks out at large splittings
    virtual double LowLimit() const {return -1;}
    virtual double HighLimit() const {return +1;}
  };

  /// \f$ \Delta m^2_{21} \f$
  const FitDmSq21 kFitDmSq21 = FitDmSq21();

  /// \f$ \rho \f$
  class FitRho: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "rho";}
    virtual std::string LatexName() const {return "#rho";}

    //Density should be greater than zero (set a ridiculously high high limit)
    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 10.0;}

  };

  /// \f$ \rho \f$
  const FitRho kFitRho = FitRho();

  //----------------------------------------------------------------------



} // namespace
