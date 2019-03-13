#pragma once

#include "CAFAna/Core/IFitVar.h"
#include "TMath.h"

namespace ana
{
  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{32} \f$
  class FitDmSq32Sterile: public IFitVar
  {
  public: 
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq32";}
    virtual std::string LatexName() const {return "#Deltam^{2}_{32}";}
  };

  /// \f$ \Delta m^2_{32} \f$
  const FitDmSq32Sterile kFitDmSq32Sterile = FitDmSq32Sterile();
  
  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{41} \f$
  class FitDmSq41Sterile: public IFitVar
  {
  public: 
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq41";}
    virtual std::string LatexName() const {return "#Deltam^{2}_{41}";}
  };

  /// \f$ \Delta m^2_{41} \f$
  const FitDmSq41Sterile kFitDmSq41Sterile = FitDmSq41Sterile();

  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{43} \f$
  class FitDmSq43Sterile: public IFitVar
  {
  public: 
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "dmsq43";}
    virtual std::string LatexName() const {return "#Deltam^{2}_{43}";}
  };
  
  /// \f$ \Delta m^2_{43} \f$
  const FitDmSq43Sterile kFitDmSq43Sterile = FitDmSq43Sterile();
  
  //----------------------------------------------------------------------

  /// \f$ \delta_{13}/\pi \f$
  class FitDelta13InPiUnitsSterile: public IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "delta13(pi)";}
    virtual std::string LatexName() const {return "#delta_{13} / #pi";}
  };

  /// \f$ \delta_{CP}/\pi \f$
  const FitDelta13InPiUnitsSterile kFitDelta13InPiUnitsSterile = FitDelta13InPiUnitsSterile();

  //----------------------------------------------------------------------

  /// \f$ \delta_{13}/\pi \f$
  class FitDelta14InPiUnitsSterile: public IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "delta14(pi)";}
    virtual std::string LatexName() const {return "#delta_{14} / #pi";}
  };

  /// \f$ \delta_{14}/\pi \f$
  const FitDelta14InPiUnitsSterile kFitDelta14InPiUnitsSterile = FitDelta14InPiUnitsSterile();

  //----------------------------------------------------------------------

  /// \f$ \delta_{24}/\pi \f$
  class FitDelta24InPiUnitsSterile: public IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "delta24(pi)";}
    virtual std::string LatexName() const {return "#delta_{24} / #pi";}
  };

  /// \f$ \delta_{24}/\pi \f$
  const FitDelta24InPiUnitsSterile kFitDelta24InPiUnitsSterile = FitDelta24InPiUnitsSterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{13} \f$
  class FitTheta13Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th13";}
    virtual std::string LatexName() const {return "#theta_{13}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{13} \f$
  const FitTheta13Sterile kFitTheta13Sterile = FitTheta13Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{13} \f$
  class FitSinSqTheta13Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ssth13";}
    virtual std::string LatexName() const {return "sin^{2}#theta_{13}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{13} \f$
  const FitSinSqTheta13Sterile kFitSinSqTheta13Sterile = FitSinSqTheta13Sterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{23} \f$
  class FitTheta23Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th23";}
    virtual std::string LatexName() const {return "#theta_{23}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{23} \f$
  const FitTheta23Sterile kFitTheta23Sterile = FitTheta23Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{23} \f$
  class FitSinSqTheta23Sterile: public IConstrainedFitVar
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
  const FitSinSqTheta23Sterile kFitSinSqTheta23Sterile = FitSinSqTheta23Sterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{14} \f$
  class FitTheta14Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th14";}
    virtual std::string LatexName() const {return "#theta_{14}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{14} \f$
  const FitTheta14Sterile kFitTheta14Sterile = FitTheta14Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{14} \f$
  class FitSinSqTheta14Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ssth14";}
    virtual std::string LatexName() const {return "sin^{2}#theta_{14}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{14} \f$
  const FitSinSqTheta14Sterile kFitSinSqTheta14Sterile = FitSinSqTheta14Sterile();


 //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{14} \f$
  class FitSinSq2Theta14Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2th14";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{14}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{14} \f$
  const FitSinSq2Theta14Sterile kFitSinSq2Theta14Sterile = FitSinSq2Theta14Sterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{24} \f$
  class FitTheta24Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th24";}
    virtual std::string LatexName() const {return "#theta_{24}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/4;}
  };

  /// \f$ \theta_{24} \f$
  const FitTheta24Sterile kFitTheta24Sterile = FitTheta24Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{24} \f$
  class FitSinSqTheta24Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ssth24";}
    virtual std::string LatexName() const {return "sin^{2}#theta_{24}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{24} \f$
  const FitSinSqTheta24Sterile kFitSinSqTheta24Sterile = FitSinSqTheta24Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{24} \f$
  class FitSinSq2Theta24Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2th24";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{24}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{24} \f$
  const FitSinSq2Theta24Sterile kFitSinSq2Theta24Sterile = FitSinSq2Theta24Sterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{34} \f$
  class FitTheta34Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th34";}
    virtual std::string LatexName() const {return "#theta_{34}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{34} \f$
  const FitTheta34Sterile kFitTheta34Sterile = FitTheta34Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{34} \f$
  class FitSinSqTheta34Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ssth34";}
    virtual std::string LatexName() const {return "sin^{2}#theta_{34}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{34} \f$
  const FitSinSqTheta34Sterile kFitSinSqTheta34Sterile = FitSinSqTheta34Sterile();

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{34} \f$
  class FitSinSq2Theta34Sterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "ss2th34";}
    virtual std::string LatexName() const {return "sin^{2}2#theta_{34}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{34} \f$
  const FitSinSq2Theta34Sterile kFitSinSq2Theta34Sterile = FitSinSq2Theta34Sterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{13} \f$
  class FitTheta13InDegreesSterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th13(degrees)";}
    virtual std::string LatexName() const {return "#theta_{13}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{13} \f$
  const FitTheta13InDegreesSterile kFitTheta13InDegreesSterile = FitTheta13InDegreesSterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{23} \f$
  class FitTheta23InDegreesSterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th23(degrees)";}
    virtual std::string LatexName() const {return "#theta_{23}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{23} \f$
  const FitTheta23InDegreesSterile kFitTheta23InDegreesSterile = FitTheta23InDegreesSterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{14} \f$
  class FitTheta14InDegreesSterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th14(degrees)";}
    virtual std::string LatexName() const {return "#theta_{14}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{14} \f$
  const FitTheta14InDegreesSterile kFitTheta14InDegreesSterile = FitTheta14InDegreesSterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{24} \f$
  class FitTheta24InDegreesSterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th24(degrees)";}
    virtual std::string LatexName() const {return "#theta_{24}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 45;}
  };

  /// \f$ \theta_{24} \f$
  const FitTheta24InDegreesSterile kFitTheta24InDegreesSterile = FitTheta24InDegreesSterile();

  //----------------------------------------------------------------------

  /// \f$ \theta_{34} \f$
  class FitTheta34InDegreesSterile: public IConstrainedFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalculatorAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalculatorAdjustable* osc, double val) const;
    virtual std::string ShortName() const {return "th34(degrees)";}
    virtual std::string LatexName() const {return "#theta_{34}";}

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{34} \f$
  const FitTheta34InDegreesSterile kFitTheta34InDegreesSterile = FitTheta34InDegreesSterile();

  //----------------------------------------------------------------------


} // namespace
