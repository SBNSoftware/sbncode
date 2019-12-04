#ifndef IOSCCALCULATOR_H
#define IOSCCALCULATOR_H

#include "TMD5.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file    IOscCalculator.h                                            //
//                                                                      //
// \brief   General interface to oscillation calculators                //
// \author  Christopher Backhouse - bckhouse@caltech.edu                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

/// Oscillation probability calculators
namespace osc
{
  /// General interface to oscillation calculators
  class IOscCalculator
  {
  public:
    virtual ~IOscCalculator() {}

    virtual IOscCalculator* Copy() const = 0;

    /// E in GeV; flavors as PDG codes (so, neg==>antinu)
    virtual double P(int flavBefore, int flavAfter, double E) = 0;

    /// \brief Use to check two calculators are in the same state
    ///
    /// \return Null means not implemented for this calculator
    virtual TMD5* GetParamsHash() const {return 0;}
  };

  /// Pass neutrinos through unchanged
  class NoOscillations: public IOscCalculator
  {
  public:
    virtual IOscCalculator* Copy() const override {return new NoOscillations;}

    virtual double P(int from, int to, double /*E*/) override
    {
      if(from == to || to == 0) return 1;
      return 0;
    }

    /// Always compare equal to self
    virtual TMD5* GetParamsHash() const override
    {
      TMD5* ret = new TMD5;
      const char* txt = "NoOscillations";
      ret->Update((unsigned char*)txt, strlen(txt));
      ret->Final();
      return ret;
    }
  };

  /// General interface to any calculator that lets you set the parameters
  class IOscCalculatorAdjustable : public IOscCalculator
  {
  public:
    virtual IOscCalculatorAdjustable* Copy() const = 0;

    // These setters are left unimplemented here, since calculators may want
    // to compute additional values when these are set.
    virtual void SetL     (double L     ) = 0;
    virtual void SetRho   (double rho   ) = 0;
    virtual void SetDmsq21(double dmsq21) = 0;
    virtual void SetDmsq32(double dmsq32) = 0;
    virtual void SetTh12  (double th12  ) = 0;
    virtual void SetTh13  (double th13  ) = 0;
    virtual void SetTh23  (double th23  ) = 0;
    virtual void SetdCP   (double dCP   ) = 0;

    virtual double GetL     () const { return fL      ; }
    virtual double GetRho   () const { return fRho    ; }
    virtual double GetDmsq21() const { return fDmsq21 ; }
    virtual double GetDmsq32() const { return fDmsq32 ; }
    virtual double GetTh12  () const { return fTh12   ; }
    virtual double GetTh13  () const { return fTh13   ; }
    virtual double GetTh23  () const { return fTh23   ; }
    virtual double GetdCP   () const { return fdCP    ; }

  protected:
    /// \brief This is only a safe implementation if your calculator only
    /// depends on these eight parameters
    ///
    /// \param txt A string to uniquely identify your calculator class
    TMD5* GetParamsHashDefault(const std::string& txt) const;

    // Set by the user. Generally useful to derived classes
    double fRho; // density (g/cm^3)
    double fL; // baseline (km)
    double fDmsq21;
    double fDmsq32;
    double fTh12;
    double fTh13;
    double fTh23;
    double fdCP;
  };

} // namespace

#endif
