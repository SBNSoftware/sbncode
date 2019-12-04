#ifndef OSCCALCULATORCPT_H
#define OSCCALCULATORCPT_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file    OscCalculatorCPT.h                                          //
//                                                                      //
// \brief   Oscillation caculator contaning separate calculators for    //
//          neutrinos and anti-neutrinos                                //
// \author  Joseph Lozier - jlozier@caltech.edu                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "OscLib/IOscCalculator.h"
#include "OscLib/OscCalculatorPMNSOpt.h"
#include <map>
#include <utility>

namespace ana{
  class SigmaDelta;
}

namespace osc
{

  /// Simple neutrino/anti-neutrino enum
  enum class ENuSign { kNu, kNuBar };

  // Because of the way Surface calculation is parallelized,
  // some SigmaDelta objects need to store data in the calc.
  // The map to parameters is keyed on the location of the helper object.
  using SDMap = std::map< const ana::SigmaDelta*, std::pair<double, double> >;

  /** Oscillation calculator implementing CPT-asymmetry.  Holds two oscillation
      calculators as data members: one for neutrinos, the other for
      anti-neutrinos.  Interface allows changing oscillation parameter for each
      calculator independently or simultaneously. **/
  class OscCalculatorCPT: public IOscCalculatorAdjustable
  {

  friend class ana::SigmaDelta;

  public:

    /// default, uses OscCalculatorPMNSOpt for sub-calculators
    OscCalculatorCPT();

    /** can pass custom calcs, OscCalculatorCPT will then own them
        and delete them when it is destroyed **/
    OscCalculatorCPT(IOscCalculatorAdjustable* calc,
                     IOscCalculatorAdjustable* barcalc,
                     SDMap sigdel={} );

    ~OscCalculatorCPT() override;

    IOscCalculatorAdjustable* Copy() const override;

    double P(int flavBefore, int flavAfter, double E) override;

    // symmetric setters
    void SetL     (double L     ) override {SetL(L, ENuSign::kNu);
                                            SetL(L, ENuSign::kNuBar);}
    void SetRho   (double rho   ) override {SetRho(rho, ENuSign::kNu);
                                            SetRho(rho, ENuSign::kNuBar);}
    void SetDmsq21(double dmsq21) override {SetDmsq21(dmsq21,ENuSign::kNu);
                                            SetDmsq21(dmsq21,ENuSign::kNuBar);}
    void SetDmsq32(double dmsq32) override {SetDmsq32(dmsq32,ENuSign::kNu);
                                            SetDmsq32(dmsq32,ENuSign::kNuBar);}
    void SetTh12  (double th12  ) override {SetTh12(th12, ENuSign::kNu);
                                            SetTh12(th12, ENuSign::kNuBar);}
    void SetTh13  (double th13  ) override {SetTh13(th13, ENuSign::kNu);
                                            SetTh13(th13, ENuSign::kNuBar);}
    void SetTh23  (double th23  ) override {SetTh23(th23, ENuSign::kNu);
                                            SetTh23(th23, ENuSign::kNuBar);}
    void SetdCP   (double dCP   ) override {SetdCP(dCP, ENuSign::kNu);
                                            SetdCP(dCP, ENuSign::kNuBar);}

    // asymmetric setters
    virtual void SetL     (double, ENuSign);
    virtual void SetRho   (double, ENuSign);
    virtual void SetDmsq21(double, ENuSign);
    virtual void SetDmsq32(double, ENuSign);
    virtual void SetTh12  (double, ENuSign);
    virtual void SetTh13  (double, ENuSign);
    virtual void SetTh23  (double, ENuSign);
    virtual void SetdCP   (double, ENuSign);

    // symmetric getters
    double GetL     () const override ;
    double GetRho   () const override ;
    double GetDmsq21() const override ;
    double GetDmsq32() const override ;
    double GetTh12  () const override ;
    double GetTh13  () const override ;
    double GetTh23  () const override ;
    double GetdCP   () const override ;

    // asymmetric setters
    virtual double GetL     (ENuSign) const ;
    virtual double GetRho   (ENuSign) const ;
    virtual double GetDmsq21(ENuSign) const ;
    virtual double GetDmsq32(ENuSign) const ;
    virtual double GetTh12  (ENuSign) const ;
    virtual double GetTh13  (ENuSign) const ;
    virtual double GetTh23  (ENuSign) const ;
    virtual double GetdCP   (ENuSign) const ;

    TMD5* GetParamsHash() const override;

  protected:

    // one calc for neutrinos, one for anti-neutrinos
    IOscCalculatorAdjustable* fCalc;
    IOscCalculatorAdjustable* fBarCalc;

    mutable SDMap fSigDel;

  };

  //---------------------------------------------------------------------------

  const OscCalculatorCPT* DowncastToCPT(const IOscCalculatorAdjustable* osc);
  OscCalculatorCPT* DowncastToCPT(IOscCalculatorAdjustable* osc);

} //namespace

#endif
