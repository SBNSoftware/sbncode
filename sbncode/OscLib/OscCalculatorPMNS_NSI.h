#ifndef OSC_OSCCALCULATORPMNS_NSI_H
#define OSC_OSCCALCULATORPMNS_NSI_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalculatorPMNS_NSI.h                                        //
//                                                                      //
// Adapt the PMNS_NSI calculator to standard interface                  //
// <bckhouse@caltech.edu>						//
// Modifications by MAAceroO <marioacero@mail.uniatlantico.edu.co>      //
//////////////////////////////////////////////////////////////////////////

#include "IOscCalculator.h"
#include "PMNS_NSI.h"

namespace osc
{
  /// \brief Optimized version of \ref OscCalculatorPMNS
  ///
  /// Adapt the \ref PMNS_NSI calculator to standard interface
  class OscCalculatorPMNS_NSI: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorPMNS_NSI();
    virtual ~OscCalculatorPMNS_NSI();

    IOscCalculatorAdjustable* Copy() const override;

    virtual double P(int flavBefore, int flavAfter, double E) override;

    // Standard oscillation parameters
    virtual void SetL     (double L     ) override {fPropDirty = true; fL      = L;}
    virtual void SetRho   (double rho   ) override {fPropDirty = true; fRho    = rho;}
    virtual void SetDmsq21(double dmsq21) override {fDmDirty   = true; fDmsq21 = dmsq21;}
    virtual void SetDmsq32(double dmsq32) override {fDmDirty   = true; fDmsq32 = dmsq32;}
    virtual void SetTh12  (double th12  ) override {fMixDirty  = true; fTh12   = th12;}
    virtual void SetTh13  (double th13  ) override {fMixDirty  = true; fTh13   = th13;}
    virtual void SetTh23  (double th23  ) override {fMixDirty  = true; fTh23   = th23;}
    virtual void SetdCP   (double dCP   ) override {fMixDirty  = true; fdCP    = dCP;}

    // Non-Standard Interactions parameters (three real -diagonal- and thee complex -off-diagonal: norm and phase-)
    virtual void SetEps_ee      (double eps_ee     ){fEpsDirty = true; fEps_ee      = eps_ee;}
    virtual void SetEps_emu     (double eps_emu    ){fEpsDirty = true; fEps_emu     = eps_emu;}
    virtual void SetEps_etau    (double eps_etau   ){fEpsDirty = true; fEps_etau    = eps_etau;}
    virtual void SetEps_mumu    (double eps_mumu   ){fEpsDirty = true; fEps_mumu    = eps_mumu;}
    virtual void SetEps_mutau   (double eps_mutau  ){fEpsDirty = true; fEps_mutau   = eps_mutau;}
    virtual void SetEps_tautau  (double eps_tautau ){fEpsDirty = true; fEps_tautau  = eps_tautau;}
    virtual void SetDelta_emu   (double Delta_emu  ){fEpsDirty = true; fDelta_emu   = Delta_emu;}
    virtual void SetDelta_etau  (double Delta_etau ){fEpsDirty = true; fDelta_etau  = Delta_etau;}
    virtual void SetDelta_mutau (double Delta_mutau){fEpsDirty = true; fDelta_mutau = Delta_mutau;}

     // Setting the state (oscillation parameters)
    void SetState(std::vector<double> state);

    // Getters
    double GetEps_ee()      const { return fEps_ee; }
    double GetEps_emu()     const { return fEps_emu; }
    double GetEps_etau()    const { return fEps_etau; }
    double GetEps_mumu()    const { return fEps_mumu; }
    double GetEps_mutau()   const { return fEps_mutau; }
    double GetEps_tautau()  const { return fEps_tautau; }
    double GetDelta_emu()   const { return fDelta_emu; }
    double GetDelta_etau()  const { return fDelta_etau; }
    double GetDelta_mutau() const { return fDelta_mutau; }

     // Getting the state (oscillation parameters)
    std::vector<double> GetState() const;

  protected:
    PMNS_NSI fPMNS_NSI;

    bool fMixDirty;
    bool fDmDirty;
    bool fPropDirty;
    bool fEpsDirty;
    double fPrevE;
    int fPrevAnti;

    double fEps_ee;
    double fEps_mumu;
    double fEps_tautau;
    double fEps_emu;
    double fEps_etau;
    double fEps_mutau;
    double fDelta_emu;
    double fDelta_etau;
    double fDelta_mutau;

  };

  const OscCalculatorPMNS_NSI* DowncastToNSI(const IOscCalculator* calc);
  OscCalculatorPMNS_NSI* DowncastToNSI(IOscCalculator* calc);

} // namespace

#endif
