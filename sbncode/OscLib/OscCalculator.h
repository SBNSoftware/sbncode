#ifndef OSCCALCULATOR_H
#define OSCCALCULATOR_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file   OscCalculator.h                                              //
//                                                                      //
// \brief  Class with methods for calculating all things                //
//         related to oscillation probabilities.                        //
// \author <rbpatter@caltech.edu>					//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalculator.h"

class TF1;

namespace osc
{

  class OscCalculator : public IOscCalculatorAdjustable
  {
  public:
    OscCalculator();
    virtual ~OscCalculator();

    virtual IOscCalculatorAdjustable* Copy() const override;

    /// E in GeV; flavors as PDG codes (so, neg==>antinu)
    double P(int flavBefore, int flavAfter, double E) override;

    double P_me(double E, bool antinu=false);
    double P_mm(double E, bool antinu=false);
    double P_mt(double E, bool antinu=false);
    double P_ee(double E, bool antinu=false);
    double P_em(double E, bool antinu=false);
    double P_et(double E, bool antinu=false);
    double P_te(double E, bool antinu=false);
    double P_tm(double E, bool antinu=false);
    double P_tt(double E, bool antinu=false);

    double P_null(double, bool) { return 0; }

    void SetL     (double L     ) override { fL      = L;      fUpdated = false; }
    void SetRho   (double rho   ) override { fRho    = rho ? rho : 1e-10; fUpdated = false;}
    void SetDmsq21(double dmsq21) override { fDmsq21 = dmsq21; fUpdated = false; }
    void SetDmsq32(double dmsq32) override { fDmsq32 = dmsq32; fUpdated = false; }
    void SetTh12  (double th12  ) override { fTh12   = th12;   fUpdated = false; }
    void SetTh13  (double th13  ) override { fTh13   = th13;   fUpdated = false; }
    void SetTh23  (double th23  ) override { fTh23   = th23;   fUpdated = false; }
    void SetdCP   (double dCP   ) override { fdCP    = dCP;    fUpdated = false; }

    // Get a TF1 for a give channel's P(E).  Reconfigurations of
    // the osc parameters do not require a new TF1.  (The TF1 just
    // accesses the same underlying functions.)  Thus, you need
    // a new OscCalculator if you want two configurations for
    // the same channel.  Having multiple channels (TF1s) from a
    // single OscCalculator is fine, though.
    //
    // NOTE: It's up to you to delete the returned object.
    TF1 *GetTF1(int flavBefore, int flavAfter);

    // You shouldn't call this ever.  It needs to be public, though,
    // so TF1 can access it.
    double P_wrapper(double *x, double *p);

    virtual TMD5* GetParamsHash() const override
    {
      return IOscCalculatorAdjustable::GetParamsHashDefault("OscCalc");
    }

  private:
    // Update derived parameters when required
    void UpdateBasic();
    void UpdateEDep(double E, bool antinu, bool fliptime);

    double P_internal_ee(double E, bool antinu, bool fliptime);
    double P_internal_me(double E, bool antinu, bool fliptime);
    double P_internal_te(double E, bool antinu, bool fliptime);
    double P_internal_mt(double E, bool antinu, bool fliptime);

    // Flags
    bool fUpdated;

    // Calculated from user parameters once (non-E-dependent)
    double fDmsq31;
    double fsin_th12;
    double fsin_th13;
    double fsin_th23;
    double fcos_th12;
    double fcos_th13;
    double fcos_th23;
    double fsin_2th12;
    double fsin_2th13;
    double fsin_2th23;
    double fcos_2th12;
    double fcos_2th13;
    double fcos_2th23;
    double fsin_sq_th12;
    double fsin_sq_th13;
    double fsin_sq_th23;
    double fcos_sq_th12;
    double fcos_sq_th13;
    double fcos_sq_th23;
    double fsin_sq_2th12;
    double fsin_sq_2th13;
    double fsin_sq_2th23;
    double fcos_sq_2th12;
    double fcos_sq_2th13;
    double fcos_sq_2th23;
    double falpha;
    double fV;

    // Calculated from user parameters every time (E-dependent)
    double fA;
    double fD;
    double fC12;
    double fC13;
    double fdCPproxy;
    double fsin_dCPproxy;
    double fcos_dCPproxy;

  };

}

#endif
