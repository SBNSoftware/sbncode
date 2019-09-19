#ifndef OSC_OSCCALCULATORSTERILEBEAM_H
#define OSC_OSCCALCULATORSTERILEBEAM_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalculatorSterileBeam.h                                     //
//                                                                      //
// Adapt the PMNS_Sterile calculator to standard interface              //
// <aurisaam@ucmail.uc.edu>                                             //
// <kasettisivaprasad@gmail.com>				     	//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "IOscCalculator.h"
#include "OscCalculatorSterile.h"
//#include "INoOscillations.h"
#include "PMNS_Sterile.h"

namespace osc
{
  /// \brief Adapt the PMNS_Sterile calculator to standard interface
  ///
  /// Adapt the \ref PMNS_Sterile calculator (3+N with matter effects) to standard interface
  class OscCalculatorSterileBeam: public OscCalculatorSterile //, public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorSterileBeam();
    virtual ~OscCalculatorSterileBeam();

    OscCalculatorSterileBeam(const OscCalculatorSterileBeam& calc);

    virtual IOscCalculatorAdjustable* Copy() const override;

    std::string kBeamMode;

    virtual void SetKaonScale(double scale);
    virtual void SetPionScale(double scale);
    virtual void SetMuonScale(double scale);

    double GetKaonScale() const;
    double GetPionScale() const;
    double GetMuonScale() const;

    virtual TMD5* GetParamsHash() const override;

  protected:

    double fKaonscale;
    double fPionscale;
    double fMuonscale;

/*  PMNS_Sterile* fPMNS_Sterile;
    int    fNFlavors;
    double fRho;
    bool   fDirty;
    double fPrevE;
    int    fPrevAnti;
    int    fPrevFlavBefore;
*/

  };

  const OscCalculatorSterileBeam* DowncastToSterileBeam(const IOscCalculator* calc);
  OscCalculatorSterileBeam* DowncastToSterileBeam(IOscCalculator* calc);

} // namespace

#endif
