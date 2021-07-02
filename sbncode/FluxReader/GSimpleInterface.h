#ifndef _GSIMPLEINTERFACE_H_
#define _GSIMPLEINTERFACE_H_

#include "FluxInterface.h"

#include "Tools/Flux/GNuMIFlux.h"
#include "Tools/Flux/GSimpleNtpFlux.h"

class TTree;
class TFile;

namespace fluxr {
  class GSimpleInterface : public FluxInterface
    {
    public:
      GSimpleInterface();
      ~GSimpleInterface();

      const Long64_t GetEntries()                    {return fNEntries;};
      const int      GetRun()                        {return fRun;};
      const void     SetRun(int run)                 {fRun = run;};
      const float    GetPOT()                        {return fPOT;};
      const TLorentzVector GetNuPosition()           {return fNuPos;};
      const TLorentzVector GetNuMomentum()           {return fNuMom;};

      void SetRootFile(TFile* rootFileName);
      bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

    private:
      TTree*                        fFluxTree;
      TTree*                        fMetaTree;
      genie::flux::GSimpleNtpEntry* fGSimpleEntry;
      genie::flux::GSimpleNtpNuMI*  fGSimpleNuMI;
      genie::flux::GSimpleNtpAux*   fGSimpleAux;
      genie::flux::GSimpleNtpMeta*  fGSimpleMeta;
      Long64_t                      fNEntries;
      int                           fRun;
      float                         fPOT;
      TLorentzVector                fNuPos;
      TLorentzVector                fNuMom;
  };

}

#endif // _GSIMPLEINTERFACE_H_
