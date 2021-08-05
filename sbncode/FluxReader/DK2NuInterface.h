#ifndef _DK2NUINTERFACE_H_
#define _DK2NUINTERFACE_H_
#include "FluxInterface.h"
#include <vector>

class TTree;
class TFile;

#include "TLorentzRotation.h"
#include "TRandom3.h"

#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/NuChoice.h"

#include "fhiclcpp/ParameterSet.h"

namespace fluxr {
  class DK2NuInterface : public FluxInterface
    {
    public:
      DK2NuInterface();
      ~DK2NuInterface();

      const Long64_t GetEntries()                    {return fNEntries;};
      const int      GetRun()                        {return fRun;};
      const void     SetRun(int run)                 {fRun = run;};
      const float    GetPOT()                        {return fPOT;};
      const TLorentzVector GetNuPosition()           {return fNuPos;};
      const TLorentzVector GetNuMomentum()           {return fNuMom;};

      void SetRootFile(TFile* rootFile);
      bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

      bsim::Dk2Nu* GetDk2Nu() {return fDk2Nu;};
      bsim::NuChoice* GetNuChoice() {return fNuChoice;};

      void Init(fhicl::ParameterSet const & ps);
      void User2BeamPos(const TLorentzVector& usrxyz, TLorentzVector& beamxyz) const;
      void Beam2UserPos(const TLorentzVector& beamxyz, TLorentzVector& usrxyz) const;
      void Beam2UserP4(const TLorentzVector& beamp4, TLorentzVector& usrp4) const;
      TVector3 AnglesToAxis(double theta, double phi);

    private:
      TTree*                      fDk2NuTree;
      TTree*                      fDkMetaTree;
      bsim::Dk2Nu*                fDk2Nu;
      bsim::DkMeta*               fDkMeta;
      bsim::NuChoice*             fNuChoice;
      Long64_t                    fNEntries;
      int                         fRun;
      float                       fPOT;

      TLorentzVector              fNuPos;
      TLorentzVector              fNuMom;

      TRotation fBeamRotXML, fTempRot;
      TLorentzRotation fBeamRot, fBeamRotInv;
      TVector3 fBeamPosXML;
      TLorentzVector fBeamZero;
      TVector3 fFluxWindowPtUser[3];
      TLorentzVector fFluxWindowBase, fFluxWindowDir1, fFluxWindowDir2;
      TVector3 fWindowNormal;
      Double_t fFluxWindowLen1, fFluxWindowLen2;
      Double_t fWindowArea;
      TRandom3 fRnd;
  };

}

#endif // _DK2NUINTERFACE_H_
