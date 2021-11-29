#ifndef _BOONEINTERFACE_H_
#define _BOONEINTERFACE_H_

#include "FluxInterface.h"

// #include "Tools/Flux/GNuMIFlux.h"
// #include "Tools/Flux/GSimpleNtpFlux.h"

//boone ntuple structure
struct BooNENtuple {
  float beamwgt;          /// Magic Weight : CrossSect * BooBeamNT
  int ntp;                /// 1,2,3,4 : nue, nuebar, numu, numubar
  int npart;              /// number of particles in the chain
                          ///    npart-1 == proton
                          ///    0 == neutrino
  int id[20];             /// id of each particle in before chain, array length 'npart'
  float ini_pos[20][3];   /// 3-pos of particle in before chain, array length 'npart'
  float ini_mom[20][3];   /// 3-mom of particle in before chain, array length 'npart'
  float ini_eng[20];      /// E of particle in before chain, array length 'npart'
  float ini_t[20];        /// "decay" time of particle (starting from proton)
                          ///             in before chain, array length 'npart'
  float fin_mom[20][3];   /// final 3-mom of particle in before chain, array length 'npart'
  float fin_pol[20][3];   /// final 3-polarization of particle in before chain, array length 'npart'

};

//boone beam ntuple structure
struct BeamNtuple {
  float tank_pos_beam[3]; /// 3-position of the flux window;
  float targ_pos_beam[3]; /// Frame conversion from beam to flux frame
  float pot;              /// pot in the file

  float windowbase[3];
  float windowdir1[3];
  float windowdir2[3];
};

class TTree;
class TFile;

namespace fluxr {
  class BooNEInterface : public FluxInterface
    {
    public:
      BooNEInterface();
      ~BooNEInterface();

      const Long64_t GetEntries()                    {return fNEntries;};
      const int      GetRun()                        {return fRun;};
      const void     SetRun(int run)                 {fRun = run;};
      const float    GetPOT()                        {return fPOT;};
      const TLorentzVector GetNuPosition()           {return fNuPos;};
      const TLorentzVector GetNuMomentum()           {return fNuMom;};

      void SetRootFile(TFile* rootFileName);
      bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

    private:
      TTree*         fBNBtree;
      TTree*         fWindowtree;
      BooNENtuple    fBooneNtp;
      BeamNtuple     fBeamNtp;
      double         fMaxWeight;
      double         fMinWeight;
      double         fSumWeights;

      Long64_t       fNEntries;
      int            fRun;
      float          fPOT;
      TLorentzVector fNuPos;
      TLorentzVector fNuMom;
  };

}

#endif // _BOONEINTERFACE_H_
