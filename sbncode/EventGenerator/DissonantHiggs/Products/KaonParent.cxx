#include "KaonParent.h"
#include <cassert>
#include "../ParticleData.h"

bool evgen::ldm::MakeKaonParent(const simb::MCFlux &flux, evgen::ldm::KaonParent &ret) {
  // set the particle codes
  switch (flux.fndecay) {
    case 1 /*K0L -> nue pi- e+ */:
    case 2 /*K0L -> nuebar pi+ e-*/:
    case 3 /* K0L -> numu pi- mu+*/:
    case 4 /*K0L -> numubar pi+ mu-*/:
      ret.kaon_pdg = 130;
      ret.pion_pdg = 111;
      break;
    case 5  /*K+  -> numu mu+*/:
    case 6  /*K+  -> nue pi0 e+*/:
    case 7  /*K+  -> numu pi0 mu+*/:
      ret.kaon_pdg = 321;
      ret.pion_pdg = 211;
      break;
    case 8  /*K-  -> numubar mu-*/:
    case 9  /*K-  -> nuebar pi0 e-*/:
    case 10 /*K-  -> numubar pi0 mu-*/:
      ret.kaon_pdg = -321;
      ret.pion_pdg = -211;
      break;
    default:
      return false; // not a kaon decay
  }

  TVector3 pos3 = TVector3(flux.fvx, flux.fvy, flux.fvz);
  float time = flux.fxpoint; /* README: the MCFlux for some reason does not have any time variable, so I have chosen to canibalize this one, 
                              which according to documentation is just for debugging. (I am very sorry). */

  ret.pos.SetVect(pos3);
  ret.pos.SetT(time);

  ret.mom.SetVectM(TVector3(flux.fpdpx, flux.fpdpy, flux.fpdpz), evgen::ldm::PDATA->GetParticle(ret.kaon_pdg)->Mass()); 

  ret.weight = flux.fnimpwt;

  return true;
}

