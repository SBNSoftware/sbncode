#include "KaonParent.h"
#include <cassert>

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
      ret.pion_pdg = 212;
      break;
    case 8  /*K-  -> numubar mu-*/:
    case 9  /*K-  -> nuebar pi0 e-*/:
    case 10 /*K-  -> numubar pi0 mu-*/:
      ret.kaon_pdg = -321;
      ret.pion_pdg = -212;
      break;
    default:
      return false; // not a kaon decay
  }

  float branching_ratio = 0.;
  // set the branching ratio
  switch (flux.fndecay) {
    case 1 /*K0L -> nue pi- e+ */:
      branching_ratio = 1.;
      break;
    case 2 /*K0L -> nuebar pi+ e-*/:
      branching_ratio = 1.;
      break;
    case 3 /* K0L -> numu pi- mu+*/:
      branching_ratio = 1.;
      break;
    case 4 /*K0L -> numubar pi+ mu-*/:
      branching_ratio = 1.;
      break;
    case 5  /*K+  -> numu mu+*/:
      branching_ratio = 1.;
      break;
    case 6  /*K+  -> nue pi0 e+*/:
      branching_ratio = 1.;
      break;
    case 7  /*K+  -> numu pi0 mu+*/:
      branching_ratio = 1.;
      break;
    case 8  /*K-  -> numubar mu-*/:
      branching_ratio = 1.;
      break;
    case 9  /*K-  -> nuebar pi0 e-*/:
      branching_ratio = 1.;
      break;
    case 10 /*K-  -> numubar pi0 mu-*/:
      branching_ratio = 1.;
      break;
    default:
      return false; // not a kaon decay
  }

  TVector3 pos3 = TVector3(flux.fvx, flux.fvy, flux.fvz);
  float time = flux.fxpoint; /* README: the MCFlux for some reason does not have any time variable, so I have chosen to canibalize this one, 
                              which according to documentation is just for debugging. (I am very sorry). */

  ret.pos.SetVect(pos3);
  ret.pos.SetT(time);

  ret.mom.SetZ(flux.fpdpz);
  ret.mom.SetX(flux.fpdpy);
  ret.mom.SetY(flux.fpdpx);
  ret.mom.SetE(flux.fppenergy);

  // make sure mass makes sense
  assert(ret.mom.M() > 0.4);

  // set the weight -- divides out the branching ratio
  ret.weight = flux.fnimpwt / branching_ratio;

  return true;
}

