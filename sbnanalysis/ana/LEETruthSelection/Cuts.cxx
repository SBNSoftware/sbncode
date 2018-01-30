#include <algorithm>
#include <cassert>
#include "Config.h"
#include "Util.h"
#include "Cuts.h"

namespace ana {
  namespace lee_truth_selection {

int GetNp(std::vector<PIDParticle>& p) {
  return std::count_if(p.begin(), p.end(), PDGTest({2212}));
}


int GetNtracks(std::vector<PIDParticle>& p) {
  return std::count_if(p.begin(), p.end(), PDGTest({2212, 211}));
}


int GetNlep(std::vector<PIDParticle>& p, int lpdg) { 
  return std::count_if(p.begin(), p.end(), PDGTest({lpdg}));
}


bool TestSelection(std::vector<PIDParticle>& p, int lpdg, EventType t) {
  // Count protons, tracks, and the chosen lepton type
  int np = GetNp(p);
  int nl = GetNlep(p, lpdg);
  int n_trk = GetNtracks(p);

  bool pass_1l1p = nl == 1 && np == 1 && nl + np == p.size();
  bool pass_1lnp = nl == 1 && np >= 1 && nl + np == p.size();
  bool pass_1lnt = nl == 1 && n_trk >= 1 && nl + n_trk == p.size();
  bool pass_1l0p = nl == 1 && np == 0 && nl + np == p.size();

  switch(t) {
    case kAny:
      return pass_1l1p || pass_1lnp || pass_1lnt || pass_1l0p; 
    case k0p:
      return pass_1l0p;
    case kNp:
      return pass_1lnp;
    case k1p:
      return pass_1l1p;
    case kNtrk:
      return pass_1lnt;
  }

  std::cerr << "TestSelection: Unknown EventType" << std::endl;
  assert(false);
  return false;
}


bool GoodObject(bool isFromNuVertex, bool isPrimaryProcess, int pdg, float ke) {
  return isFromNuVertex && isPrimaryProcess && KineticEnergyThreshold(pdg, ke);
}


bool KineticEnergyThreshold(int pdg, float ke) {
  if (abs(pdg) ==  13 || abs(pdg) == 11 ||
      abs(pdg) == 211 || abs(pdg) == 22) {
    return ke > 20.0;
  }

  if (pdg == 2212) {
    return ke > 40.0;
  }

  if (pdg == 111) {
    return true;
  }

  // we shouldn't be looking at any other particles
  assert(false);
  return false;
}

  }  // namespace lee_truth_selection
}  // namespace ana

