#include "FoldEMDaughters.h"  
#include <map>

bool isShowerDaughterProcess(caf::g4_process_ p) {
  return (
    (p == caf::kG4conv) ||
    (p == caf::kG4muPairProd) ||
    (p == caf::kG4hPairProd) ||
    (p == caf::kG4compt) ||
    (p == caf::kG4eBrem) ||
    (p == caf::kG4muBrems) ||
    (p == caf::kG4hBrems) ||
    (p == caf::kG4phot) ||
    (p == caf::kG4photonNuclear) ||
    (p == caf::kG4muIoni) ||
    (p == caf::kG4eIoni) ||
    (p == caf::kG4hIoni) ||
    (p == caf::kG4ionIoni) ||
    (p == caf::kG4annihil));
}

bool isShowerDaughter(const caf::SRTrueParticle &p) {
  return isShowerDaughterProcess(p.start_process);
}

std::vector<caf::SRTrueParticle> caf::FoldEMShowerDaughters(const std::vector<caf::SRTrueParticle> &ps) {
  std::vector<caf::SRTrueParticle> ret;

  // Speed gain -- construct a map from track ID to particle for fast lookup
  std::map<int, const caf::SRTrueParticle *> idmap;
  for (const caf::SRTrueParticle &p: ps) {
    idmap[p.G4ID] = &p;
  }

  // Do the rollup
  for (const caf::SRTrueParticle &p: ps) {
    if (!isShowerDaughter(p)) continue;

    ret.push_back(p);
    
    // Not a Shower-Daughter!
    // For these particles -- collect up the visible energy of the daughters
    for (unsigned d: p.daughters) {
      if (!idmap.count((int)d)) continue;
      const caf::SRTrueParticle &d_p = *idmap.at((int)d);
      if (!isShowerDaughter(d_p)) continue;

      for (unsigned i_c = 0; i_c < 2; i_c++) {
        for (unsigned i_p = 0; i_p < 3; i_p++) {
          ret.back().plane[i_c][i_p].visE += d_p.plane[i_c][i_p].visE;
          ret.back().plane[i_c][i_p].nhit += d_p.plane[i_c][i_p].nhit;
        }
      }
    }
  }

  return ret;
  
}
