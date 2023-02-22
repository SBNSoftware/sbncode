// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "sbncode/CAFMaker/ICAFSelectionTool.h"

#include "TMath.h"

namespace caf {

class CAFSelectionBlindEnergy: public ICAFSelectionTool {
public:

  CAFSelectionBlindEnergy(const fhicl::ParameterSet &p);
  ~CAFSelectionBlindEnergy() {}

  bool Select(StandardRecord &rec) override;

  // override PreScale computation to blind POT
  float GetPrescale() const override;

private:
  double GetBlindPOTScale() const;

  int fPOTBlindSeed;
};

CAFSelectionBlindEnergy::CAFSelectionBlindEnergy(const fhicl::ParameterSet &p):
  ICAFSelectionTool(p),
  fPOTBlindSeed(p.get<int>("POTBlindSeed", 655277))
{}

float CAFSelectionBlindEnergy::GetPrescale() const {
  float prescale_factor = (fPreScale > 0) ? fPreScale : 1.;
  return prescale_factor / GetBlindPOTScale();
}

bool CAFSelectionBlindEnergy::Select(StandardRecord &rec) {
  // Set blinding on
  rec.hdr.isblind = true;

  //simple cuts for trk and shower variables
  //blind events with a potential lepton with momentum > 0.6 that starts in fiducial volume
  for (caf::SRPFP& pfp: rec.reco.pfp) {
    const caf::SRVector3D start = pfp.trk.start;
    if ( ((start.x < -71.1 - 25 && start.x > -369.33 + 25 ) ||
          (start.x > 71.1 + 25 && start.x < 369.33 - 25 )) &&
         (start.y > -181.7 + 25 && start.y < 134.8 - 25 ) &&
         (start.z  > -895.95 + 30 && start.z < 895.95 - 50)) {

      if (pfp.trk.mcsP.fwdP_muon > 0.6) {
        pfp.trk.mcsP.fwdP_muon = TMath::QuietNaN();
      }
      if (pfp.trk.rangeP.p_muon > 0.6) {
        pfp.trk.rangeP.p_muon = TMath::QuietNaN();
      }
    }
  }

  //Note shower energy may not be currently very functional
  for (caf::SRPFP& pfp: rec.reco.pfp) {
    const caf::SRVector3D start = pfp.shw.start;
    if ( ((start.x < -71.1 - 25 && start.x > -369.33 + 25 ) ||
          (start.x > 71.1 + 25 && start.x < 369.33 - 25 )) &&
         (start.y > -181.7 + 25 && start.y < 134.8 - 25 ) &&
         (start.z  > -895.95 + 30 && start.z < 895.95 - 50)) {
      if (pfp.shw.bestplane_energy > 0.6) {
        pfp.shw.bestplane_energy = TMath::QuietNaN();
        pfp.shw.plane[0].energy = TMath::QuietNaN();
        pfp.shw.plane[1].energy = TMath::QuietNaN();
        pfp.shw.plane[2].energy = TMath::QuietNaN();
      }
    }
  }

  // And for slices, check vertex in FV and then check tracks and showers
  for (caf::SRSlice& slc: rec.slc) {
    const caf::SRVector3D vtx = slc.vertex;
    if ( ((vtx.x < -71.1 - 25 && vtx.x > -369.33 + 25 ) ||
          (vtx.x > 71.1 + 25 && vtx.x < 369.33 - 25 )) &&
         (vtx.y > -181.7 + 25 && vtx.y < 134.8 - 25 ) &&
         (vtx.z  > -895.95 + 30 && vtx.z < 895.95 - 50)) {

      for (caf::SRPFP& pfp: slc.reco.pfp) {
        if (pfp.trk.mcsP.fwdP_muon > 0.6) {
          pfp.trk.mcsP.fwdP_muon = TMath::QuietNaN();
        }
        if (pfp.trk.rangeP.p_muon > 0.6) {
          pfp.trk.rangeP.p_muon = TMath::QuietNaN();
        }
      }
      for (caf::SRPFP& pfp: slc.reco.pfp) {
        if (pfp.shw.bestplane_energy > 0.6) {
          pfp.shw.bestplane_energy = TMath::QuietNaN();
          pfp.shw.plane[0].energy = TMath::QuietNaN();
          pfp.shw.plane[1].energy = TMath::QuietNaN();
          pfp.shw.plane[2].energy = TMath::QuietNaN();
        }
      }
    }
  }

  return true;
}

double CAFSelectionBlindEnergy::GetBlindPOTScale() const {
  std::string bstring = std::to_string(fPOTBlindSeed);
  int slen = bstring.length();
  std::string s1 = bstring.substr(0,int(slen/2));
  std::string s2 = bstring.substr(int(slen/2));
  double rat = stod(s1)/stod(s2);
  while (abs(rat)>1){
    rat = -1 * (abs(rat) - 1);
  }
  return 1 + rat*0.3;
}


DEFINE_ART_CLASS_TOOL(CAFSelectionBlindEnergy)

} // end namespace caf
