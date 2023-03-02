// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "sbncode/CAFMaker/ICAFSelectionTool.h"

namespace caf {

class CAFSelectionNtrkNproton: public ICAFSelectionTool {
public:

  CAFSelectionNtrkNproton(const fhicl::ParameterSet &p);
  ~CAFSelectionNtrkNproton() {}

  bool Select(StandardRecord &rec) override;

private:
  unsigned fNTrack;
  unsigned fNProton;
};

CAFSelectionNtrkNproton::CAFSelectionNtrkNproton(const fhicl::ParameterSet &p):
  ICAFSelectionTool(p),
  fNTrack(p.get<unsigned>("NTrack", 3)),
  fNProton(p.get<unsigned>("NProton", 1))
{}

bool CAFSelectionNtrkNproton::Select(StandardRecord &rec) {
  // clear the normal reco branch
  rec.reco.pfp.clear();

  std::vector<caf::SRSlice> outslc;
  // run selection per-slice
  for (const caf::SRSlice& slc: rec.slc) {
    // count stuff
    unsigned n_proton = 0;
    unsigned n_trk_proton = 0;

    for (const caf::SRPFP &pfp: slc.reco.pfp) {
      bool is_track = pfp.trackScore > 0.5;
      const SRTrkChi2PID &pid = pfp.trk.chi2pid[2];
      bool is_proton = (pid.pid_ndof > 0) && (pid.chi2_proton < 60) && (pid.chi2_muon > 35);
      n_trk_proton += (is_track | is_proton);
      n_proton += is_proton;
    }
    if (n_proton >= fNProton && n_trk_proton >= fNTrack && !slc.is_clear_cosmic) {
      outslc.push_back(slc);
    }
  }

  rec.slc = outslc;

  // Save all slices so we normalize correctly
  return true;
}

DEFINE_ART_CLASS_TOOL(CAFSelectionNtrkNproton)

} // end namespace sbn
