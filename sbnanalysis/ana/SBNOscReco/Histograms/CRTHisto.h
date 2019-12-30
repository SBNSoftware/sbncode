#ifndef _sbnanalysis_CRTHisto_hh
#define _sbnanalysis_CRTHisto_hh

#include <vector>
#include <map>

#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

#include "../Data/RecoEvent.h"
#include "HistoList.h"
#include "../MultiThread/THShared.h"

#include "TFile.h"

namespace ana {
 namespace SBNOsc {
/**
 * Histograms to be filled per track
 */
struct CRTHistos : public HistoList {

  /**
 * Initialize this set of histograms
 * \param postfix The postfix to add to all histogram names
 */
  void Initialize(const std::string &postfix, const std::vector<double> &tagger_volume);

  /**
 * Fill all of the histograms in this class with a track
 * \param track The track to fill
 * \param true_tracks The list of true particles in this event
 */
  void Fill(const numu::CRTHit &hit);

  void Fill(const sbnd::crt::CRTHit &hit);
  void Get(TFile &f, const std::string &postfix);

public:
  TH2Shared crt_hits_xy;
  TH2Shared crt_hits_xz;
  TH2Shared crt_hits_yz;
};
  } // namespace SBNOSc
} // namespace ana

#endif

