/**
 * \file OpSimHitPrinter.cc
 *
 *
 * Author:
 */

#include <iostream>
#include <array>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/Experiment.hh"
#include "core/ProviderManager.hh"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "fhiclcpp/ParameterSet.h"
#include "../OpHitFinder/opHitFinderSBND.hh"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class OpSimHitPrinter
 * \brief Electron neutrino event selection
 */
class OpSimHitPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  OpSimHitPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fOpDetHitTag = config ? config->get<std::string>("OpDetHitTag", "ophit") : "ophit";
    MakeOpHits = config ? config->get<bool>("MakeOpHits", false) : false;
    _op_hit_maker = (MakeOpHits) ? new opdet::opHitFinderSBND(config->get<fhicl::ParameterSet>("OpHitMaker"), fProviderManager->GetDetectorClocksProvider()) :
        NULL;
    event_ind = 0;
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
    std::vector<recob::OpHit> op_hits;
    if (MakeOpHits) {
      const std::vector<raw::OpDetWaveform> &op_waveforms = *ev.getValidHandle<std::vector<raw::OpDetWaveform>>("opdaq");
      op_hits = _op_hit_maker->MakeHits(op_waveforms);
    }
    else {
      op_hits = *ev.getValidHandle<std::vector<recob::OpHit>>(fOpDetHitTag);
    }

    std::cout << "New event!" << std::endl;
    for (const recob::OpHit &hit: op_hits) {
      std::cout << "Op hit:" << " channel: " << hit.OpChannel() << " time:" << hit.PeakTime() << std::endl;
    }
    return false; 
  }

protected:
  std::string fOpDetHitTag;
  opdet::opHitFinderSBND *_op_hit_maker;
  bool MakeOpHits;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::OpSimHitPrinter)

