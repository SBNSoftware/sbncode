/**
 * \file OpDetBackTrackerPrinter.cc
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
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class OpDetBackTrackerPrinter
 * \brief Electron neutrino event selection
 */
class OpDetBackTrackerPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  OpDetBackTrackerPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fSimPhotonTag = config ? config->get<std::string>("SimPhotonTag", "largeant") : "largeant";
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
    std::cout << "New Event!\n";
    const std::vector<sim::OpDetBacktrackerRecord> &records = *ev.getValidHandle<std::vector<sim::OpDetBacktrackerRecord>>(fSimPhotonTag);
    for (const sim::OpDetBacktrackerRecord &record: records) {
      record.Dump(std::cout, "    ");
    }
    return false; 
  }

protected:
  std::string fSimPhotonTag;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::OpDetBackTrackerPrinter)

