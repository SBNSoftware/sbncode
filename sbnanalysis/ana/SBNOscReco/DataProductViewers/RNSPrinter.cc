/**
 * \file SimPhotonPrinter.cc
 *
 *
 * Author:
 */

#include <iostream>
#include <vector>
#include <string>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/Experiment.hh"
#include "core/ProviderManager.hh"

#include "canvas/Persistency/Common/RNGsnapshot.h"


namespace ana {
  namespace SBNOsc {

/**
 * \class SimPhotonPrinter
 * \brief Electron neutrino event selection
 */
class SimPhotonPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  SimPhotonPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
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
    std::cout << "NEW EVENT\n";
    std::vector<gallery::Handle<std::vector<art::RNGsnapshot>>> rngs;
    ev.getManyByType(rngs);
    for (const gallery::Handle<std::vector<art::RNGsnapshot>> &r: rngs) {
      std::cout << "NEW HANDLE:\n";
      for (const art::RNGsnapshot &s: *r) {
        std::cout << "RNG: " << s.label() << " engine: " << s.ekind() << " state: ";
        for (unsigned u: s.state()) std::cout << u << " ";
        std::cout << std::endl;
      }
    }
    return false;
  }

protected:
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::SimPhotonPrinter)

