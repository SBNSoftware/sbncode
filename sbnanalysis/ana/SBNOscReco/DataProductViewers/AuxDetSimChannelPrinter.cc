/**
 * \file AuxDetSimChannelPrinter.cc
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

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class AuxDetSimChannelPrinter
 * \brief Electron neutrino event selection
 */
class AuxDetSimChannelPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  AuxDetSimChannelPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
     fTag = config ? config->get<std::string>("Tag", "largeant") : "largeant";
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
    const std::vector<sim::AuxDetSimChannel> &channels = *ev.getValidHandle<std::vector<sim::AuxDetSimChannel>>(fTag);
    for (const sim::AuxDetSimChannel &ch: channels) {
      std::cout << "Channel ID: " << ch.AuxDetID() << std::endl;
      std::cout << "Sensative ID: " << ch.AuxDetSensitiveID() << std::endl;
      std::cout << "NHits: " << ch.AuxDetIDEs().size() << std::endl << std::endl;
    }
    return false; 
  }

protected:
  std::string fTag;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::AuxDetSimChannelPrinter)

