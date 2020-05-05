/**
 * \file WirePrinter.cc
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

#include "lardataobj/RecoBase/Wire.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class WirePrinter
 * \brief Electron neutrino event selection
 */
class WirePrinter : public core::SelectionBase {
public:
  /** Constructor. */
  WirePrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fTag = config ? config->get<std::string>("Tag", "caldata") : "caldata";
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
    const std::vector<recob::Wire> &wires = *ev.getValidHandle<std::vector<recob::Wire>>(fTag);

    for (const recob::Wire &w: wires) {
      std::cout << "Wire on channel: " << w.Channel() << " view: " << w.View() << std::endl;
      for (auto const &range: w.SignalROI().get_ranges()) {
        std::cout << "New Range: ";
        for (auto const &val: range) {
        std::cout << val << " ";
        }
        std::cout << std::endl;
      }
    }

    return false; 
  }

protected:
  std::string fTag;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::WirePrinter)

