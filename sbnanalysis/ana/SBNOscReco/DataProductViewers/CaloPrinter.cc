/**
 * \file CaloPrinter.cc
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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class CaloPrinter
 * \brief Electron neutrino event selection
 */
class CaloPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  CaloPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fTag = config ? config->get<std::string>("Tag", "pandoraCalo") : "pandoraCalo";
    event_ind = 0;
    z = new TH1D("calo point z", "calo point z", 1000, -1000, 1000);
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
    fOutputFile->cd();
    z->Write();
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
    const std::vector<anab::Calorimetry> calos = *ev.getValidHandle<std::vector<anab::Calorimetry>>(fTag);
    std::cout << "New Event!\n";
    for (const anab::Calorimetry& calo: calos) {
      std::cout << "New Calo!\n";
      for (recob::tracking::Point_t xyz: calo.XYZ()) {
        std::cout << xyz.X() << " " << xyz.Y() << " " << xyz.Z() << std::endl; 
        if (xyz.Z() < -100.) std::cout << "Low!\n";
        z->Fill(xyz.Z());
      }
    }
    return false; 
  }

protected:
  std::string fTag;
  TH1D *z;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::CaloPrinter)

