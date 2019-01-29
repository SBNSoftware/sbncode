#include <iostream>
#include <vector>
#include <TH2D.h>
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "AnalysisMain.h"
#include "core/Event.hh"

namespace ana {
  namespace SelectionAnalysis {

    AnalysisMain::AnalysisMain() : 
      fEventCounter(0),
      fGeneratorLabel(""){}

    void AnalysisMain::Initialize(fhicl::ParameterSet* config) {
      std::cout << " Initialising " << std::endl;
      if (config) {
        fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("SelectionAnalysis");

        fGeneratorLabel = pconfig.get<std::string>("GeneratorLabel", "");
        std::cout << " Generator label: " << fGeneratorLabel << std::endl;
      }
    }

    void AnalysisMain::Finalize() {
      std::cout << "------------------------------------" << std::endl;
      std::cout << " Number of events: " << fEventCounter << std::endl;
    }

    bool AnalysisMain::ProcessEvent(
      const gallery::Event& ev,
      const std::vector<Event::Interaction>& truth,
      std::vector<Event::RecoInteraction>& reco) {
    
      if (fEventCounter % 10 == 0)
        std::cout << "AnalysisMain: Processing event " << fEventCounter << std::endl;

      fEventCounter++;

      /*
      // Grab a data product from the event
      auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);

      // Iterate through the neutrinos
      for (size_t i=0; i<mctruths.size(); i++) {
        auto const& mctruth = mctruths.at(i);

        // Fill neutrino vertex position histogram
        if (mctruth.NeutrinoSet()) 
          fNuVertexXZHist->Fill(mctruth.GetNeutrino().Nu().Vx(), mctruth.GetNeutrino().Nu().Vz());
      }
      */
      return true;
    }
  }  // namespace SelectionAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::SelectionAnalysis::AnalysisMain)
