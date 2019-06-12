#include <iostream>
#include <vector>
#include <TH2D.h>
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "AnalysisMain.h"
#include "core/Event.hh"
#include "core/ProviderManager.hh"
#include "util/Interaction.hh"

namespace ana {
  namespace SelectionAnalysis {

    AnalysisMain::AnalysisMain() : 
      fEventCounter(0),
      fTruthEntryCounter(0),
      fRecoEntryCounter(0){}

    void AnalysisMain::Initialize(fhicl::ParameterSet* config) {
      std::cout << " Initialising " << std::endl;
    }

    void AnalysisMain::Finalize() {
      std::cout << "------------------------------------" << std::endl;
      std::cout << " Number of events       : " << fEventCounter      << std::endl;
      std::cout << " Number of truth entries: " << fTruthEntryCounter << std::endl;
      std::cout << " Number of reco entries : " << fRecoEntryCounter  << std::endl;
      fOutputFile->cd();
    }

    bool AnalysisMain::ProcessEvent(
      const gallery::Event& ev,
      const std::vector<Event::Interaction>& truth,
      const std::vector<Event::RecoInteraction>& reco) {

      int start = static_cast<int>(time(NULL));

      if (fEventCounter % 10 == 0)
        std::cout << "AnalysisMain: Processing event " << fEventCounter << std::endl;

      fEventCounter++;
/*
      std::cout << " Number of truth objects : " << truth.size() << std::endl;
      std::cout << " Number of reco objects  : " << reco.size() << std::endl;
      if(truth.size() >= 1){
        for(unsigned int i = 0; i< truth.size(); ++i){
          std::cout << " Neutrino PDG                      : " << truth[i].neutrino.pdg << std::endl;
          std::cout << " Number of true particles          : " << truth[i].finalstate.size() << std::endl;
          std::cout << " Number of reconstructed particles : " << reco[i].recofinalstate.size() << std::endl;
        }
      }
*/
      for(unsigned int i = 0; i < reco.size(); ++i)
        fRecoEntryCounter++;

      for(unsigned int i = 0; i < truth.size(); ++i){
        fTruthEntryCounter++;
        return true;
      }
    }
  }  // namespace SelectionAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::SelectionAnalysis::AnalysisMain)
