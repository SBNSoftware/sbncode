#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "TestDuringMerge.h"
#include "ExampleTools.h"
#include "core/Event.hh"

namespace ana {
  namespace SelectionAnalysis {

    TestDuringMerge::TestDuringMerge() : ObjectBuilder() {}


    void TestDuringMerge::Initialize(Json::Value* config) {
    } // Initialise


    void TestDuringMerge::Finalize() {
    } // Finalise


    bool TestDuringMerge::ProcessEvent(
        const gallery::Event& ev,
        std::vector<Event::RecoInteraction>& reco) {
      return true;
    } // Process Event
  }  // namespace SelectionAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::SelectionAnalysis::TestDuringMerge)

