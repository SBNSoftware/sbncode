#include <iostream>
#include <vector>
#include <TH2D.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ObjectBuilder.h"
#include "ana/ExampleAnalysis/ExampleTools.h"

namespace ana {
  namespace SelectionTool {

    ObjectBuilder::ObjectBuilder() : SelectionBase() {}


    //void ObjectBuilder::Initialize(Json::Value* config) {
    void ObjectBuilder::Initialize(fhicl::ParameterSet const &p) {
      // Initialise the objects we are accessing

      // Handle labels                                                                
      fGeneratorLabel;            = p.get<art::InputTag>("TruthLabel");              
      fGeantLabel;                = p.get<art::InputTag>("G4Label");                 
      fPandoraLabel;              = p.get<art::InputTag>("PandoraLabel");            
      fRecoTrackLabel;            = p.get<art::InputTag>("RecoTrackLabel");          
      fRecoShowerLabel;           = p.get<art::InputTag>("RecoShowerLabel");         
      fRecoTrackCalorimetryLabel; = p.get<art::InputTag>("RecoTrackCalorimetryLabel");
      fRecoTrackParticleidLabel;  = p.get<art::InputTag>("RecoTrackParticleIDLabel");
      fHitLabel;                  = p.get<art::InputTag>("HitLabel");                
    }


    void ObjectBuilder::Finalize() {
      // Output our histograms to the ROOT file
      //fOutputFile->cd();
    }


    bool ObjectBuilder::ProcessEvent(const gallery::Event& ev, std::vector<Event::RecoInteraction>& reco) {

      // Example of using library data
      ana::ExampleAnalysis::hello();

      return true;
    } // ProcessEvent function

  }  // namespace SelectionTool
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::SelectionTool::ObjectBuilder)

