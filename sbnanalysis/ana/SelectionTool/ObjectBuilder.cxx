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


    void ObjectBuilder::Initialize(fhicl::ParameterSet const &p) {
      // Handle labels                                                                
      m_generator_label              = p.get<std::string>("TruthLabel");              
      m_geant_label                  = p.get<std::string>("G4Label");                 
      m_pandora_label                = p.get<std::string>("PandoraLabel");            
      m_reco_track_label             = p.get<std::string>("RecoTrackLabel");          
      m_reco_shower_label            = p.get<std::string>("RecoShowerLabel");         
      m_reco_track_calorimetry_label = p.get<std::string>("RecoTrackCalorimetryLabel");
      m_reco_track_particleid_label  = p.get<std::string>("RecoTrackParticleIDLabel");
      m_hit_label                    = p.get<std::string>("HitLabel");                
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

