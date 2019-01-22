#include <iostream>
#include <vector>
#include <TH2D.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
//#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ObjectBuilder.h"

namespace util {
  namespace SelectionTool {

    ObjectBuilder::ObjectBuilder() : SelectionBase() {}

    void ObjectBuilder::Initialize(fhicl::ParameterSet *p) {
      // Initialise the objects we are accessing
      if(!p) {
        std::cerr << "No parameters declared" << std::endl;
        //mf::LogWarning("ObjectBuilder") << "No parameters declared";
        return;
      }
      
      fhicl::ParameterSet pconfig = p->get<fhicl::ParameterSet>("objectBuilder");
      
      // Handle labels                                                                
      fGeneratorLabel            = pconfig.get<std::string>("TruthLabel");              
      fGeantLabel                = pconfig.get<std::string>("G4Label");                 
      fPandoraLabel              = pconfig.get<std::string>("PandoraLabel");            
      fRecoTrackLabel            = pconfig.get<std::string>("RecoTrackLabel");          
      fRecoShowerLabel           = pconfig.get<std::string>("RecoShowerLabel");         
      fRecoTrackCalorimetryLabel = pconfig.get<std::string>("RecoTrackCalorimetryLabel");
      fRecoTrackParticleidLabel  = pconfig.get<std::string>("RecoTrackParticleIDLabel");
      fHitLabel                  = pconfig.get<std::string>("HitLabel");                
    } // Initialise function

    void ObjectBuilder::Finalize() {
      // Output our histograms to the ROOT file
      //fOutputFile->cd();
    } // Finalise function

    bool ObjectBuilder::ProcessEvent(const gallery::Event& ev) {
      std::cout << " Truth label : " << fGeneratorLabel << std::endl;
      return true;
    } // ProcessEvent function
  }  // namespace SelectionTool
}  // namespace util

// This line must be included for all selections!
//DECLARE_SBN_PROCESSOR(util::SelectionTool::ObjectBuilder)

