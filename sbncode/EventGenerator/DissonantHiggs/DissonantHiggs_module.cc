////////////////////////////////////////////////////////////////////////
// Class:       DissonantHiggs
// Plugin Type: producer (art v3_02_06)
// File:        DissonantHiggs_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "Tools/IHiggsFlux.h"
#include "Tools/IRayTrace.h"
#include "Tools/IHiggsDecay.h"
#include "Products/HiggsFlux.h"

#include <memory>

namespace evgen {
  namespace ldm {
    class DissonantHiggs;
  }
}


class evgen::ldm::DissonantHiggs : public art::EDProducer {
public:
  explicit DissonantHiggs(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DissonantHiggs(DissonantHiggs const&) = delete;
  DissonantHiggs(DissonantHiggs&&) = delete;
  DissonantHiggs& operator=(DissonantHiggs const&) = delete;
  DissonantHiggs& operator=(DissonantHiggs&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  std::unique_ptr<evgen::ldm::IHiggsFlux> fFluxTool;
  std::unique_ptr<evgen::ldm::IRayTrace> fRayTool;
  std::unique_ptr<evgen::ldm::IHiggsDecay> fDecayTool;
};


evgen::ldm::DissonantHiggs::DissonantHiggs(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  produces< std::vector<simb::MCTruth> >();
  // produces< std::vector<simb::MCFlux>  >();
  // produces< sumdata::RunData, art::InRun >();
  // produces< sumdata::POTSummary, art::InSubRun >();
  // produces< art::Assns<simb::MCTruth, simb::MCFlux> >();
  // produces< std::vector<sim::BeamGateInfo> >();

  // bring in the tools
  fFluxTool = art::make_tool<IHiggsFlux>(p.get<fhicl::ParameterSet>("Flux"));
  fRayTool = art::make_tool<IRayTrace>(p.get<fhicl::ParameterSet>("RayTrace"));
  fDecayTool = art::make_tool<IHiggsDecay>(p.get<fhicl::ParameterSet>("Decay"));

  // setup the random number engine
  art::ServiceHandle<rndm::NuRandomService> seedSvc;

  CLHEP::HepJamesRandom *decayEngine = new CLHEP::HepJamesRandom;
  seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(decayEngine), "DissonantHiggs_Decay");
  fDecayTool->SetEngine(decayEngine);

}

void evgen::ldm::DissonantHiggs::produce(art::Event& e)
{
  std::unique_ptr<std::vector<simb::MCTruth>> mctruthColl;

  // get the next Higgs Scalar  
  bool produced = false;
  while (!produced) {
    evgen::ldm::HiggsFlux flux = fFluxTool->Sample();
    std::vector<TVector3> intersection = fRayTool->IntersectDetector(flux);
    if (intersection.size() == 2) {
      std::pair<simb::MCTruth, float> thisTruth = fDecayTool->Decay(flux, intersection[0], intersection[1]);
      // only take positive weight decays
      if (thisTruth.second > 0.) {
        mctruthColl->push_back(thisTruth.first);
        produced = true;
      }
    }
  }

  e.put(std::move(mctruthColl));
}

DEFINE_ART_MODULE(evgen::ldm::DissonantHiggs)
