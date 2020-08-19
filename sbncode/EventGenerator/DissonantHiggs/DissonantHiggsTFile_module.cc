////////////////////////////////////////////////////////////////////////
// Class:       DissonantHiggsTFile
// Plugin Type: analyzer (art v3_02_06)
// File:        DissonantHiggsTFile_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "Tools/IKaonGen.h"
#include "Tools/IHiggsFlux.h"
#include "Tools/IRayTrace.h"
#include "Tools/IHiggsDecay.h"
#include "Products/HiggsFlux.h"
#include "Products/HiggsDecay.h"
#include "Products/KaonParent.h"

#include "TTree.h"

#include <memory>

namespace evgen {
  namespace ldm {
    class DissonantHiggsTFile;
  }
}


class evgen::ldm::DissonantHiggsTFile : public art::EDAnalyzer {
public:
  explicit DissonantHiggsTFile(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DissonantHiggsTFile(DissonantHiggsTFile const&) = delete;
  DissonantHiggsTFile(DissonantHiggsTFile&&) = delete;
  DissonantHiggsTFile& operator=(DissonantHiggsTFile const&) = delete;
  DissonantHiggsTFile& operator=(DissonantHiggsTFile&&) = delete;

  // Required functions.
  void analyze(const art::Event& e) override;

  bool Deweight(float weight, float max_weight);

private:
  std::unique_ptr<evgen::ldm::IKaonGen> fGenTool;
  std::unique_ptr<evgen::ldm::IHiggsFlux> fFluxTool;
  std::unique_ptr<evgen::ldm::IRayTrace> fRayTool;
  std::unique_ptr<evgen::ldm::IHiggsDecay> fDecayTool;
  TTree *fTree;
  HiggsDecay *fDecay;

  CLHEP::HepJamesRandom *fEngine;
};

evgen::ldm::DissonantHiggsTFile::DissonantHiggsTFile(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  // bring in the tools
  fGenTool = art::make_tool<IKaonGen>(p.get<fhicl::ParameterSet>("KaonGen"));
  fFluxTool = art::make_tool<IHiggsFlux>(p.get<fhicl::ParameterSet>("Flux"));
  fRayTool = art::make_tool<IRayTrace>(p.get<fhicl::ParameterSet>("RayTrace"));
  fDecayTool = art::make_tool<IHiggsDecay>(p.get<fhicl::ParameterSet>("Decay"));

  // setup the random number engine
  art::ServiceHandle<rndm::NuRandomService> seedSvc;
  fEngine = new CLHEP::HepJamesRandom;
  seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), "DissonantHiggsTFile");
    
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("dissonant_higgs", "dissonant_higgs");
  fDecay = new evgen::ldm::HiggsDecay();
  fTree->Branch("decay", &fDecay);
}

bool evgen::ldm::DissonantHiggsTFile::Deweight(float weight, float max_weight) {
  assert(max_weight >= weight);
  float rand = CLHEP::RandFlat::shoot(fEngine, 0, max_weight);
  return rand <= weight;
}

void evgen::ldm::DissonantHiggsTFile::analyze(const art::Event& e)
{
  // get the next Higgs Scalar  
  while (1) {
    simb::MCFlux kaon = fGenTool->GetNext();

    evgen::ldm::KaonParent kaonp;
    bool is_kaon = MakeKaonParent(kaon, kaonp);
    std::cout << "Flux is kaon (" << is_kaon << "). Weight: " << kaonp.weight << ". Produced with energy: " << kaonp.mom.E() 
              << " P=(" << kaonp.mom.Px() << ", " << kaonp.mom.Py() << ", " << kaonp.mom.Pz() << ") At: "
              << kaonp.pos.X() << ", " << kaonp.pos.Y() << ", " << kaonp.pos.Z() << ")" << std::endl;

    bool success;

    evgen::ldm::HiggsFlux higgs;
    float flux_weight;
    success = fFluxTool->MakeFlux(kaon, higgs, flux_weight) && Deweight(flux_weight, fFluxTool->MaxWeight());
    if (!success) continue;

    std::cout << "New higgs. E=" << higgs.mom.E() << " At: (" << higgs.pos.X() << ", " << higgs.pos.Y() << ", " << higgs.pos.Z() << ")" << std::endl;
    std::vector<TVector3> intersection;
    float ray_weight;
    success = fRayTool->IntersectDetector(higgs, intersection, ray_weight) && Deweight(ray_weight, fRayTool->MaxWeight());
    if (!success) continue;
      
    std::cout << "Intersects Detector." << std::endl;

    simb::MCTruth truth;
    float decay_weight;
    success = fDecayTool->Decay(higgs, intersection[0], intersection[1], truth, decay_weight) && Deweight(decay_weight, fDecayTool->MaxWeight()); 
    if (!success) continue;

    *fDecay = evgen::ldm::BuildDecay(higgs, truth, fGenTool->GetPOT() / (fFluxTool->ConstantWeight() * fRayTool->ConstantWeight() * fDecayTool->ConstantWeight()));
    std::cout << "Decays At: (" << fDecay->decay.X() << ", " << fDecay->decay.Y() << ", " << fDecay->decay.Z() << ")" << std::endl;
    std::cout << "Into A E=" << fDecay->daughterA_mom.E() << " Px=" << fDecay->daughterA_mom.Px() << std::endl;
    std::cout << "Into B E=" << fDecay->daughterB_mom.E() << " Px=" << fDecay->daughterB_mom.Px() << std::endl;

    break;
  }
  fTree->Fill();
}

DEFINE_ART_MODULE(evgen::ldm::DissonantHiggsTFile)
