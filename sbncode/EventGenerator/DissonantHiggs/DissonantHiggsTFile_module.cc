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
#include "Products/DissonantHiggs.h"
#include "Products/HiggsFlux.h"
#include "Products/HiggsDecay.h"
#include "Products/KaonParent.h"

#include "TTree.h"

#include <memory>
#include <chrono>

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

  bool Deweight(double &weight, double max_weight);

  ~DissonantHiggsTFile() noexcept {
    std::cout << "GenTool called (" << fNCalls[0] << ") times. Total duration (" << fNTime[0] << ") ms. Duration per call (" << (fNTime[0] / fNCalls[0]) << ") ms.\n";
    std::cout << "FlxTool called (" << fNCalls[1] << ") times. Total duration (" << fNTime[1] << ") ms. Duration per call (" << (fNTime[1] / fNCalls[1]) << ") ms.\n";
    std::cout << "RayTool called (" << fNCalls[2] << ") times. Total duration (" << fNTime[2] << ") ms. Duration per call (" << (fNTime[2] / fNCalls[2]) << ") ms.\n";
    std::cout << "DcyTool called (" << fNCalls[3] << ") times. Total duration (" << fNTime[3] << ") ms. Duration per call (" << (fNTime[3] / fNCalls[3]) << ") ms.\n";
  }

private:
  bool fDoDeweight;

  std::unique_ptr<evgen::ldm::IKaonGen> fGenTool;
  std::unique_ptr<evgen::ldm::IHiggsFlux> fFluxTool;
  std::unique_ptr<evgen::ldm::IRayTrace> fRayTool;
  std::unique_ptr<evgen::ldm::IHiggsDecay> fDecayTool;
  TTree *fTree;
  DissonantHiggs *fHiggs;

  CLHEP::HepJamesRandom *fEngine;

  std::array<uint64_t, 4> fNCalls;
  std::array<double, 4> fNTime;
};

evgen::ldm::DissonantHiggsTFile::DissonantHiggsTFile(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fDoDeweight = p.get<bool>("Deweight", false);

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
  fHiggs = new evgen::ldm::DissonantHiggs();
  fTree->Branch("higgs", &fHiggs);

  fNCalls = {0, 0, 0, 0};
  fNTime = {0., 0., 0., 0.};
}

bool evgen::ldm::DissonantHiggsTFile::Deweight(double &weight, double max_weight) {
  if (!fDoDeweight || max_weight < 0) { //  don't do deweighting procedure
    return true;
  }
  // do deweighting
  assert(max_weight >= weight);
  double rand = CLHEP::RandFlat::shoot(fEngine, 0, max_weight);
  double test = weight;

  // update the weight value
  weight = max_weight;
  return rand <= test;
}

void evgen::ldm::DissonantHiggsTFile::analyze(const art::Event& e)
{
  // get the next Higgs Scalar  
  while (1) {
    simb::MCFlux kaon = fGenTool->GetNext();

    evgen::ldm::KaonParent kaonp;
    fNCalls[0] ++;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    bool is_kaon = MakeKaonParent(kaon, kaonp);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = t2 - t1;
    fNTime[0] = duration.count();

    // (void) is_kaon;
    if (is_kaon) { 
     std::cout << "Flux is kaon (" << is_kaon << "). Weight: " << kaonp.weight << ". Produced with energy: " << kaonp.mom.E() 
             << " M=" << kaonp.mom.M() << " P=(" << kaonp.mom.Px() << ", " << kaonp.mom.Py() << ", " << kaonp.mom.Pz() << ") At: ("
             << kaonp.pos.X() << ", " << kaonp.pos.Y() << ", " << kaonp.pos.Z() << ")" << std::endl;
    }

    bool success;

    evgen::ldm::HiggsFlux higgs;
    double flux_weight;

    fNCalls[1] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fFluxTool->MakeFlux(kaon, higgs, flux_weight) && Deweight(flux_weight, fFluxTool->MaxWeight());
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[1] += duration.count();

    if (!success) continue;

    std::cout << "New higgs. E=" << higgs.mom.E() << " At: (" << higgs.pos.X() << ", " << higgs.pos.Y() << ", " << higgs.pos.Z() << ")" << std::endl;
    std::cout << "P=(" << higgs.mom.Px() << ", " << higgs.mom.Py() << ", " << higgs.mom.Pz() << ")" << std::endl;
    std::cout << "Flux weight: " << flux_weight << std::endl;

    std::array<TVector3, 2> intersection;
    double ray_weight;

    fNCalls[2] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fRayTool->IntersectDetector(higgs, intersection, ray_weight) && Deweight(ray_weight, fRayTool->MaxWeight());
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[2] += duration.count();

    // if (!success) {
    //   simb::MCTruth truth;
    //   *fHiggs = evgen::ldm::BuildDecay(higgs, kaonp, truth, 
    //		intersection, 0., 0., 0., 0.);
    //   fTree->Fill();
    // }


    if (!success) continue;
      
    std::cout << "Ray weight: " << ray_weight << std::endl;

    simb::MCTruth truth;
    evgen::ldm::HiggsDecay decay;
    double decay_weight;

    fNCalls[3] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fDecayTool->Decay(higgs, intersection[0], intersection[1], decay, truth, decay_weight) && Deweight(decay_weight, fDecayTool->MaxWeight()); 
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[3] += duration.count();

    if (!success) continue;

    std::cout << "Decay weight: " << decay_weight << std::endl;

    *fHiggs = evgen::ldm::BuildHiggs(higgs, kaonp, decay, truth, 
      intersection,
      flux_weight,
      ray_weight,
      decay_weight,
      fGenTool->GetPOT()
    );
    std::cout << "Decays At: (" << fHiggs->decay_pos.X() << ", " << fHiggs->decay_pos.Y() << ", " << fHiggs->decay_pos.Z() << ")" << std::endl;
    std::cout << "Into A E=" << fHiggs->daughterA_mom.E() << " Px=" << fHiggs->daughterA_mom.Px() << std::endl;
    std::cout << "Into B E=" << fHiggs->daughterB_mom.E() << " Px=" << fHiggs->daughterB_mom.Px() << std::endl;


    break;
  }
  fTree->Fill();
}

DEFINE_ART_MODULE(evgen::ldm::DissonantHiggsTFile)
