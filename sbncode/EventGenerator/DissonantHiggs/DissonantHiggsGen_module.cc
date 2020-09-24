////////////////////////////////////////////////////////////////////////
// Class:       DissonantHiggsGen
// Plugin Type: producer (art v3_02_06)
// File:        DissonantHiggsGen_module.cc
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
    class DissonantHiggsGen;
  }
}


class evgen::ldm::DissonantHiggsGen : public art::EDProducer {
public:
  explicit DissonantHiggsGen(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DissonantHiggsGen(DissonantHiggsGen const&) = delete;
  DissonantHiggsGen(DissonantHiggsGen&&) = delete;
  DissonantHiggsGen& operator=(DissonantHiggsGen const&) = delete;
  DissonantHiggsGen& operator=(DissonantHiggsGen&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // produce per run and per subrun stuff
  void beginRun(art::Run& run);
  void endSubRun(art::SubRun& sr);

  bool Deweight(double &weight, double max_weight);

  ~DissonantHiggsGen() noexcept {
    std::cout << "GenTool called (" << fNCalls[0] << ") times. Total duration (" << fNTime[0] << ") ms. Duration per call (" << (fNTime[0] / fNCalls[0]) << ") ms.\n";
    std::cout << "FlxTool called (" << fNCalls[1] << ") times. Total duration (" << fNTime[1] << ") ms. Duration per call (" << (fNTime[1] / fNCalls[1]) << ") ms.\n";
    std::cout << "RayTool called (" << fNCalls[2] << ") times. Total duration (" << fNTime[2] << ") ms. Duration per call (" << (fNTime[2] / fNCalls[2]) << ") ms.\n";
    std::cout << "DcyTool called (" << fNCalls[3] << ") times. Total duration (" << fNTime[3] << ") ms. Duration per call (" << (fNTime[3] / fNCalls[3]) << ") ms.\n";
  }

private:
  bool fDoDeweight;
  double fSubRunPOT;

  std::unique_ptr<evgen::ldm::IKaonGen> fGenTool;
  std::unique_ptr<evgen::ldm::IHiggsFlux> fFluxTool;
  std::unique_ptr<evgen::ldm::IRayTrace> fRayTool;
  std::unique_ptr<evgen::ldm::IHiggsDecay> fDecayTool;
  TTree *fTree;

  CLHEP::HepJamesRandom *fEngine;

  std::array<uint64_t, 4> fNCalls;
  std::array<double, 4> fNTime;
};

evgen::ldm::DissonantHiggsGen::DissonantHiggsGen(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  fDoDeweight = p.get<bool>("Deweight", false);
  fSubRunPOT = 0.;

  // bring in the tools
  fGenTool = art::make_tool<IKaonGen>(p.get<fhicl::ParameterSet>("KaonGen"));
  fFluxTool = art::make_tool<IHiggsFlux>(p.get<fhicl::ParameterSet>("Flux"));
  fRayTool = art::make_tool<IRayTrace>(p.get<fhicl::ParameterSet>("RayTrace"));
  fDecayTool = art::make_tool<IHiggsDecay>(p.get<fhicl::ParameterSet>("Decay"));

  // All the standard generator outputs
  produces< std::vector<simb::MCTruth> >();
  produces< std::vector<simb::MCFlux>  >();
  produces< sumdata::RunData, art::InRun >();
  produces< sumdata::POTSummary, art::InSubRun >();
  produces< art::Assns<simb::MCTruth, simb::MCFlux> >();
  produces< std::vector<sim::BeamGateInfo> >();

  // also save info pertinent to the scalar
  produces< std::vector<evgen::ldm::DissonantHiggs> >();

  // setup the random number engine
  art::ServiceHandle<rndm::NuRandomService> seedSvc;
  fEngine = new CLHEP::HepJamesRandom;
  seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), "DissonantHiggsGen");
    
  fNCalls = {0, 0, 0, 0};
  fNTime = {0., 0., 0., 0.};
}

void evgen::ldm::DissonantHiggsGen::beginRun(art::Run& run) {
  art::ServiceHandle<geo::Geometry const> geo;
  run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
}

void evgen::ldm::DissonantHiggsGen::endSubRun(art::SubRun& sr) {
  auto p = std::make_unique<sumdata::POTSummary>();
  p->totpot = fSubRunPOT;
  p->totgoodpot = fSubRunPOT;

  sr.put(std::move(p));

  fSubRunPOT = 0.;
}

bool evgen::ldm::DissonantHiggsGen::Deweight(double &weight, double max_weight) {
  if (!fDoDeweight || max_weight < 0) { //  don't do deweighting procedure
    return true;
  }

  // guard
  if (max_weight < weight) {
    std::cout << "ERROR: weight (" << weight << ") with max weight (" << max_weight << "). Reconfiguration needed.\n";
  }
  assert(max_weight >= weight);

  // do deweighting
  double rand = CLHEP::RandFlat::shoot(fEngine, 0, max_weight);
  double test = weight;

  // update the weight value
  weight = max_weight;
  return rand <= test;
}

void evgen::ldm::DissonantHiggsGen::produce(art::Event& evt)
{
  std::unique_ptr<std::vector<simb::MCFlux>> mcfluxColl(new std::vector<simb::MCFlux>);
  std::unique_ptr<std::vector<simb::MCTruth>> mctruthColl(new std::vector<simb::MCTruth>);
  std::unique_ptr<art::Assns<simb::MCTruth, simb::MCFlux>> truth2fluxAssn(new art::Assns<simb::MCTruth, simb::MCFlux>);
  std::unique_ptr<std::vector<sim::BeamGateInfo>> beamgateColl(new std::vector<sim::BeamGateInfo>);

  std::unique_ptr<std::vector<evgen::ldm::DissonantHiggs>> higgsColl(new std::vector<evgen::ldm::DissonantHiggs>);

  art::PtrMaker<simb::MCFlux> MCFluxPtrMaker {evt};
  art::PtrMaker<simb::MCTruth> MCTruthPtrMaker {evt};

  // TODO: pileup? For now, don't worry

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

    if (!success) continue;
      
    std::cout << "Ray weight: " << ray_weight << std::endl;

    evgen::ldm::HiggsDecay decay;
    double decay_weight;

    fNCalls[3] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fDecayTool->Decay(higgs, intersection[0], intersection[1], decay, decay_weight) && Deweight(decay_weight, fDecayTool->MaxWeight()); 
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[3] += duration.count();

    if (!success) continue;

    std::cout << "Decay weight: " << decay_weight << std::endl;

    // get the POT
    double thisPOT = fGenTool->GetPOT();

    // if we are de-weighting, then the scaling all gets put into the POT variable
    if (fDoDeweight) {
      thisPOT = thisPOT / (flux_weight * ray_weight * decay_weight);
      flux_weight = 1.;
      ray_weight = 1.;
      decay_weight = 1.;
    }

    fSubRunPOT += thisPOT;

    // build the output objects
    evgen::ldm::DissonantHiggs dhiggs = evgen::ldm::BuildHiggs(higgs, decay, 
      intersection,
      flux_weight,
      ray_weight,
      decay_weight,
      thisPOT
    );

    higgsColl->push_back(dhiggs);

    simb::MCTruth mctruth;

    simb::MCParticle daughterA(0, dhiggs.daughterA_pdg, "primary", -1, dhiggs.daughterA_mom.M());
    daughterA.AddTrajectoryPoint(dhiggs.decay_pos, dhiggs.daughterA_mom);
    simb::MCParticle daughterB(0, dhiggs.daughterB_pdg, "primary", -1, dhiggs.daughterB_mom.M());
    daughterB.AddTrajectoryPoint(dhiggs.decay_pos, dhiggs.daughterB_mom);

    mctruth.Add(daughterA);
    mctruth.Add(daughterB);

    mctruthColl->push_back(mctruth);

    mcfluxColl->push_back(kaon);

    art::Ptr<simb::MCTruth> MCTruthPtr = MCTruthPtrMaker(mctruthColl->size() - 1); 
    art::Ptr<simb::MCFlux> MCFluxPtr = MCFluxPtrMaker(mcfluxColl->size() - 1); 
    truth2fluxAssn->addSingle(MCTruthPtr, MCFluxPtr);

    // TODO: implement for real
    sim::BeamGateInfo gate;
    beamgateColl->push_back(gate);

    break;
  }

  evt.put(std::move(mctruthColl));
  evt.put(std::move(mcfluxColl));
  evt.put(std::move(truth2fluxAssn));
  evt.put(std::move(beamgateColl));

  evt.put(std::move(higgsColl));
}

DEFINE_ART_MODULE(evgen::ldm::DissonantHiggsGen)
