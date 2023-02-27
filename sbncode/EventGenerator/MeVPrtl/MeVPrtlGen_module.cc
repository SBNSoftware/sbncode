////////////////////////////////////////////////////////////////////////
// Class:       MeVPrtlGen
// Plugin Type: producer (art v3_02_06)
// File:        MeVPrtlGen_module.cc
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

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlDecay.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MesonParent.h"

#include "Tools/IMesonGen.h"
#include "Tools/IMeVPrtlFlux.h"
#include "Tools/IRayTrace.h"
#include "Tools/IMeVPrtlDecay.h"

#include "TTree.h"

#include <memory>
#include <chrono>

namespace evgen {
  namespace ldm {
    class MeVPrtlGen;
  }
}


class evgen::ldm::MeVPrtlGen : public art::EDProducer {
public:
  explicit MeVPrtlGen(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MeVPrtlGen(MeVPrtlGen const&) = delete;
  MeVPrtlGen(MeVPrtlGen&&) = delete;
  MeVPrtlGen& operator=(MeVPrtlGen const&) = delete;
  MeVPrtlGen& operator=(MeVPrtlGen&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // produce per run and per subrun stuff
  void beginRun(art::Run& run) override ;
  void endSubRun(art::SubRun& sr) override ;

  bool Deweight(double &weight, double &max_weight);

  ~MeVPrtlGen() noexcept {
    std::cout << "GenTool called (" << fNCalls[0] << ") times. Total duration (" << fNTime[0] << ") ms. Duration per call (" << (fNTime[0] / fNCalls[0]) << ") ms.\n";
    std::cout << "FlxTool called (" << fNCalls[1] << ") times. Total duration (" << fNTime[1] << ") ms. Duration per call (" << (fNTime[1] / fNCalls[1]) << ") ms.\n";
    std::cout << "RayTool called (" << fNCalls[2] << ") times. Total duration (" << fNTime[2] << ") ms. Duration per call (" << (fNTime[2] / fNCalls[2]) << ") ms.\n";
    std::cout << "DcyTool called (" << fNCalls[3] << ") times. Total duration (" << fNTime[3] << ") ms. Duration per call (" << (fNTime[3] / fNCalls[3]) << ") ms.\n";

    if (fEngine) delete fEngine;
    if (fMeVPrtl) delete fMeVPrtl;
  }

private:
  bool fProduce;
  bool fAnaOutput;
  bool fVerbose;

  bool fDoDeweight;
  double fSubRunPOT;

  std::unique_ptr<evgen::ldm::IMesonGen> fGenTool;
  std::unique_ptr<evgen::ldm::IMeVPrtlFlux> fFluxTool;
  std::unique_ptr<evgen::ldm::IRayTrace> fRayTool;
  std::unique_ptr<evgen::ldm::IMeVPrtlDecay> fDecayTool;

  double fGenMaxWeight;
  double fFluxMaxWeight;
  double fRayMaxWeight;
  double fDecayMaxWeight;
  double fRayDecayMaxWeight;

  TTree *fTree;
  MeVPrtlTruth *fMeVPrtl;

  CLHEP::HepJamesRandom *fEngine;

  std::array<uint64_t, 4> fNCalls;
  std::array<double, 4> fNTime;
};

evgen::ldm::MeVPrtlGen::MeVPrtlGen(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  // nullify pointers
  fMeVPrtl = NULL;
  fEngine = NULL;

  fProduce = p.get<bool>("Produce", true);
  fAnaOutput = p.get<bool>("AnaOutput", false);
  fVerbose = p.get<bool>("Verbose", true);

  fDoDeweight = p.get<bool>("Deweight", false);
  fSubRunPOT = 0.;

  // Update constants
  if (p.has_key("Constants")) Constants::Configure(p.get<fhicl::ParameterSet>("Constants"));

  // bring in the tools
  fGenTool = art::make_tool<IMesonGen>(p.get<fhicl::ParameterSet>("MesonGen"));
  fFluxTool = art::make_tool<IMeVPrtlFlux>(p.get<fhicl::ParameterSet>("Flux"));
  fRayTool = art::make_tool<IRayTrace>(p.get<fhicl::ParameterSet>("RayTrace"));
  fDecayTool = art::make_tool<IMeVPrtlDecay>(p.get<fhicl::ParameterSet>("Decay"));

  fGenMaxWeight = fGenTool->MaxWeight();
  if (fVerbose) std::cout << "Gen max weight: " << fGenMaxWeight << std::endl;

  fFluxMaxWeight = fFluxTool->MaxWeight();
  if (fVerbose) std::cout << "Flux max weight: " << fFluxMaxWeight << std::endl;

  fRayMaxWeight = fRayTool->MaxWeight();
  if (fVerbose) std::cout << "Ray max weight: " << fRayMaxWeight << std::endl;

  fDecayMaxWeight = fDecayTool->MaxWeight();
  if (fVerbose) std::cout << "Decay max weight: " << fDecayMaxWeight << std::endl;

  fRayDecayMaxWeight = fDecayMaxWeight*fRayMaxWeight;

  if (fProduce) {
    // All the standard generator outputs
    produces< std::vector<simb::MCTruth> >();
    produces< std::vector<simb::MCFlux>  >();
    produces< sumdata::RunData, art::InRun >();
    produces< sumdata::POTSummary, art::InSubRun >();
    produces< art::Assns<simb::MCTruth, simb::MCFlux> >();
    produces< std::vector<sim::BeamGateInfo> >();

    // also save info pertinent to the scalar
    produces< std::vector<evgen::ldm::MeVPrtlTruth> >();
  }

  // setup the random number engine
  art::ServiceHandle<rndm::NuRandomService> seedSvc;
  fEngine = new CLHEP::HepJamesRandom;
  seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), "MeVPrtlGen");

  if (fAnaOutput) {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("mevprtl_gen", "mevprtl_gen");
    fMeVPrtl = new evgen::ldm::MeVPrtlTruth();
    fTree->Branch("mevprtl", &fMeVPrtl);
  }

  fNCalls = {0, 0, 0, 0};
  fNTime = {0., 0., 0., 0.};
}

void evgen::ldm::MeVPrtlGen::beginRun(art::Run& run) {
  art::ServiceHandle<geo::Geometry const> geo;
  if (fProduce) run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
}

void evgen::ldm::MeVPrtlGen::endSubRun(art::SubRun& sr) {
  auto p = std::make_unique<sumdata::POTSummary>();
  p->totpot = fSubRunPOT;
  p->totgoodpot = fSubRunPOT;

  if (fProduce) sr.put(std::move(p));

  fSubRunPOT = 0.;
}

bool evgen::ldm::MeVPrtlGen::Deweight(double &weight, double &max_weight) {
  
  if (!fDoDeweight || max_weight < 0) { //  don't do deweighting procedure
    return true;
  }

  // Guard against bad max weight
  //
  // There is some question of whether to crash or handle this gracefully.
  // Note that there is some literature that proposes that doing rejection
  // sampling with an adaptive maximum weight (as done below) is valid.
  // https://www.jstor.org/stable/4140534?seq=1
  //
  // However, this is really only true for a larger N than what is typically done
  // in a larsoft job (10-150). Still, it's probably best not to crash.
  if (max_weight < weight) {
    std::cerr << "ERROR: weight (" << weight << ") with max weight (" << max_weight << "). Reconfiguration needed.\n";
    std::cout << "ERROR: weight (" << weight << ") with max weight (" << max_weight << "). Reconfiguration needed.\n";
    std::cout << "Updating max_weight to new value!\n";
    max_weight = weight;
    return true;
  }

  // do deweighting
  double rand = CLHEP::RandFlat::shoot(fEngine, 0, max_weight);
  double test = weight;

  // update the weight value
  weight = max_weight;
  
  return rand <= test;
}

void evgen::ldm::MeVPrtlGen::produce(art::Event& evt)
{
  std::unique_ptr<std::vector<simb::MCFlux>> mcfluxColl(new std::vector<simb::MCFlux>);
  std::unique_ptr<std::vector<simb::MCTruth>> mctruthColl(new std::vector<simb::MCTruth>);
  std::unique_ptr<art::Assns<simb::MCTruth, simb::MCFlux>> truth2fluxAssn(new art::Assns<simb::MCTruth, simb::MCFlux>);
  std::unique_ptr<std::vector<sim::BeamGateInfo>> beamgateColl(new std::vector<sim::BeamGateInfo>);

  std::unique_ptr<std::vector<evgen::ldm::MeVPrtlTruth>> mevprtlColl(new std::vector<evgen::ldm::MeVPrtlTruth>);

  // TODO: pileup? For now, don't worry

  // get the next MeVPrtl Truth
  while (1) {

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    simb::MCFlux meson = fGenTool->GetNext();
    fNCalls[0] ++;
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = t2 - t1;
    fNTime[0] = duration.count();

    evgen::ldm::MesonParent mesonp(meson);
    bool is_meson = mesonp.meson_pdg != 0;

    // (void) is_meson;
    if (fVerbose){
      if (is_meson) {
       std::cout << "Flux is meson (" << is_meson << "). Weight: " << mesonp.weight << ". Produced with energy: " << mesonp.mom.E()
               << " M=" << mesonp.mom.M() << " P=(" << mesonp.mom.Px() << ", " << mesonp.mom.Py() << ", " << mesonp.mom.Pz() << ") At: ("
               << mesonp.pos.X() << ", " << mesonp.pos.Y() << ", " << mesonp.pos.Z() << ")" << std::endl;
      }
    }

    bool success;

    evgen::ldm::MeVPrtlFlux flux;
    double flux_weight;

    fNCalls[1] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fFluxTool->MakeFlux(meson, flux, flux_weight) && Deweight(flux_weight, fFluxMaxWeight);
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[1] += duration.count();

    if (!success) continue;
    
    if (fVerbose){
      std::cout << "New flux. E=" << flux.mom.E() << " At: (" << flux.pos.X() << ", " << flux.pos.Y() << ", " << flux.pos.Z() << ")" << std::endl;
      std::cout << "P=(" << flux.mom.Px() << ", " << flux.mom.Py() << ", " << flux.mom.Pz() << ")" << std::endl;
      std::cout << "Flux weight: " << flux_weight << std::endl;
    }

    std::array<TVector3, 2> intersection;
    double ray_weight;

    fNCalls[2] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fRayTool->IntersectDetector(flux, intersection, ray_weight);
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[2] += duration.count();

    if (!success) continue;
    if (fVerbose) std::cout << "Ray weight: " << ray_weight << std::endl;

    evgen::ldm::MeVPrtlDecay decay;
    double decay_weight;

    fNCalls[3] ++;
    t1 = std::chrono::high_resolution_clock::now();
    success = fDecayTool->Decay(flux, intersection[0], intersection[1], decay, decay_weight);
    t2 = std::chrono::high_resolution_clock::now();
    duration = t2 - t1;
    fNTime[3] += duration.count();

    if (!success) continue;

    if (fVerbose) std::cout << "Decay weight: " << decay_weight << std::endl;

    // Deweight the ray and decay weights together because they have some anti-correlation
    double ray_decay_weight = ray_weight * decay_weight;
    success = Deweight(ray_decay_weight, fRayDecayMaxWeight);

    if (!success) continue;

    if (fVerbose) std::cout << "RayDecay weight: " << ray_decay_weight << std::endl;
    if (fVerbose) std::cout << "PASSED!\n";

    // get the POT
    double thisPOT = fGenTool->GetPOT();

    // if we are de-weighting, then the scaling all gets put into the POT variable
    if (fDoDeweight) {
      thisPOT = thisPOT / (flux_weight * ray_decay_weight);
      flux_weight = 1.;
      ray_weight = 1.;
      decay_weight = 1.;
      ray_decay_weight = 1.;
    }

    fSubRunPOT += thisPOT;

    // build the output objects
    evgen::ldm::MeVPrtlTruth mevprtl_truth(flux, decay,
      intersection,
      flux_weight,
      ray_weight,
      decay_weight,
      thisPOT
    );

    mevprtlColl->push_back(mevprtl_truth);

    simb::MCTruth mctruth;

    // Add the "Neutrino" as the 0th MCParticle
    // This hopefully (???) won't do anything too bad and will give us
    // the chance to use the neutrino energy in other things
    simb::MCParticle fakenu(0, meson.fntype, "primary", -1, 0, -1/* don't track */);
    fakenu.AddTrajectoryPoint(mevprtl_truth.decay_pos, TLorentzVector(0, 0, flux.equiv_enu, flux.equiv_enu));
    mctruth.Add(fakenu);
    mctruth.SetNeutrino(-1, -1, -1, -1, -1, -1, 
                        -1., -1., -1., -1.); 

    for (unsigned i_d = 0; i_d < mevprtl_truth.daughter_mom.size(); i_d++) {
      TLorentzVector daughter4p(mevprtl_truth.daughter_mom[i_d], mevprtl_truth.daughter_e[i_d]);
      simb::MCParticle d(0, mevprtl_truth.daughter_pdg[i_d], "primary", -1, daughter4p.M());
      d.AddTrajectoryPoint(mevprtl_truth.decay_pos, daughter4p);
      mctruth.Add(d);
    }

    // TODO:
    //
    // Flux systematic uncertainties are often evaluated on the neutrino energy.
    // We need to figure out how to translate this for the case of heavy particles
    // with different production kinematics. For now, we could save the neutrino energy
    // so that the flux uncertainties "work" at some level.
    //
    // However, the existing MCTruth object has no way to "just" set a neutrino energy.
    // This is __very__ very annoying.

    mctruthColl->push_back(mctruth);
    mcfluxColl->push_back(meson);

    // Make the associations only if we are producing stuff
    // Otherwise this crashes
    if (fProduce) {
      art::PtrMaker<simb::MCFlux> MCFluxPtrMaker {evt};
      art::PtrMaker<simb::MCTruth> MCTruthPtrMaker {evt};

      art::Ptr<simb::MCTruth> MCTruthPtr = MCTruthPtrMaker(mctruthColl->size() - 1);
      art::Ptr<simb::MCFlux> MCFluxPtr = MCFluxPtrMaker(mcfluxColl->size() - 1);
      truth2fluxAssn->addSingle(MCTruthPtr, MCFluxPtr);
    }

    // TODO: implement for real
    sim::BeamGateInfo gate;
    beamgateColl->push_back(gate);

    if (fAnaOutput) {
      *fMeVPrtl = mevprtl_truth;
      fTree->Fill();
    }

    break;
  }

  if (fProduce) {
    evt.put(std::move(mctruthColl));
    evt.put(std::move(mcfluxColl));
    evt.put(std::move(truth2fluxAssn));
    evt.put(std::move(beamgateColl));

    evt.put(std::move(mevprtlColl));
  }

}

DEFINE_ART_MODULE(evgen::ldm::MeVPrtlGen)
