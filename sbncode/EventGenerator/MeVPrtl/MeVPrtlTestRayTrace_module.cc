////////////////////////////////////////////////////////////////////////
// Class:       MeVPrtlTestRayTrace
// Plugin Type: analyzer (art v3_02_06)
// File:        MeVPrtlTestRayTrace_module.cc
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

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlDecay.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MesonParent.h"

#include "Tools/IMesonGen.h"
#include "Tools/IMeVPrtlFlux.h"
#include "Tools/IRayTrace.h"
#include "Tools/IMeVPrtlDecay.h"

#include "Tools/Constants.h"

#include "TTree.h"

#include <memory>
#include <chrono>

namespace evgen {
  namespace ldm {
    class MeVPrtlTestRayTrace;
  }
}


class evgen::ldm::MeVPrtlTestRayTrace : public art::EDAnalyzer {
public:
  explicit MeVPrtlTestRayTrace(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MeVPrtlTestRayTrace(MeVPrtlTestRayTrace const&) = delete;
  MeVPrtlTestRayTrace(MeVPrtlTestRayTrace&&) = delete;
  MeVPrtlTestRayTrace& operator=(MeVPrtlTestRayTrace const&) = delete;
  MeVPrtlTestRayTrace& operator=(MeVPrtlTestRayTrace&&) = delete;

  // Required functions.
  void analyze(const art::Event& e) override;

  ~MeVPrtlTestRayTrace() noexcept {}

private:
  std::unique_ptr<evgen::ldm::IMesonGen> fGenTool;
  std::unique_ptr<evgen::ldm::IMeVPrtlFlux> fFluxTool;
  std::vector<std::unique_ptr<evgen::ldm::IRayTrace>> fRayTools;

  unsigned fNCall;
  unsigned fNEvt;

  TTree *fTree;
  std::vector<double> fBranchWeights;
  MeVPrtlTruth *fMeVPrtl;
};

evgen::ldm::MeVPrtlTestRayTrace::MeVPrtlTestRayTrace(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  // bring in the tools
  fGenTool = art::make_tool<IMesonGen>(p.get<fhicl::ParameterSet>("MesonGen"));
  fFluxTool = art::make_tool<IMeVPrtlFlux>(p.get<fhicl::ParameterSet>("Flux"));

  for (auto const &rayconfig: p.get<std::vector<fhicl::ParameterSet>>("RayTraces")) {
    fRayTools.push_back(art::make_tool<IRayTrace>(rayconfig));
    fBranchWeights.push_back(0);
  }

  fNCall = p.get<unsigned>("NCall", 10000);
  fNEvt = 0;

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("testraytrace", "testraytrace");
  fMeVPrtl = new evgen::ldm::MeVPrtlTruth();
  fTree->Branch("mevprtl", &fMeVPrtl);
  for (unsigned iray = 0; iray < fRayTools.size(); iray++) {
    fTree->Branch((std::string(fRayTools[iray]->Name()) + "_wgt").c_str(), &fBranchWeights[iray]);
  }

}

void evgen::ldm::MeVPrtlTestRayTrace::analyze(const art::Event& evt)
{
  // get the next MeVPrtl Truth 
  while (1) {
    simb::MCFlux meson = fGenTool->GetNext();

    evgen::ldm::MesonParent mesonp(meson);
    bool is_kaon = mesonp.isKaon();

    // (void) is_kaon;
    if (is_kaon) { 
     std::cout << "Flux is kaon (" << is_kaon << "). Weight: " << mesonp.weight << ". Produced with energy: " << mesonp.mom.E() 
             << " M=" << mesonp.mom.M() << " P=(" << mesonp.mom.Px() << ", " << mesonp.mom.Py() << ", " << mesonp.mom.Pz() << ") At: ("
             << mesonp.pos.X() << ", " << mesonp.pos.Y() << ", " << mesonp.pos.Z() << ")" << std::endl;
    }

    bool success;

    evgen::ldm::MeVPrtlFlux flux;
    double flux_weight;

    success = fFluxTool->MakeFlux(meson, flux, flux_weight);
    if (!success) continue;

    std::cout << "New flux. E=" << flux.mom.E() << " At: (" << flux.pos.X() << ", " << flux.pos.Y() << ", " << flux.pos.Z() << ")" << std::endl;
    std::cout << "P=(" << flux.mom.Px() << ", " << flux.mom.Py() << ", " << flux.mom.Pz() << ")" << std::endl;
    std::cout << "Flux weight: " << flux_weight << std::endl;


    // See if an intersection is possible
    double costh_crit = minKinematicCosTheta(flux.mmom.M(), flux.sec.M(), flux.mom.M(), flux.mmom.E());
    TVector3 det(0., 0., 0.); // detector should be near origin
    double costh = flux.mmom.Vect().Unit().Dot((det - flux.pos.Vect()).Unit());
    std::cout << "COSTH CRIT: " << costh_crit << " DETECTOR COSTH: " << costh << std::endl;
    if (costh < costh_crit) continue;

    std::cout << "CALLING RAY TOOLS\n";
    fNEvt ++;
    unsigned iray = 0;
    std::vector<double> weightsum(fRayTools.size(), 0.);
    std::vector<double> time(fRayTools.size(), 0.);
    for (auto const &r: fRayTools) {
      std::array<TVector3, 2> intersection;
      double ray_weight;
      std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

      for (unsigned icall = 0; icall < fNCall; icall++) {
        success = r->IntersectDetector(flux, intersection, ray_weight);
        if (!success) ray_weight = 0.;
        weightsum[iray] += ray_weight;
      }

      std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration = t2 - t1;
      time[iray] = duration.count();

      iray ++;
    }

    std::cout << "DONE CALLING RAY TOOLS\n";
    for (unsigned iray = 0; iray < fRayTools.size(); iray++) {
      std::cout << fRayTools[iray]->Name() << " " << (time[iray] / fNCall) << " " << (weightsum[iray] / fNCall) << std::endl;
    }

    evgen::ldm::MeVPrtlDecay decay;
    std::array<TVector3, 2> intersection;
    *fMeVPrtl = evgen::ldm::MeVPrtlTruth(flux, decay, intersection, flux_weight, 1., 1., 0.);
    for (unsigned iray = 0; iray < fRayTools.size(); iray++) {
      fBranchWeights[iray] = weightsum[iray] / fNCall;
    }

    fTree->Fill();

    break;
  }
}

DEFINE_ART_MODULE(evgen::ldm::MeVPrtlTestRayTrace)
