////////////////////////////////////////////////////////////////////////
// Class:       DThetaKinkTree
// Plugin Type: analyzer (art v3_05_01)
// File:        DThetaKinkTree_module.cc
//
// Generated at Tue Sep  8 10:02:52 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "Products/DThetaKink.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"



namespace sbn {
  class DThetaKinkTree;
}


class sbn::DThetaKinkTree : public art::EDAnalyzer {
public:
  explicit DThetaKinkTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DThetaKinkTree(DThetaKinkTree const&) = delete;
  DThetaKinkTree(DThetaKinkTree&&) = delete;
  DThetaKinkTree& operator=(DThetaKinkTree const&) = delete;
  DThetaKinkTree& operator=(DThetaKinkTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // helpers
  void Clear();
  void FillTruth(const recob::PFParticle &particle, const simb::MCParticle &trueParticle);
  void FillKinks(const recob::PFParticle &particle, const std::vector<art::Ptr<sbn::DThetaKink>> &kinks);
  void FillMeta(const art::Event &evt, unsigned i_part);

  // config
  std::vector<std::string> fPFParticleTags;
  std::vector<std::string> fKinkTags;
  float fAngleCut;

  // output
  TTree *_tree;

  // List of kinks per-plane
  std::vector<sbn::DThetaKink> fDThetaU;
  std::vector<sbn::DThetaKink> fDThetaV;
  std::vector<sbn::DThetaKink> fDThetaY;

  std::vector<float> fTrajX;
  std::vector<float> fTrajY;
  std::vector<float> fTrajZ;

  std::vector<float> fTrajPU;
  std::vector<float> fTrajPV;
  std::vector<float> fTrajPY;
  std::vector<float> fTrajPT;

  std::vector<float> fScatterX;
  std::vector<float> fScatterY;
  std::vector<float> fScatterZ;

  std::vector<float> fScatterPU;
  std::vector<float> fScatterPV;
  std::vector<float> fScatterPY;
  std::vector<float> fScatterPT;

  int fEvt;
  int fNo;
};


sbn::DThetaKinkTree::DThetaKinkTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
    fPFParticleTags(p.get<std::vector<std::string>>("ParticleTags")),
    fKinkTags(p.get<std::vector<std::string>>("KinkTags")),
    fAngleCut(p.get<float>("AngleCut", 5.))
  // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  _tree = tfs->make<TTree>("DThetaKinkAnalyzer", "kink_tree");

  _tree->Branch("dtheta_u", &fDThetaU);
  _tree->Branch("dtheta_v", &fDThetaV);
  _tree->Branch("dtheta_y", &fDThetaY);

  _tree->Branch("traj_x", &fTrajX);
  _tree->Branch("traj_y", &fTrajY);
  _tree->Branch("traj_z", &fTrajZ);
  _tree->Branch("traj_pu", &fTrajPU);
  _tree->Branch("traj_pv", &fTrajPV);
  _tree->Branch("traj_py", &fTrajPY);
  _tree->Branch("traj_pt", &fTrajPT);

  _tree->Branch("scatter_x", &fScatterX);
  _tree->Branch("scatter_y", &fScatterY);
  _tree->Branch("scatter_z", &fScatterZ);
  _tree->Branch("scatter_pu", &fScatterPU);
  _tree->Branch("scatter_pv", &fScatterPV);
  _tree->Branch("scatter_py", &fScatterPY);
  _tree->Branch("scatter_pt", &fScatterPT);

  _tree->Branch("evt", &fEvt, "evt/i");
  _tree->Branch("no", &fNo, "no/i");
}

void sbn::DThetaKinkTree::Clear() {
  fDThetaU.clear();
  fDThetaV.clear();
  fDThetaY.clear();

  fTrajX.clear();
  fTrajY.clear();
  fTrajZ.clear();

  fTrajPU.clear();
  fTrajPV.clear();
  fTrajPY.clear();
  fTrajPT.clear();

  fScatterX.clear();
  fScatterY.clear();
  fScatterZ.clear();

  fScatterPU.clear();
  fScatterPV.clear();
  fScatterPY.clear();
  fScatterPT.clear();

  fEvt = 0;
  fNo = 0;
}

void sbn::DThetaKinkTree::FillMeta(const art::Event &evt, unsigned i_part) {
  fNo = i_part;
  fEvt = evt.event();
}

void sbn::DThetaKinkTree::FillTruth(const recob::PFParticle &particle, const simb::MCParticle &trueParticle) {
  (void) particle;

  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();

  // save the particle trajectory
  for (unsigned i_traj = 0; i_traj < trueParticle.NumberTrajectoryPoints(); i_traj ++) {
    TLorentzVector pos = trueParticle.Position(i_traj);

    fTrajX.push_back(pos.X());
    fTrajY.push_back(pos.Y());
    fTrajZ.push_back(pos.Z());

    fTrajPU.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 0, 0, 0) * geo->WirePitch());
    fTrajPV.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 1, 0, 0) * geo->WirePitch());
    fTrajPY.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 2, 0, 0) * geo->WirePitch());
    fTrajPT.push_back(pos.X());

    // also find any elastic scatters
    if (i_traj > 1) {
      // angle in degree
      float angle = trueParticle.Momentum(i_traj).Vect().Angle(trueParticle.Momentum(i_traj-1).Vect()) * 180. / M_PI;
      if (angle > fAngleCut) {
	fScatterX.push_back(pos.X());
	fScatterY.push_back(pos.Y());
	fScatterZ.push_back(pos.Z());
	
	fScatterPU.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 0, 0, 0) * geo->WirePitch());
	fScatterPV.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 1, 0, 0) * geo->WirePitch());
	fScatterPY.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 2, 0, 0) * geo->WirePitch());
	fScatterPT.push_back(pos.X());
      }
    }

  }
}

void sbn::DThetaKinkTree::FillKinks(const recob::PFParticle &particle, const std::vector<art::Ptr<sbn::DThetaKink>> &kinks) {

  for (unsigned i_kink = 0; i_kink < kinks.size(); i_kink++) {
    const sbn::DThetaKink &kink = *kinks[i_kink];
    // std::cout << "Got new kink! Plane: " << kink.wire.Plane << std::endl;
    if (kink.wire.Plane == 0) {
      fDThetaU.push_back(kink);
    }
    else if (kink.wire.Plane == 1) {
      fDThetaV.push_back(kink);
    }
    else if (kink.wire.Plane == 2) {
      fDThetaY.push_back(kink);
    }
  }
}

void sbn::DThetaKinkTree::analyze(art::Event const& evt)
{
  // input data
  std::vector<std::vector<art::Ptr<recob::PFParticle>>> particleList;
  std::vector<art::FindManyP<sbn::DThetaKink>> particleKinkList;
  std::vector<art::FindManyP<recob::Cluster>> particleClusterList;
  
  for (unsigned i = 0; i < fPFParticleTags.size(); i++) {
    art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
    evt.getByLabel(fPFParticleTags.at(i), pfparticle_handle);
    particleList.emplace_back();
    art::fill_ptr_vector(particleList.back(), pfparticle_handle);
    particleKinkList.emplace_back(particleList.back(), evt, fKinkTags.at(i));
    particleClusterList.emplace_back(particleList.back(), evt, fPFParticleTags.at(i));
  }

  // flatten
  std::vector<art::Ptr<recob::PFParticle>> particles;
  std::vector<std::vector<art::Ptr<sbn::DThetaKink>>> particleKinks;
  std::vector<std::vector<art::Ptr<recob::Cluster>>> particleClusters;
  std::vector<std::vector<art::Ptr<recob::Hit>>> particleHits;

  for (unsigned i = 0; i < particleList.size(); i++) {
    particles.insert(particles.end(), particleList[i].begin(), particleList[i].end());

    for (unsigned j = 0; j < particleList[i].size(); j++) {
      particleKinks.push_back(particleKinkList[i].at(j));
      particleClusters.push_back(particleClusterList[i].at(j));

      particleHits.emplace_back();
      art::FindManyP<recob::Hit> clusterHits(particleClusters.back(), evt, fPFParticleTags.at(i));
      for (unsigned i_clus = 0; i_clus < particleClusters.back().size(); i_clus++) {
        particleHits.back().insert(particleHits.back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
      }
    }
  }

  // also get truth stuff
  art::Handle<std::vector<simb::MCParticle>> trueParticleHandle;
  evt.getByLabel("largeant", trueParticleHandle);
  std::vector<art::Ptr<simb::MCParticle>> trueParticles;
  art::fill_ptr_vector(trueParticles, trueParticleHandle);

  // process
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    const recob::PFParticle &particle = *particles[i_part];
    const std::vector<art::Ptr<sbn::DThetaKink>> &kinks = particleKinks[i_part];
    const std::vector<art::Ptr<recob::Hit>> &hits = particleHits[i_part];

    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(hits, true);
    int match_id = matches.size() ? std::max_element(matches.begin(), matches.end(),
                                                     [](const auto &a, const auto &b) { return a.second < b.second; })->first : -1;
    // float matchE = matches.size() ? std::max_element(matches.begin(), matches.end(),
    //                                                  [](const auto &a, const auto &b) { return a.second < b.second; })->second : -1;
    art::Ptr<simb::MCParticle> matched_particle;
    for (art::Ptr<simb::MCParticle> p: trueParticles) {
      if (p->TrackId() == match_id) {
        matched_particle = p;
        break;
      }
    }

    // ignore cases with no match
    if (!matched_particle) continue;

    // fill stuff
    Clear();
    FillTruth(particle, *matched_particle);
    FillKinks(particle, kinks);
    FillMeta(evt, i_part);
    _tree->Fill();
  }
}

DEFINE_ART_MODULE(sbn::DThetaKinkTree)
