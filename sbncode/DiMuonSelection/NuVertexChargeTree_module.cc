////////////////////////////////////////////////////////////////////////
// Class:       NuVertexChargeTree
// Plugin Type: analyzer (art v3_05_01)
// File:        NuVertexChargeTree_module.cc
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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "sbncode/LArRecoProducer/Products/VertexHit.h"

namespace sbn {
  class NuVertexChargeTree;
}

static const unsigned MAX_PFP = 10;

class sbn::NuVertexChargeTree : public art::EDAnalyzer {
public:
  explicit NuVertexChargeTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuVertexChargeTree(NuVertexChargeTree const&) = delete;
  NuVertexChargeTree(NuVertexChargeTree&&) = delete;
  NuVertexChargeTree& operator=(NuVertexChargeTree const&) = delete;
  NuVertexChargeTree& operator=(NuVertexChargeTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // keep track of file index
  void respondToOpenInputFile(const art::FileBlock& fb) {
    (void) fb;
    fIFile += 1;
  }


private:
  // helpers
  void Clear();
  bool InFV(TVector3 pos);
  void FillMeta(const art::Event &evt);
  void FillNeutrino(const simb::MCTruth &nu, 
   const recob::Vertex &vert,
   const std::vector<std::vector<std::pair<int, float>>> &slice_pfp_matches,
   const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles, 
   const std::vector<const sim::GeneratedParticleInfo *> infos);
  void FillVertexHits(const std::vector<art::Ptr<sbn::VertexHit>> &hits);
  void FillSlice(const std::vector<art::Ptr<recob::PFParticle>> &slice_particles);

  // config
  std::vector<std::string> fPFParticleTags;
  std::vector<std::string> fVertexHitTags;
  std::array<float, 6> fFiducialInset;
  std::vector<geo::BoxBoundedGeo> fFiducialVolumes;
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;

  // output
  TTree *_tree;

  int fIEvt;
  int fIFile;
  int fEvt;

  float fNuE;
  int fNuMode;
  int fNuIntType;

  std::vector<float> fFSPE;
  std::vector<float> fFSPKE;
  std::vector<int> fFSPPID;
  std::vector<bool> fFSPMatched;
  std::vector<bool> fFSPMtachedPrimary;

  float fRecoVertDist;

  std::vector<float> fVertexHitCharge;
  std::vector<float> fVertexHitDist;
  std::vector<int> fVertexHitPlane;
  std::vector<int> fNVertexPFPs;
  std::vector<std::array<int, MAX_PFP>> fVertexHitPFPIDs;
  std::vector<std::array<float, MAX_PFP>> fVertexHitPFPPCvals;

  int fNSliceParticles;
  int fNSlicePrimaryTracks;
  int fNSlicePrimaryShowers;
};

sbn::NuVertexChargeTree::NuVertexChargeTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
    fPFParticleTags(p.get<std::vector<std::string>>("ParticleTags")),
    fVertexHitTags(p.get<std::vector<std::string>>("VertexHitTags")),
    fFiducialInset({p.get<float>("xmin"), p.get<float>("xmax"), p.get<float>("ymin"), p.get<float>("ymax"), p.get<float>("zmin"), p.get<float>("zmax")}) 
{
  fIEvt = 0;
  fIFile = 0;
  art::ServiceHandle<art::TFileService> tfs;

  _tree = tfs->make<TTree>("VertexEAnalyzer", "kink_tree");

  _tree->Branch("ifile", &fIFile, "ifile/i");
  _tree->Branch("ievt", &fIEvt, "ievt/i");
  _tree->Branch("evt", &fEvt, "evt/i");

  _tree->Branch("nue", &fNuE, "nue/F");
  _tree->Branch("numode", &fNuMode, "numode/i");
  _tree->Branch("nuinttype", &fNuIntType, "nuinttype/i");

  _tree->Branch("fsp_e", &fFSPE);
  _tree->Branch("fsp_ke", &fFSPKE);
  _tree->Branch("fsp_pid", &fFSPPID);
  _tree->Branch("fsp_matched", &fFSPMatched);
  _tree->Branch("fsp_matched_primary", &fFSPMtachedPrimary);

  _tree->Branch("reco_vert_dist", &fRecoVertDist, "reco_vert_dist/F");

  _tree->Branch("vertex_hit_charge", &fVertexHitCharge);
  _tree->Branch("vertex_hit_dist", &fVertexHitDist);
  _tree->Branch("vertex_hit_plane", &fVertexHitPlane);
  _tree->Branch("vertex_npfps", &fNVertexPFPs);
  _tree->Branch("vertex_hit_pfp_ids", &fVertexHitPFPIDs);
  _tree->Branch("vertex_hit_pfp_pcvals", &fVertexHitPFPPCvals);

  _tree->Branch("n_slice_particles", &fNSliceParticles, "n_slice_particles/i");
  _tree->Branch("n_slice_primary_tracks", &fNSlicePrimaryTracks, "n_slice_particles/i");
  _tree->Branch("n_slice_primary_showers", &fNSlicePrimaryShowers, "n_slice_particles/i");

  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes 
  for (auto const &cryo: geometry->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator iTPC = geometry->begin_TPC(cryo.ID()),
                                    tend = geometry->end_TPC(cryo.ID());
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    while (iTPC != tend) {
      geo::TPCGeo const& TPC = *iTPC;
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
      iTPC++;
    }
     fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: fTPCVolumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    // fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
    fFiducialVolumes.emplace_back(XMin + fFiducialInset[0], XMax - fFiducialInset[1], 
                                  YMin + fFiducialInset[2], YMax - fFiducialInset[3], 
                                  ZMin + fFiducialInset[4], ZMax - fFiducialInset[5]);
  }
}

bool sbn::NuVertexChargeTree::InFV(TVector3 pos) {
  for (auto const &FV: fFiducialVolumes) {
    if (FV.ContainsPosition(pos)) return true;
  }
  return false;
}

const simb::MCParticle *Genie2G4MCParticle(
  const simb::MCParticle &genie_part,
  const simb::MCTruth &mctruth,
  const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles,
  const std::vector<const sim::GeneratedParticleInfo *> infos) {

  const simb::MCParticle *ret = NULL;
  for (int iparticle = 0; iparticle < (int)g4_mcparticles.size(); iparticle++) {
    if (infos[iparticle]->hasGeneratedParticleIndex() &&
        (int)infos[iparticle]->generatedParticleIndex() < mctruth.NParticles() && // TODO: why is this number sometimes bigger than the number of particles?
        mctruth.GetParticle(infos[iparticle]->generatedParticleIndex()).TrackId() == genie_part.TrackId() &&
        g4_mcparticles[iparticle]->Process() == "primary" /* TODO: will have to remove this restriction to include secondary particles*/) {

      // if a genie particle re-scatters in g4 and makes more particles, then multiple g4 particles can match to a 
      // genie particle. Thus, we also check that the start location of the associated genie particle matches the g4
      // and that the pdgid matches (to be on the safe side)
      //
      // Note that this should be accounted for by requiring the process to be primary. This is a bit of redundancy.
      const simb::MCParticle& matched_genie_particle = mctruth.GetParticle(infos[iparticle]->generatedParticleIndex());
      if ((matched_genie_particle.Position().Vect() - g4_mcparticles[iparticle]->Position().Vect()).Mag() < 1e-4 &&
        matched_genie_particle.PdgCode() == g4_mcparticles[iparticle]->PdgCode()) {

        // this should only be true for one particle
        assert(ret == NULL);
        ret = g4_mcparticles[iparticle].get();
      }
    }
  }
  return ret;
}

void sbn::NuVertexChargeTree::Clear() {
  fEvt = 0;

  fNuE = 0;
  fNuMode = -1;
  fNuIntType = -1;

  fFSPE.clear();
  fFSPKE.clear();
  fFSPPID.clear();
  fFSPMatched.clear();
  fFSPMtachedPrimary.clear();

  fRecoVertDist = 0.;
  fVertexHitCharge.clear();
  fVertexHitDist.clear();
  fVertexHitPlane.clear();
  fNVertexPFPs.clear();
  fVertexHitPFPIDs.clear();
  fVertexHitPFPPCvals.clear();

  fNSliceParticles = 0;
  fNSlicePrimaryTracks = 0;
  fNSlicePrimaryShowers = 0;
}

void sbn::NuVertexChargeTree::FillMeta(const art::Event &evt) {
  fEvt = evt.event();
  fIEvt += 1;
}

void sbn::NuVertexChargeTree::FillSlice(const std::vector<art::Ptr<recob::PFParticle>> &slice_particles) {
  fNSliceParticles = slice_particles.size();
  for (unsigned i = 0; i < slice_particles.size(); i++) {
    if (slice_particles[i]->IsPrimary() && lar_pandora::LArPandoraHelper::IsTrack(slice_particles[i])) {
      fNSlicePrimaryTracks ++;
    }
    else if (slice_particles[i]->IsPrimary() && lar_pandora::LArPandoraHelper::IsShower(slice_particles[i])) {
      fNSlicePrimaryShowers ++;
    }
  }
}

void sbn::NuVertexChargeTree::FillVertexHits(const std::vector<art::Ptr<sbn::VertexHit>> &hits) {
  for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
    const sbn::VertexHit &hit = *hits[i_hit];
    unsigned plane = hit.wire.Plane;
    fVertexHitPlane.push_back(plane);
    fVertexHitCharge.push_back(hit.charge);
    fVertexHitDist.push_back(hit.proj_dist_to_vertex);
    if (hit.nearbyPFPIDs.size() > MAX_PFP) {
      std::cout << "ERROR: MORE PRIMARY PFPs (" << hit.nearbyPFPIDs.size() << ") than tree slots (" << MAX_PFP << "(\n";
    }
    fNVertexPFPs.push_back(hit.nearbyPFPIDs.size());
    std::array<int, MAX_PFP> emptyintArray {};
    std::array<float, MAX_PFP> emptyfloatArray {};
    fVertexHitPFPIDs.push_back(emptyintArray);
    fVertexHitPFPPCvals.push_back(emptyfloatArray);
    for (unsigned i_part = 0; i_part < hit.nearbyPFPIDs.size() && i_part < MAX_PFP; i_part++) {
      fVertexHitPFPIDs.back()[i_part] = hit.nearbyPFPIDs[i_part];
      fVertexHitPFPPCvals.back()[i_part] = hit.PCproj[i_part];
    }

  }
}

void sbn::NuVertexChargeTree::FillNeutrino(const simb::MCTruth &nu, 
   const recob::Vertex &vert,
   const std::vector<std::vector<std::pair<int, float>>> &slice_pfp_matches,
   const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles, 
   const std::vector<const sim::GeneratedParticleInfo *> infos) {

  art::ServiceHandle<cheat::BackTrackerService> backtracker;

  TVector3 vpos(vert.position().X(), vert.position().Y(), vert.position().Z());

  fRecoVertDist = (nu.GetNeutrino().Nu().Position().Vect() - vpos).Mag();

  fNuE = nu.GetNeutrino().Nu().E();

  // save each stable final-state particle
  for (int i_part = 0; i_part < nu.NParticles(); i_part++) {
    const simb::MCParticle &particle = nu.GetParticle(i_part);
    // is stable?
    if (particle.StatusCode() != 1) continue;

    // Find the corresponding G4 particle
    const simb::MCParticle *G4 = Genie2G4MCParticle(particle, nu, g4_mcparticles, infos);

    if (!G4) continue;

    // lookup its hits

    // save some crap
    fFSPE.push_back(particle.Momentum().E());
    fFSPKE.push_back(particle.Momentum().E() - particle.Momentum().M());
    fFSPPID.push_back(particle.PdgCode());

    // see if there is a match
    //
    // Total up the energy fromthis particle
    std::vector<const sim::IDE*> particle_ides(backtracker->TrackIdToSimIDEs_Ps(particle.TrackId()));
    float thisVisE = 0.;
    for (const sim::IDE *ide: particle_ides) thisVisE += ide->energy / 1000.;
    bool found_match = false;
    bool primary_match = false;
    for (unsigned i_pfp = 0; i_pfp < slice_pfp_matches.size(); i_pfp++) {
      for (unsigned i_match = 0; i_match < slice_pfp_matches[i_pfp].size(); i_match++) {

        if (slice_pfp_matches[i_pfp][i_match].first == G4->TrackId() &&
            (slice_pfp_matches[i_pfp][i_match].second / thisVisE > 0.5)) {
          found_match = true;     
          primary_match = i_match == 0;
          break;
        }
      }
    }

    fFSPMatched.push_back(found_match);
    fFSPMtachedPrimary.push_back(primary_match);

  }
}

void sbn::NuVertexChargeTree::analyze(art::Event const& evt)
{
  // services
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> inventory_service;
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  // input data
  std::vector<std::vector<art::Ptr<recob::Slice>>> sliceList;
  std::vector<art::FindManyP<recob::Hit>> sliceHitList;
  std::vector<art::FindManyP<recob::PFParticle>> slicePFParticleList;

  std::vector<std::vector<art::Ptr<recob::PFParticle>>> particleList;
  std::vector<art::FindManyP<recob::Cluster>> particleClusterList;
  std::vector<art::FindManyP<recob::Vertex>> particleVertexList;
  std::vector<art::FindManyP<sbn::VertexHit>> particleVertexHitList;
  
  for (unsigned i = 0; i < fPFParticleTags.size(); i++) {
    art::Handle<std::vector<recob::Slice>> slice_handle;
    evt.getByLabel(fPFParticleTags.at(i), slice_handle);
    sliceList.emplace_back();
    art::fill_ptr_vector(sliceList.back(), slice_handle);
    sliceHitList.emplace_back(sliceList.back(), evt, fPFParticleTags.at(i));
    slicePFParticleList.emplace_back(sliceList.back(), evt, fPFParticleTags.at(i));

    art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
    evt.getByLabel(fPFParticleTags.at(i), pfparticle_handle);
    particleList.emplace_back();
    art::fill_ptr_vector(particleList.back(), pfparticle_handle);
    particleClusterList.emplace_back(particleList.back(), evt, fPFParticleTags.at(i));
    particleVertexList.emplace_back(particleList.back(), evt, fPFParticleTags.at(i));
    particleVertexHitList.emplace_back(particleList.back(), evt, fVertexHitTags.at(i));
  }

  // flatten
  std::vector<art::Ptr<recob::Slice>> slices;
  std::vector<std::vector<art::Ptr<recob::PFParticle>>> sliceParticles;
  std::vector<std::vector<art::Ptr<recob::Hit>>> sliceHits;

  std::vector<art::Ptr<recob::PFParticle>> particles;
  std::vector<std::vector<art::Ptr<recob::Cluster>>> particleClusters;
  std::vector<std::vector<art::Ptr<recob::Vertex>>> particleVertexs;
  std::vector<std::vector<art::Ptr<recob::Hit>>> particleHits;
  std::vector<std::vector<art::Ptr<sbn::VertexHit>>> particleVertexHits;

  for (unsigned i = 0; i < sliceList.size(); i++) {
    slices.insert(slices.end(), sliceList[i].begin(), sliceList[i].end());
    for (unsigned j = 0; j < sliceList[i].size(); j++) {
      sliceParticles.push_back(slicePFParticleList[i].at(j));
      sliceHits.push_back(sliceHitList[i].at(j));
    }
  }

  for (unsigned i = 0; i < particleList.size(); i++) {
    particles.insert(particles.end(), particleList[i].begin(), particleList[i].end());

    for (unsigned j = 0; j < particleList[i].size(); j++) {
      particleVertexs.push_back(particleVertexList[i].at(j));
      particleVertexHits.push_back(particleVertexHitList[i].at(j));

      particleClusters.push_back(particleClusterList[i].at(j));
      particleHits.emplace_back();
      art::FindManyP<recob::Hit> clusterHits(particleClusters.back(), evt, fPFParticleTags.at(i));
      for (unsigned i_clus = 0; i_clus < particleClusters.back().size(); i_clus++) {
        particleHits.back().insert(particleHits.back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
      }
    }
  }

  // also get truth stuff
  art::Handle<std::vector<simb::MCTruth>> trueNuHandle;
  evt.getByLabel("generator", trueNuHandle);
  std::vector<art::Ptr<simb::MCTruth>> trueNus;
  art::fill_ptr_vector(trueNus, trueNuHandle);

  art::Handle<std::vector<simb::MCParticle>> trueParticleHandle;
  evt.getByLabel("largeant", trueParticleHandle);
  std::vector<art::Ptr<simb::MCParticle>> trueParticles;
  art::fill_ptr_vector(trueParticles, trueParticleHandle);

  // get associations from MCTruth information to truth 
  art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truth_to_particles(trueNus, evt, "largeant");

  // iterate over each neutrino
  for (unsigned i_nu = 0; i_nu < trueNus.size(); i_nu++) {
    const simb::MCTruth &nu = *trueNus[i_nu];
    if (!nu.NeutrinoSet()) continue; // ignore non-neutrino generators

    // require inside fiducial volume
    if (!InFV(nu.GetNeutrino().Nu().Position().Vect())) continue;

    // require CC numu
    if (nu.GetNeutrino().CCNC() != 0 || abs(nu.GetNeutrino().Nu().PdgCode()) != 14) continue;

    // try to do reco-matching
    //
    // Total up the visible energy of this neutrino
    float visE = 0.;
    for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
      const simb::MCParticle &particle = *trueParticles[i_g4];
      art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(particle.TrackId());
      if (truth == trueNus[i_nu]) {
        std::vector<const sim::IDE*> particle_ides(backtracker->TrackIdToSimIDEs_Ps(particle.TrackId()));
        for (const sim::IDE *ide: particle_ides) visE += ide->energy / 1000.;
      }
    }

    int slice_match = -1;
    // see if the visible energy from any slice matches
    for (unsigned i_slice = 0; i_slice < slices.size(); i_slice++) {
      const recob::Slice &slice = *slices[i_slice];
      (void) slice;
      float matchingE = 0.;
      float totalE = 0.;
      std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, sliceHits.at(i_slice));
      for (unsigned i_match = 0; i_match < matches.size(); i_match++) {
        totalE += matches[i_match].second / 1000.;
      }

      for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
        for (unsigned i_match = 0; i_match < matches.size(); i_match++) {
          if (matches[i_match].first == trueParticles[i_g4]->TrackId()) {
            art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(trueParticles[i_g4]->TrackId());
            if (truth == trueNus[i_nu]) {
              matchingE += matches[i_match].second / 1000.;
            }
            break;
          }
        }
      }

      if (matchingE / visE > 0.5) {
        if (matchingE / totalE > 0.5) slice_match = i_slice;
        break;
      }
    }

    if (slice_match < 0) continue;
    // We found a matching slice! 

    // Do reco stuff
    //
    // Find the primary PFParticle of the slice -- that one is the key
    std::vector<unsigned> slice_pfp_inds;
    int primary = -1;
    const std::vector<art::Ptr<recob::PFParticle>> thisSliceParticles = sliceParticles[slice_match];
    for (unsigned i_pfp = 0; i_pfp < thisSliceParticles.size(); i_pfp++) {
      int ind = -1;
      for (unsigned j_pfp = 0; j_pfp < particles.size(); j_pfp++) {
        if (thisSliceParticles[i_pfp] == particles[j_pfp]) {
          ind = j_pfp;
          break;
        }
      }
      assert(ind >= 0);
      slice_pfp_inds.push_back((unsigned)ind);

      if (thisSliceParticles[i_pfp]->IsPrimary()) {
        primary = i_pfp;
      }
    }
    if (primary < 0) continue;

    // Collect reco stuff
    const std::vector<art::Ptr<sbn::VertexHit>> &thisVertexHits = particleVertexHits.at(slice_pfp_inds[primary]);
    const recob::Vertex &vert = *particleVertexs.at(slice_pfp_inds[primary]).at(0);

    std::vector<std::vector<std::pair<int, float>>> thisSlicePFPMatches; 
    for (unsigned i_pfp = 0; i_pfp < thisSliceParticles.size(); i_pfp++) {
      std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, particleHits.at(slice_pfp_inds[i_pfp]));
      // sort
      std::sort(matches.begin(), matches.end(), 
        [](auto const &lhs, auto const &rhs) {
          return lhs.second > rhs.second;
      });
      thisSlicePFPMatches.push_back(matches);
    }

    Clear();
    FillMeta(evt);

    // save the true-neutrino info
    FillNeutrino(nu, vert, thisSlicePFPMatches, truth_to_particles.at(i_nu), truth_to_particles.data(i_nu));

    // The vertex hits!
    FillVertexHits(thisVertexHits);

    FillSlice(thisSliceParticles);

  }
}

DEFINE_ART_MODULE(sbn::NuVertexChargeTree)
