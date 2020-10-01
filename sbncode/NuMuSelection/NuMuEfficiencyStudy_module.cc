////////////////////////////////////////////////////////////////////////
// Class:       NuMuEfficiencyStudy
// Plugin Type: analyzer (art v3_03_01)
// File:        NuMuEfficiencyStudy_module.cc
//
// Generated at Wed Mar 25 12:10:37 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_08_00.
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

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/AnalysisBase/T0.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbncode/LArRecoProducer/Products/RangeP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "TH1D.h"
#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"

namespace numu {
  class NuMuEfficiencyStudy;
}

struct Histos {
  TH1D *neutrino_energy;
  TH1D *muon_momentum;
  TH1D *muon_length;
  TH1D *vertex_x;
  TH1D *vertex_y;
  TH1D *vertex_z;
};

float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                      const std::vector<geoalgo::AABox> &boxes);

const simb::MCParticle *Genie2G4MCParticle(
  const simb::MCParticle &genie_part,
  const simb::MCTruth &mctruth,
  const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles,
  const std::vector<const sim::GeneratedParticleInfo *> infos);

class numu::NuMuEfficiencyStudy : public art::EDAnalyzer {
public:
  explicit NuMuEfficiencyStudy(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMuEfficiencyStudy(NuMuEfficiencyStudy const&) = delete;
  NuMuEfficiencyStudy(NuMuEfficiencyStudy&&) = delete;
  NuMuEfficiencyStudy& operator=(NuMuEfficiencyStudy const&) = delete;
  NuMuEfficiencyStudy& operator=(NuMuEfficiencyStudy&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  Histos fTruth;
  std::vector<Histos> fReco;

  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  static const unsigned n_reco_eff = 9;

  std::array<float, 6> fFiducialInset;
  std::vector<geo::BoxBoundedGeo> fFiducialVolumes;

  bool InFV(TVector3 pos);

  void InitHistos();

  void InitHisto(Histos &h, const std::string &prefix);

  void Fill(Histos &h, const simb::MCTruth &interaction, const simb::MCParticle &G4lepton);

  float ParticleLength(const simb::MCParticle &particle);

};

bool numu::NuMuEfficiencyStudy::InFV(TVector3 pos) {
  for (auto const &FV: fFiducialVolumes) {
    if (FV.ContainsPosition(pos)) return true;
  }
  return false;
}

void numu::NuMuEfficiencyStudy::InitHistos() {
  InitHisto(fTruth, "truth");

  std::array<std::string, n_reco_eff> effnames {"has_muon", "has_reco_muon", "muon_is_fid", "has_slice", "has_pure_slice", "has_reco_slice", "has_fid_slice", "slice_has_muon", "has_fid_pure_slice"};
  for (unsigned i = 0; i < n_reco_eff; i++) {
    for (unsigned j = 0; j < n_reco_eff; j++) {
      std::string name;
      if (i == j) name = effnames[i];
      else name = effnames[i] + "_" + effnames[j];
      fReco.emplace_back();
      InitHisto(fReco.back(), name);
    }
  }
}

void numu::NuMuEfficiencyStudy::InitHisto(Histos &h, const std::string &prefix) {
  art::ServiceHandle<art::TFileService> tfs; 

  h.neutrino_energy = tfs->make<TH1D>((prefix + "_nuE").c_str(), "nuE", 100, 0, 4);
  h.muon_momentum = tfs->make<TH1D>((prefix + "_muP").c_str(), "muP", 100, 0, 4);
  h.muon_length = tfs->make<TH1D>((prefix + "_muL").c_str(), "muL", 100, 0, 600);

  float XMin = std::min_element(fFiducialVolumes.begin(), fFiducialVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
  float XMax = std::max_element(fFiducialVolumes.begin(), fFiducialVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
  float YMin = std::min_element(fFiducialVolumes.begin(), fFiducialVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
  float YMax = std::max_element(fFiducialVolumes.begin(), fFiducialVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
  float ZMin = std::min_element(fFiducialVolumes.begin(), fFiducialVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();
  float ZMax = std::max_element(fFiducialVolumes.begin(), fFiducialVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();
                 
  h.vertex_x = tfs->make<TH1D>((prefix + "_v_x").c_str(), "v_x", 200, XMin, XMax); 
  h.vertex_y = tfs->make<TH1D>((prefix + "_v_y").c_str(), "v_y", 200, YMin, YMax); 
  h.vertex_z = tfs->make<TH1D>((prefix + "_v_z").c_str(), "v_z", 200, ZMin, ZMax); 
}

numu::NuMuEfficiencyStudy::NuMuEfficiencyStudy(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fFiducialInset({p.get<float>("xmin"), p.get<float>("xmax"), p.get<float>("ymin"), p.get<float>("ymax"), p.get<float>("zmin"), p.get<float>("zmax")}) 
  // More initializers here.
{

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

    fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
    fFiducialVolumes.emplace_back(XMin + fFiducialInset[0], XMax - fFiducialInset[1], 
                                  YMin + fFiducialInset[2], YMax - fFiducialInset[3], 
                                  ZMin + fFiducialInset[4], ZMax - fFiducialInset[5]);
  }

  InitHistos();
}

float numu::NuMuEfficiencyStudy::ParticleLength(const simb::MCParticle &particle) {
  std::vector<geoalgo::AABox> aa_volumes;
  for (const geo::BoxBoundedGeo &AV: fActiveVolumes) {
    if (AV.ContainsPosition(particle.Position().Vect())) {
      aa_volumes.emplace_back(AV.MinX(), AV.MinY(), AV.MinZ(), AV.MaxX(), AV.MaxY(), AV.MaxZ());
      break;
    }
  }

  float length = 0.;
  for (unsigned i = 1; i < particle.NumberTrajectoryPoints(); i++) {
    length += ContainedLength(particle.Position(i-1).Vect(), particle.Position(i).Vect(), aa_volumes);
  }

  return length;
}

void numu::NuMuEfficiencyStudy::Fill(Histos &h, const simb::MCTruth &truth, const simb::MCParticle &G4lepton) {
  h.neutrino_energy->Fill(truth.GetNeutrino().Nu().E());
  h.muon_momentum->Fill(truth.GetNeutrino().Lepton().Momentum().Vect().Mag());
  h.muon_length->Fill(ParticleLength(G4lepton));

  h.vertex_x->Fill(truth.GetNeutrino().Nu().Position().X());
  h.vertex_y->Fill(truth.GetNeutrino().Nu().Position().Y());
  h.vertex_z->Fill(truth.GetNeutrino().Nu().Position().Z());
}

void numu::NuMuEfficiencyStudy::analyze(art::Event const& ev)
{
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> inventory_service;

  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(ev);

  art::ValidHandle<std::vector<simb::MCTruth>> mctruth_handle = ev.getValidHandle<std::vector<simb::MCTruth>>("generator");
  const std::vector<simb::MCTruth> &mctruth_list = *mctruth_handle;
  std::vector<art::Ptr<simb::MCTruth>> mctruth_ptrs;
  art::fill_ptr_vector(mctruth_ptrs, mctruth_handle);

  art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truth_to_particles(mctruth_handle, ev, "largeant");
  const std::vector<simb::MCParticle> &g4_particles = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  std::vector<int> truth_to_numucc;
  std::vector<simb::MCTruth> numucc;
  std::vector<simb::MCParticle> G4muons;
  // iterate over the true interactions
  for (unsigned i = 0; i < mctruth_list.size(); i++) {
    const simb::MCTruth &truth = mctruth_list[i];
    TVector3 vertex(truth.GetNeutrino().Nu().Position().X(), truth.GetNeutrino().Nu().Position().Y(), truth.GetNeutrino().Nu().Position().Z());
    // fiducial numu CC interactions
    if (InFV(vertex) && truth.GetNeutrino().CCNC() == 0 && abs(truth.GetNeutrino().Nu().PdgCode()) == 14) {
      // find the muon in G4
      const simb::MCParticle *G4muon = Genie2G4MCParticle(truth.GetNeutrino().Lepton(), truth, truth_to_particles.at(i), truth_to_particles.data(i));
      Fill(fTruth, truth, *G4muon);
      numucc.push_back(truth);
      G4muons.push_back(*G4muon);
      truth_to_numucc.push_back(numucc.size()-1);
    }
    else truth_to_numucc.push_back(-1);
  }

  // no true numu CC? Return
  if (numucc.size() == 0) return;

  // get the deposited energy for each muon
  std::vector<float> muonVisE(G4muons.size(), 0.);
  for (unsigned i = 0; i < G4muons.size(); i++) {
    std::vector<const sim::IDE*> particle_ides(backtracker->TrackIdToSimIDEs_Ps(G4muons[i].TrackId()));
    for (const sim::IDE *ide: particle_ides) muonVisE[i] += ide->energy / 1000.;
  }

  // get the deposited energy from primary particles in the neutrino
  std::vector<float> numuVisE(numucc.size(), 0.);
  for (unsigned i_g4 = 0; i_g4 < g4_particles.size(); i_g4++) {
    const simb::MCParticle &particle = g4_particles[i_g4];
    if (particle.Process() != "primary") continue;
    art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(particle.TrackId());
    for (unsigned i_truth = 0; i_truth < mctruth_ptrs.size(); i_truth++) {
      if (truth == mctruth_ptrs[i_truth]) {
        if (truth_to_numucc[i_truth] != -1) {
          std::vector<const sim::IDE*> particle_ides(backtracker->TrackIdToSimIDEs_Ps(particle.TrackId()));
          for (const sim::IDE *ide: particle_ides) numuVisE[truth_to_numucc[i_truth]] += ide->energy / 1000.;
        }
        break;
      }
    }
  } 

  // now collect the reco information 
  //
  // does a muon track exist?
  std::vector<int> matching_track(G4muons.size(), -1);

  art::ValidHandle<std::vector<recob::PFParticle>> particle_handle = ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");
  art::FindManyP<recob::Track> tracks(particle_handle, ev, "pandoraTrack");
  art::FindManyP<larpandoraobj::PFParticleMetadata> metadatas(particle_handle, ev, "pandora");
  art::FindManyP<recob::Vertex> particles_to_vertex(particle_handle, ev, "pandora");
  for (unsigned i_part = 0; i_part < particle_handle->size(); i_part++) {
    if (!tracks.at(i_part).size()) continue;
    art::FindManyP<recob::Hit> fmHits(tracks.at(i_part), ev, "pandoraTrack");
    const std::vector<art::Ptr<recob::Hit>> &this_hits = fmHits.at(0);

    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, this_hits);

    // sort by matching energy
    std::sort(matches.begin(), matches.end(), 
      [] (auto const &lhs, auto const &rhs) { return lhs.second > rhs.second; });

    for (unsigned i_cc = 0; i_cc < G4muons.size(); i_cc++) {
      if (matches.size() && matches[0].first == G4muons[i_cc].TrackId() && matches[0].second / muonVisE[i_cc] > 0.5) {
        matching_track[i_cc] = particle_handle->at(i_part).Self();
      }
    }
  }

  // also find a matching slice
  std::vector<int> matching_slice(numucc.size(), -1);
  std::vector<float> matching_purity(numucc.size(), -1);

  art::ValidHandle<std::vector<recob::Slice>> slice_handle = ev.getValidHandle<std::vector<recob::Slice>>("pandora"); 
  art::FindManyP<recob::Hit> slice_hits(slice_handle, ev, "pandora");
  art::FindManyP<recob::PFParticle> slice_particles(slice_handle, ev, "pandora");

  for (unsigned i_slice = 0; i_slice < slice_handle->size(); i_slice++) {
    std::vector<float> matchingE(numuVisE.size(), 0.);
    float totalE = 0.;

    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, slice_hits.at(i_slice));
    for (unsigned i_match = 0; i_match < matches.size(); i_match++) {
      totalE += matches[i_match].second / 1000.;
    }

    for (unsigned i_g4 = 0; i_g4 < g4_particles.size(); i_g4++) {
      for (unsigned i_match = 0; i_match < matches.size(); i_match++) {
        if (matches[i_match].first == g4_particles[i_g4].TrackId()) {
          if (g4_particles[i_g4].Process() == "primary") {
            art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(g4_particles[i_g4].TrackId());
            for (unsigned i_truth = 0; i_truth < mctruth_ptrs.size(); i_truth++) {
              if (truth == mctruth_ptrs[i_truth]) {
                if (truth_to_numucc[i_truth] != -1) {
                  matchingE[truth_to_numucc[i_truth]] += matches[i_match].second / 1000.;
                }
                break;
              }
            }
          } 
          break;
        }
      }
    }

    unsigned best_match = std::distance(matchingE.begin(), std::max_element(matchingE.begin(), matchingE.end()));
    if (matchingE[best_match] / numuVisE[best_match] > 0.5) {
      matching_slice[best_match] = i_slice;
      matching_purity[best_match] = matchingE[best_match] / totalE; 
    }
  }

  // gather information for efficiency calculation
  for (unsigned i = 0; i < numucc.size(); i++) {
    bool has_muon_track = matching_track[i] != -1;
    bool has_slice = matching_slice[i] != -1;

    bool has_pure_slice = has_slice && matching_purity[i] > 0.5;

    bool muon_is_reco = false;
    bool muon_is_fid = false;
    if (has_muon_track) {
      const larpandoraobj::PFParticleMetadata &meta = *metadatas.at(matching_track[i]).at(0);
      auto const &properties = meta.GetPropertiesMap();
      std::cout << "Muon Track Properties:\n";
      for (auto pair: properties) std::cout << pair.first << std::endl;
      if (properties.count("IsClearCosmic")) {
        muon_is_reco = false;
      }
      else {
        muon_is_reco = true;
      }

      auto muon_vert = tracks.at(matching_track[i]).at(0)->Start();
      TVector3 muon_vert_v(muon_vert.X(), muon_vert.Y(), muon_vert.Z());
            
      muon_is_fid = InFV(muon_vert_v);
    }

    bool slice_is_reco = false;
    bool slice_is_fiducial = false;
    bool slice_has_muon = false;
    if (has_slice) {
      const std::vector<art::Ptr<recob::PFParticle>> &this_slice_particles = slice_particles.at(matching_slice[i]); 
      for (art::Ptr<recob::PFParticle> pfp: this_slice_particles) {
        if (pfp->IsPrimary()) {
          const larpandoraobj::PFParticleMetadata &meta = *metadatas.at(pfp.key()).at(0);
          auto const &properties = meta.GetPropertiesMap();
          std::cout << "Particle Properties:\n";
          for (auto pair: properties) std::cout << pair.first << std::endl;
          if (properties.count("IsNeutrino") && pfp->PdgCode() == 14) slice_is_reco = true;
          if (slice_is_reco) {
            const recob::Vertex &vert = *particles_to_vertex.at(pfp.key()).at(0);
            TVector3 vect(vert.position().X(), vert.position().Y(), vert.position().Z());
            slice_is_fiducial = InFV(vect);
          }
        }
        if (has_muon_track && matching_track[i] == (int)pfp->Self()) {
          slice_has_muon = true;
        }
      }
    }

    bool has_fid_pure_slice  = slice_is_fiducial  && has_pure_slice;

    std::array<bool, n_reco_eff> pass {has_muon_track, muon_is_reco, muon_is_fid, has_slice, has_pure_slice, slice_is_reco, slice_is_fiducial, slice_has_muon, has_fid_pure_slice};
    for (unsigned j = 0; j < n_reco_eff; j++) {
      for (unsigned k = 0; k < n_reco_eff; k++) {
        if (pass[j] && pass[k]) {
          Fill(fReco[j * n_reco_eff + k], numucc[i], G4muons[i]);
        }
      }
    }
  }
  
}

float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes) {
  static const geoalgo::GeoAlgo algo;

  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  // construct individual points
  geoalgo::Point_t p0(v0);
  geoalgo::Point_t p1(v1);

  // construct line segment
  geoalgo::LineSegment line(p0, p1);

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.Contain(p0) + box.Contain(p1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist)
    if (n_contained == 1) {
      auto intersections = algo.Intersection(line, box);
      // Because of floating point errors, it can sometimes happen
      // that there is 1 contained point but no "Intersections"
      // if one of the points is right on the edge
      if (intersections.size() == 0) {
        // determine which point is on the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        assert(p0_edge || p1_edge);
        // contained one is on edge -- can treat both as not contained
        //
        // In this case, no length
        if ((p0_edge && box.Contain(p0)) || (box.Contain(p1) && p1_edge))
          continue;
        // un-contaned one is on edge -- treat both as contained
        else if ((p0_edge && box.Contain(p1)) || (box.Contain(p0) && p1_edge)) {
	  length = (v1 - v0).Mag();
	  break;
        }
        else {
          assert(false); // bad
        }
      }
      // floating point errors can also falsely cause 2 intersection points
      //
      // in this case, one of the intersections must be very close to the 
      // "contained" point, so the total contained length will be about
      // the same as the distance between the two intersection points
      else if (intersections.size() == 2) {
        length += (intersections.at(0).ToTLorentzVector().Vect() - intersections.at(1).ToTLorentzVector().Vect()).Mag();
        continue;
      }
      // "Correct"/ideal case -- 1 intersection point
      else if (intersections.size() == 1) {
        // get TVector at intersection point
        TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
        length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() ); 
      }
      else assert(false); // bad
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      if (!(intersections.size() == 0 || intersections.size() == 2)) {
        // more floating point error fixes...
        //
        // figure out which points are near the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        // and which points are near the intersection
        TVector3 vint = intersections.at(0).ToTLorentzVector().Vect();

        bool p0_int = (v0 - vint).Mag() < tol;
        bool p1_int = (v1 - vint).Mag() < tol;
        // exactly one of them should produce the intersection
        assert((p0_int && p0_edge) != (p1_int && p1_edge));
        // void variables when assert-ions are turned off
        (void) p0_int; (void) p1_int;

        // both close to edge -- full length is contained
        if (p0_edge && p1_edge) {
          length += (v0 - v1).Mag();
        }
        // otherwise -- one of them is not on an edge, no length is contained
        else {}
      }
      // assert(intersections.size() == 0 || intersections.size() == 2);
      else if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;
}//ContainedLength

DEFINE_ART_MODULE(numu::NuMuEfficiencyStudy)

const simb::MCParticle *Genie2G4MCParticle(
  const simb::MCParticle &genie_part,
  const simb::MCTruth &mctruth,
  const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles,
  const std::vector<const sim::GeneratedParticleInfo *> infos) {

  // only stable final state particles are propogated by G4
  if (genie_part.StatusCode() != 1) return NULL;

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

