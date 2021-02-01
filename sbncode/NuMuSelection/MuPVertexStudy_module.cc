////////////////////////////////////////////////////////////////////////
// Class:       MuPVertexStudy
// Plugin Type: analyzer (art v3_03_01)
// File:        MuPVertexStudy_module.cc
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
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "TH1D.h"
#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"

namespace numu {
  class MuPVertexStudy;
}

struct TrackHistos {
  TH1D *cos_open_angle;
  TH1D *muon_length;
  TH1D *proton_length;
  TH1D *muon_completion;
  TH1D *proton_completion;
  TH1D *muon_purity;
  TH1D *proton_purity;
};

struct TrueHistos {
  TH1D *cos_open_angle;
  TH1D *muon_length;
  TH1D *proton_length;
};

float Completion(const simb::MCParticle &particle, float particleE, const std::vector<std::pair<int, float>> &matches, const std::vector<simb::MCParticle> &particles);
float Purity(const simb::MCParticle &particle, float totalE, const std::vector<std::pair<int, float>> &matches, const std::vector<simb::MCParticle> &particles);
float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                      const std::vector<geoalgo::AABox> &boxes);

class numu::MuPVertexStudy : public art::EDAnalyzer {
public:
  explicit MuPVertexStudy(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuPVertexStudy(MuPVertexStudy const&) = delete;
  MuPVertexStudy(MuPVertexStudy&&) = delete;
  MuPVertexStudy& operator=(MuPVertexStudy const&) = delete;
  MuPVertexStudy& operator=(MuPVertexStudy&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  //TH1D *nPFParticle;
  TrueHistos fTrueHistos;
  std::array<TrackHistos, 3> fTrackHistos; 
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  void InitTrue();
  void InitReco(unsigned i, const std::string &prefix);

  void FillReco(unsigned i, const simb::MCParticle &muon, const simb::MCParticle &proton, float muon_completion, float proton_completion, float muon_purity, float proton_purity);

  void FillTrue(const simb::MCParticle &muon, const simb::MCParticle &proton);
  float ParticleLength(const simb::MCParticle &particle);

};

void numu::MuPVertexStudy::InitTrue() {
  art::ServiceHandle<art::TFileService> tfs; 

  fTrueHistos.cos_open_angle = tfs->make<TH1D>("CosOpenAngle", "", 64, -3.2, 3.2);
  fTrueHistos.muon_length = tfs->make<TH1D>("mu_length", "", 100, 0., 500.);
  fTrueHistos.proton_length = tfs->make<TH1D>("p_length", "", 100, 0., 200.);
}

void numu::MuPVertexStudy::InitReco(unsigned i, const std::string &prefix) {
  art::ServiceHandle<art::TFileService> tfs; 

  fTrackHistos[i].cos_open_angle = tfs->make<TH1D>((prefix + "reco_CosOpenAngle").c_str(), "", 64, -3.2, 3.2);
  fTrackHistos[i].muon_length = tfs->make<TH1D>((prefix + "reco_mu_length").c_str(), "", 100, 0., 500.);
  fTrackHistos[i].proton_length = tfs->make<TH1D>((prefix + "reco_p_length").c_str(), "", 100, 0., 200.);

  fTrackHistos[i].muon_completion = tfs->make<TH1D>((prefix + "reco_mu_completion").c_str(), "", 100, 0., 1.);
  fTrackHistos[i].proton_completion = tfs->make<TH1D>((prefix + "reco_p_completion").c_str(), "", 100, 0., 1.);
  fTrackHistos[i].muon_purity = tfs->make<TH1D>((prefix + "reco_mu_purity").c_str(), "", 100, 0., 1.);
  fTrackHistos[i].proton_purity = tfs->make<TH1D>((prefix + "reco_p_purity").c_str(), "", 100, 0., 1.);

}


numu::MuPVertexStudy::MuPVertexStudy(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  InitTrue();

  std::vector<std::string> prefixes {"both", "mu", "p"};
  for (unsigned i = 0; i < 3; i++) {
    InitReco(i, prefixes[i]);
  }

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
  }

}

float numu::MuPVertexStudy::ParticleLength(const simb::MCParticle &particle) {
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

void numu::MuPVertexStudy::FillTrue(const simb::MCParticle &muon, const simb::MCParticle &proton) {
  fTrueHistos.muon_length->Fill(ParticleLength(muon));
  fTrueHistos.proton_length->Fill(ParticleLength(proton));

  fTrueHistos.cos_open_angle->Fill(muon.Momentum().Vect().Dot(proton.Momentum().Vect()) / (muon.Momentum().Vect().Mag() * proton.Momentum().Vect().Mag()));
}

void numu::MuPVertexStudy::FillReco(unsigned i, const simb::MCParticle &muon, const simb::MCParticle &proton, float muon_completion, float proton_completion, float muon_purity, float proton_purity) {

  fTrackHistos[i].muon_completion->Fill(muon_completion);
  fTrackHistos[i].proton_completion->Fill(proton_completion);
  fTrackHistos[i].muon_purity->Fill(muon_purity);
  fTrackHistos[i].proton_purity->Fill(proton_purity);

  fTrackHistos[i].muon_length->Fill(ParticleLength(muon));
  fTrackHistos[i].proton_length->Fill(ParticleLength(proton));

  fTrackHistos[i].cos_open_angle->Fill(muon.Momentum().Vect().Dot(proton.Momentum().Vect()) / (muon.Momentum().Vect().Mag() * proton.Momentum().Vect().Mag()));
}

void numu::MuPVertexStudy::analyze(art::Event const& ev)
{
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(ev);
  art::ServiceHandle<cheat::BackTrackerService> backtracker;

  const std::vector<simb::MCParticle> &mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  const simb::MCParticle *muon = NULL;
  const simb::MCParticle *proton = NULL;

  // get the muon and proton
  for (unsigned i = 0; i < mcparticle_list.size(); i++) {
    if (mcparticle_list[i].PdgCode() == 13 && mcparticle_list[i].Process() == "primary") {
      assert(muon == NULL);
      muon = &mcparticle_list[i];
    }
    if (mcparticle_list[i].PdgCode() == 2212 && mcparticle_list[i].Process() == "primary") {
      assert(proton == NULL);
      proton = &mcparticle_list[i];
    }
  }

  assert(muon != NULL);
  assert(proton != NULL);

  FillTrue(*muon, *proton);

  float muonE = 0.;
  std::vector<const sim::IDE*> muon_ides(backtracker->TrackIdToSimIDEs_Ps(muon->TrackId()));
  for (const sim::IDE *ide: muon_ides) muonE += ide->energy;

  float protonE = 0.;
  std::vector<const sim::IDE*> proton_ides(backtracker->TrackIdToSimIDEs_Ps(proton->TrackId()));
  for (const sim::IDE *ide: proton_ides) protonE += ide->energy;

  const auto &particle_handle = ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");
  const std::vector<recob::PFParticle> &particle_list = *particle_handle;
  art::FindManyP<recob::Cluster> particles_to_clusters(particle_handle, ev, "pandora");

  for (unsigned i_part = 0; i_part < particle_list.size(); i_part++) {
    const std::vector<art::Ptr<recob::Cluster>> clusters = particles_to_clusters.at(i_part); 
    art::FindManyP<recob::Hit> clusters_to_hits(clusters, ev, "pandora");

    std::vector<art::Ptr<recob::Hit>> hits;
    for (unsigned i_clus = 0; i_clus < clusters.size(); i_clus++) {
      const std::vector<art::Ptr<recob::Hit>> &this_hits = clusters_to_hits.at(i_clus);
      hits.insert(hits.end(), this_hits.begin(), this_hits.end());
    }

    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, hits);
    float totalE = CAFRecoUtils::TotalHitEnergy(clock_data, hits);

    float muon_completion = Completion(*muon, muonE, matches, mcparticle_list);
    float proton_completion = Completion(*proton, protonE, matches, mcparticle_list);

    float muon_purity = Purity(*muon, totalE, matches, mcparticle_list);
    float proton_purity = Purity(*proton, totalE, matches, mcparticle_list);

    if (muon_completion > 0.5 && proton_completion > 0.5) {
      FillReco(0, *muon, *proton, muon_completion, proton_completion, muon_purity, proton_purity);
    }
    else if (muon_completion > 0.5) {
      FillReco(1, *muon, *proton, muon_completion, proton_completion, muon_purity, proton_purity);
    }
    else if (proton_completion > 0.5) {
      FillReco(2, *muon, *proton, muon_completion, proton_completion, muon_purity, proton_purity);
    }
  } 
}

float Completion(const simb::MCParticle &particle, float particleE, const std::vector<std::pair<int, float>> &matches, const std::vector<simb::MCParticle> &particles) {
  for (auto const &pair: matches) {
    if (pair.first == particle.TrackId()) {
      return pair.second / particleE;
    } 
  }
  return 0.;
}

float Purity(const simb::MCParticle &particle, float totalE, const std::vector<std::pair<int, float>> &matches, const std::vector<simb::MCParticle> &particles) {
  float matchE = 0.;
  for (auto const &pair: matches) {
    bool matches = false;
    if (pair.first == particle.TrackId()) matches = true;
    else {
      for (int i_d = 0; i_d < particle.NumberDaughters(); i_d++) {
        if (pair.first == particle.Daughter(i_d)) {
          bool showerDaughter = true;
          for (unsigned i_p = 0; i_p < particles.size(); i_p++) {
            if (pair.first == particles[i_p].TrackId()) {
              showerDaughter = particles[i_p].Process() == "muIoni";
              break;
            }
          }
          matches = showerDaughter;
          break;
        }
      }
    }
    if (matches) matchE += pair.second;
  }
  return matchE / totalE;
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

DEFINE_ART_MODULE(numu::MuPVertexStudy)
