////////////////////////////////////////////////////////////////////////
// Class:       VertexChargeVacuumFinder
// Plugin Type: producer (art v3_02_06)
// File:        VertexChargeVacuum_module.cc
// Author:      grayputnam@uchicago.edu
//
// Art module for finding large charge depositions near a reconstructed
// vertex. Operates on input recob::Hit's. Produces sbn::VertexHit objects
// which conatin summary information useful for downstream reco/analysis.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/Calorimetry/INormalizeCharge.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Utilities/make_tool.h"

#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

namespace sbn {
  class VertexChargeVacuum;
}


class sbn::VertexChargeVacuum : public art::EDProducer {
public:
  explicit VertexChargeVacuum(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VertexChargeVacuum(VertexChargeVacuum const&) = delete;
  VertexChargeVacuum(VertexChargeVacuum&&) = delete;
  VertexChargeVacuum& operator=(VertexChargeVacuum const&) = delete;
  VertexChargeVacuum& operator=(VertexChargeVacuum&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;


private:
  // config
  art::InputTag fPFParticleLabel;
  art::InputTag fTrackLabel;
  float fHitVacuumRadius;
  bool fUseTrackSPRecovery;
  bool fCorrectSCE;
  bool fPositionsAreSCECorrected;
  bool fSelectNeutrino;
  std::vector<fhicl::ParameterSet> fNormToolConfig;
  std::vector<std::unique_ptr<INormalizeCharge>> fNormTools;

  // private data 
  calo::CalorimetryAlg fCaloAlg;

  // helpers
  geo::Point_t PositionAtWires(const geo::Point_t &p, const geo::GeometryCore *geo, const spacecharge::SpaceCharge *sce);
  geo::Point_t PositionAbsolute(const geo::Point_t &p, const geo::GeometryCore *geo, const spacecharge::SpaceCharge *sce);
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0);
};


sbn::VertexChargeVacuum::VertexChargeVacuum(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabel(p.get<std::string>("PFParticleLabel", "pandora")),
    fTrackLabel(p.get<std::string>("TrackLabel", "pandoraTrack")),
    fHitVacuumRadius(p.get<float>("HitVacuumRadius")),
    fUseTrackSPRecovery(p.get<bool>("UseTrackSPRecovery")),
    fCorrectSCE(p.get<bool>("CorrectSCE")),
    fPositionsAreSCECorrected(p.get<bool>("PositionsAreSCECorrected")),
    fSelectNeutrino(p.get<bool>("SelectNeutrino")),
    fNormToolConfig(p.get<std::vector<fhicl::ParameterSet>>("NormTools", {})),
    fCaloAlg(p.get<fhicl::ParameterSet >("CaloAlg"))
{

  for (const fhicl::ParameterSet &p: fNormToolConfig) {
    fNormTools.push_back(art::make_tool<INormalizeCharge>(p));
  }

  produces<std::vector<sbn::VertexHit>>();
  produces<art::Assns<recob::Slice, sbn::VertexHit>>();
  produces<art::Assns<recob::Hit, sbn::VertexHit>>();
  produces<art::Assns<recob::Vertex, sbn::VertexHit>>();
}

std::array<float, 2> HitVector(const recob::Hit &hit, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  float wire_distance = hit.WireID().Wire;
  // convert to cm
  float wire_distance_cm = wire_distance * geo->WirePitch();
  // and the time difference
  float time_distance = hit.PeakTime();
  // convert to cm
  float time_distance_cm = dprop.ConvertTicksToX(time_distance, hit.WireID());

  return {wire_distance_cm, time_distance_cm};
}

std::array<float, 2> VertexVector(const recob::Vertex &vert, const geo::PlaneID &plane, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  return {(float)(geo->WireCoordinate(vert.position(), plane) * geo->WirePitch()), (float)vert.position().X()};
}

float Vert2HitDistance(const recob::Hit &hit, const recob::Vertex &vert, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  std::array<float, 2> vert_v = VertexVector(vert, hit.WireID(), geo, dprop);
  std::array<float, 2> hit_v = HitVector(hit, geo, dprop);

  return sqrt((vert_v[0] - hit_v[0]) * (vert_v[0] - hit_v[0]) + (vert_v[1] - hit_v[1]) * (vert_v[1] - hit_v[1]));
}

// local helper function: Get the mean direction of a set of hits away from a vertex (projected on a plane)
std::array<float, 2> HitDirection(const std::vector<art::Ptr<recob::Hit>> &hits, const recob::Vertex &vert,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  if (!hits.size()) return {0., 0.};

  geo::Vector_t avg(0., 0., 0.);

  std::array<float, 2> vert_v = VertexVector(vert, hits[0]->WireID(), geo, dprop);
  geo::Point_t vert_p(vert_v[0], vert_v[1], 0.);

  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> hit_v = HitVector(*h, geo, dprop);
    geo::Point_t hit_p(hit_v[0], hit_v[1], 0.);
    avg += (hit_p - vert_p).Unit(); 
  }

  avg = avg.Unit();

  return {(float)avg.X(), (float)avg.Y()};
}

float TrackDirectionParallel(const recob::Track &trk, const geo::PlaneID &plane,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  double angleToVert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();
  double cosgamma = std::abs(std::sin(angleToVert)*trk.StartDirection().y() + std::cos(angleToVert)*trk.StartDirection().z());

  float ret = sqrt(cosgamma * cosgamma + trk.StartDirection().x() * trk.StartDirection().x());

  return ret;
  
}

geo::Point_t PlaceHitAlongTrack(const recob::Track &trk, const recob::Vertex &vert, const recob::Hit &hit,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  // project the vertex onto the hit Plane
  std::array<float, 2> v_plane = VertexVector(vert, hit.WireID(), geo, dprop); 

  // also get the 2d hit location
  std::array<float, 2> h_plane = HitVector(hit, geo, dprop);

  float t_dir_parallel_plane = TrackDirectionParallel(trk, hit.WireID(), geo, dprop);

  // Get the distance from the vertex to the hit on the plane.
  //
  // The track "should" be along the same direction as the vertex to the hit in the plane,
  // so this distance is the same as the distance along the track
  float plane_dist = (h_plane[0] - v_plane[0]) + (h_plane[1] - v_plane[1]);
 
  // Use this distance to place the hit along the track trajectory in 3D
  float dist_3d = plane_dist / t_dir_parallel_plane;

  // std::cout << "Placing Hit\n";
  // std::cout << "Hit wire: " << hit.WireID() << std::endl;
  // std::cout << "Hit plane coord: " << h_plane[0] << " " << h_plane[1] << std::endl;
  // std::cout << "Vtx plane coord: " << v_plane[0] << " " << v_plane[1] << std::endl;
  // std::cout << "Trk dir: " << trk.StartDirection().x() << " " << trk.StartDirection().y() << " " << trk.StartDirection().z() << std::endl;
  // std::cout << "Vtx pos: " << vert.position().x() << " " << vert.position().y() << " " << vert.position().z() << std::endl;

  // std::cout << "Trk parallel dir: " << t_dir_parallel_plane << std::endl;
  // std::cout << "Hit plane dist: " << plane_dist << std::endl;
  // std::cout << "Hit 3d dist: " << dist_3d << std::endl;

  geo::Vector_t trk_dir(trk.StartDirection().x(), trk.StartDirection().y(), trk.StartDirection().z());

  return vert.position() + trk_dir * dist_3d;
}

geo::Point_t sbn::VertexChargeVacuum::PositionAtWires(const geo::Point_t &p, const geo::GeometryCore *geo, const spacecharge::SpaceCharge *sce) {
  geo::TPCID tpc = geo->FindTPCAtPosition(p);

  if (tpc && fPositionsAreSCECorrected) return sbn::GetLocationAtWires(sce, geo, p, tpc);
  return p; 
}

geo::Point_t sbn::VertexChargeVacuum::PositionAbsolute(const geo::Point_t &p, const geo::GeometryCore *geo, const spacecharge::SpaceCharge *sce) {
  geo::TPCID tpc = geo->FindTPCAtPosition(p);

  if (tpc && !fPositionsAreSCECorrected) return sbn::GetLocation(sce, p, tpc);
  return p;
}

double sbn::VertexChargeVacuum::Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
    
  double ret = dQdx;

  // Normtools configured -- use those
  if (fNormTools.size()) {
    for (auto const &nt: fNormTools) {
      ret = nt->Normalize(ret, e, h, location, direction, t0);
    }
  }
  // Otherwise, fix using configured electron lifetime
  else {
    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const dprop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

    ret = ret * fCaloAlg.LifetimeCorrection(clock_data, dprop, h.PeakTime(), 0.);
  }

  return ret;
}

void sbn::VertexChargeVacuum::produce(art::Event& evt)
{
  // output stuff
  std::unique_ptr<std::vector<sbn::VertexHit>> outVHit(new std::vector<sbn::VertexHit>);
  std::unique_ptr<art::Assns<recob::Slice, sbn::VertexHit>> assn(new art::Assns<recob::Slice, sbn::VertexHit>);
  std::unique_ptr<art::Assns<recob::Vertex, sbn::VertexHit>> vtxAssn(new art::Assns<recob::Vertex, sbn::VertexHit>);
  std::unique_ptr<art::Assns<recob::Hit, sbn::VertexHit>> hitAssn(new art::Assns<recob::Hit, sbn::VertexHit>);

  art::PtrMaker<sbn::VertexHit> vhitPtrMaker {evt};

  // collect services
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  // If we are not correcting for SCE, then blank out the service
  if (!fCorrectSCE) sce = nullptr;

  for (auto const &nt: fNormTools)
    nt->setup(evt);
  
  // get the PFParticle's and the associated data
  art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
  evt.getByLabel(fPFParticleLabel, pfparticle_handle);

  std::vector<art::Ptr<recob::PFParticle>> pfparticles;
  art::fill_ptr_vector(pfparticles, pfparticle_handle);

  art::FindManyP<recob::Vertex> pfparticleVertices(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::Cluster> pfparticleClusters(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::Slice> pfparticleSlices(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::Track> pfparticleTracks(pfparticles, evt, fTrackLabel);

  // organize the PFPlist into a map
  std::map<unsigned, art::Ptr<recob::PFParticle>> id_to_pfp;
  for (unsigned i = 0; i < pfparticles.size(); i++) {
    id_to_pfp[pfparticles[i]->Self()] = pfparticles[i];
  }

  // get all the hits
  // art::Handle<std::vector<recob::Hit>> hit_handle;
  // evt.getByLabel(fHitLabel, hit_handle);

  // std::vector<art::Ptr<recob::Hit>> hits;
  // art::fill_ptr_vector(hits, hit_handle);

  // map to space-points
  // art::FindManyP<recob::SpacePoint> hitSPs(hits, evt, fPFParticleLabel);

  // iterate over the "primary" PFParticles
  for (unsigned i_pfp = 0; i_pfp < pfparticles.size(); i_pfp++) {
    const recob::PFParticle &pfp = *pfparticles[i_pfp];
    if (!pfp.IsPrimary()) continue;
    // Ignore PFP's with no vertex
    if (!pfparticleVertices.at(i_pfp).size()) continue;

    // If configured, require this to be a PFP neutrino
    unsigned pfpPDGC = std::abs(pfp.PdgCode());
    if(fSelectNeutrino &&
        (pfpPDGC != 12) && (pfpPDGC != 14) && (pfpPDGC != 16) ) continue;

    // we found a primary PFP! Get its vertex.
    const art::Ptr<recob::Vertex> &vtx_ptr = pfparticleVertices.at(i_pfp).at(0); 
    const recob::Vertex &vert = *vtx_ptr;

    // The presence of space charge creates two different positions: the position as
    // "seen" by the wires (post-space charge) and the true position (pre-space charge)
    //
    // In some cases we want the wire-position and in other we want the true-position.
    // To be explicit, create two different vertexes for these two different cases
    recob::Vertex vert_absolute(PositionAbsolute(vert.position(), geo, sce),
            vert.covariance(), 
            vert.chi2(),
            vert.ndof(),
            vert.ID());

    recob::Vertex vert_atwires(PositionAtWires(vert.position(), geo, sce),
            vert.covariance(), 
            vert.chi2(),
            vert.ndof(),
            vert.ID());

    // also get all the daughter PFParticles
    const std::vector<size_t> &daughters = pfp.Daughters();
    std::vector<art::Ptr<recob::PFParticle>> daughterPFPs;
    for (size_t d: daughters) {
      daughterPFPs.push_back(id_to_pfp.at(d));
    }

    // look up the hits of each daughter
    std::array<std::vector<std::vector<art::Ptr<recob::Hit>>>, 3> daughterPlaneHits;
    for (const art::Ptr<recob::PFParticle> &d: daughterPFPs) {
      for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
        daughterPlaneHits[i_plane].emplace_back();
      }
      const std::vector<art::Ptr<recob::Cluster>> &d_clusters = pfparticleClusters.at(d.key());
      art::FindManyP<recob::Hit> d_cluster_hits(d_clusters, evt, fPFParticleLabel);
      for (unsigned i = 0; i < d_clusters.size(); i++) {
        const std::vector<art::Ptr<recob::Hit>> &this_cluster_hits = d_cluster_hits.at(i);
        daughterPlaneHits[d_clusters[i]->Plane().Plane].back().insert(
          daughterPlaneHits[d_clusters[i]->Plane().Plane].back().end(),
          this_cluster_hits.begin(), this_cluster_hits.end());
      }
    }

    // Get the Slice associated with the primary PFP
    art::Ptr<recob::Slice> thisSlc = pfparticleSlices.at(i_pfp).at(0);
    // look up the hits
    art::FindManyP<recob::Hit> thisSlcHits({thisSlc}, evt, fPFParticleLabel);
    const std::vector<art::Ptr<recob::Hit>> &hits = thisSlcHits.at(0);

    // work on each plane
    for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
      // vacuum up all the hits within the radius
      std::vector<art::Ptr<recob::Hit>> nearbyHits;
      for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
        const recob::Hit &hit = *hits[i_hit];
        if (hit.WireID().Plane == i_plane) {
          if (Vert2HitDistance(hit, vert_atwires, geo, dprop) < fHitVacuumRadius) {
            nearbyHits.push_back(hits[i_hit]);
          }
        }
      }

      // and find the hit SP's
      art::FindManyP<recob::SpacePoint> hitSPs(nearbyHits, evt, fPFParticleLabel);

      // Compute all needed information for each hit
      for (unsigned i_hit = 0; i_hit < nearbyHits.size(); i_hit++) {
        const recob::Hit &hit = *nearbyHits[i_hit];

	sbn::VertexHit vhit;
        vhit.wire = hit.WireID();
	vhit.proj_dist_to_vertex = Vert2HitDistance(hit, vert_atwires, geo, dprop);
        vhit.vtxw = geo->WireCoordinate(vert_atwires.position(), hit.WireID());
        vhit.vtxx = vert_atwires.position().x();
        vhit.vtxXYZ = vert_absolute.position();

        // Compute the charge, using everything we have
        vhit.charge = fCaloAlg.ElectronsFromADCArea(Normalize(hit.Integral(), evt, hit, 
                                                              vert_absolute.position() /* close enuff to hit */, 
                                                              geo::Vector_t() /* no direction available */,
                                                              0. /* TODO: add T0*/),  
                                   hit.WireID().Plane);


        // lookup the spacepoint location
        const std::vector<art::Ptr<recob::SpacePoint>> &hit_sp = hitSPs.at(i_hit);

        bool has_xyz = false;
        int spID = -1;
        geo::Point_t spXYZ;
        // Space-Point!
        if (hit_sp.size()) {
          const recob::SpacePoint sp = *hit_sp.at(0);
          spID = sp.ID();
          spXYZ = PositionAbsolute(sp.position(), geo, sce);
          has_xyz = true;
        }
        // No Space-Point. If configured, see if we can look up a point along the assigned track
        else if (fUseTrackSPRecovery) {
          unsigned plane = hit.WireID().Plane;
          art::Ptr<recob::PFParticle> matchingPFP;
          for (unsigned i_pfp_chk = 0; i_pfp_chk < daughterPlaneHits[plane].size(); i_pfp_chk++) {
            for (unsigned i_hit_chk = 0; i_hit_chk < daughterPlaneHits[plane][i_pfp_chk].size(); i_hit_chk++) {
              if (daughterPlaneHits[plane][i_pfp_chk][i_hit_chk] == nearbyHits[i_hit]) {
                matchingPFP = daughterPFPs[i_pfp_chk];
                break;
              }
            }
          }
          if (matchingPFP) {
            const std::vector<art::Ptr<recob::Track>> &pfptrack = pfparticleTracks.at(matchingPFP.key()); 
            if (pfptrack.size()) {
              const recob::Track &thisTrack = *pfptrack.at(0);
              spXYZ = PositionAbsolute(PlaceHitAlongTrack(thisTrack, vert_atwires, hit, geo, dprop), geo, sce); 
              has_xyz = true;
            }
          }
        }

        if (has_xyz) {
          vhit.spID = spID;
          vhit.spXYZ = spXYZ;

          geo::Vector_t dir = (spXYZ - vert_absolute.position()).Unit();

          // Compute the pitch. Since we have used the corrected vertex and Space-Point position, the
          // pt and dir here are space-charge corrected regardless of the input configuration
          vhit.pitch = sbn::GetPitch(geo, sce, spXYZ, dir, hit.View(), hit.WireID(), fCorrectSCE, true); 

          vhit.dqdx = vhit.charge / vhit.pitch;

          // Same here -- input position is already corrected
          float EField = sbn::GetEfield(dprop, sce, spXYZ, hit.WireID(), false);

          vhit.dedx = fCaloAlg.dEdx_AREA(clock_data, dprop, vhit.dqdx, hit.PeakTime(), hit.WireID().Plane, 0., EField);
        }
        else {
          vhit.spID = -1;
          vhit.pitch = -1;
          vhit.dqdx = -1;
          vhit.dedx = -1;
        }

        // Save!
        outVHit->push_back(vhit);
        art::Ptr<sbn::VertexHit> thisVHitPtr = vhitPtrMaker(outVHit->size()-1);
        assn->addSingle(thisSlc, thisVHitPtr);
        vtxAssn->addSingle(vtx_ptr, thisVHitPtr);
        hitAssn->addSingle(nearbyHits[i_hit], thisVHitPtr);

      } // end iterate over hits
    } // end iterate over planes
  } // end iterate over pfparticle's 

  evt.put(std::move(outVHit));
  evt.put(std::move(vtxAssn));
  evt.put(std::move(assn));
  evt.put(std::move(hitAssn));

}

DEFINE_ART_MODULE(sbn::VertexChargeVacuum)
