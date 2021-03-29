////////////////////////////////////////////////////////////////////////
// Class:       VertexChargeVacuumFinder
// Plugin Type: producer (art v3_02_06)
// File:        VertexChargeVacuum_module.cc
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
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

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

  // private data 
  calo::CalorimetryAlg fCaloAlg;
};


sbn::VertexChargeVacuum::VertexChargeVacuum(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabel(p.get<std::string>("PFParticleLabel", "pandora")),
    fTrackLabel(p.get<std::string>("TrackLabel", "pandoraTrack")),
    fHitVacuumRadius(p.get<float>("HitVacuumRadius")),
    fUseTrackSPRecovery(p.get<bool>("UseTrackSPRecovery")),
    fCaloAlg(p.get<fhicl::ParameterSet >("CaloAlg"))
{

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

  TVector3 avg(0., 0., 0.);

  std::array<float, 2> vert_v = VertexVector(vert, hits[0]->WireID(), geo, dprop);
  TVector3 vert_p(vert_v[0], vert_v[1], 0.);

  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> hit_v = HitVector(*h, geo, dprop);
    TVector3 hit_p(hit_v[0], hit_v[1], 0.);
    avg += (hit_p - vert_p).Unit(); 
  }

  avg = avg.Unit();

  return {(float)avg.X(), (float)avg.Y()};
}

std::array<float, 2> TrackDirection(const recob::Track &trk, const geo::PlaneID &plane,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {

  float angleToVert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();
  float cosgamma = std::abs(std::sin(angleToVert)*trk.StartDirection().y() + std::cos(angleToVert)*trk.StartDirection().z());

  std::array<float, 2> ret {cosgamma, (float)trk.StartDirection().x()};

  // normalize
  float ret_norm = sqrt(ret[0] * ret[0] + ret[1] * ret[1]);

  return {ret[0] / ret_norm, ret[1] / ret_norm};
}

float TrackDirectionPerp(const recob::Track &trk, const geo::PlaneID &plane,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  double angleToVert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();
  double cosgamma = std::abs(std::sin(angleToVert)*trk.StartDirection().y() + std::cos(angleToVert)*trk.StartDirection().z());

  float ret = sqrt(1 - cosgamma * cosgamma - trk.StartDirection().x() * trk.StartDirection().x());

  return ret;
}

TVector3 PlaceHitAlongTrack(const recob::Track &trk, const recob::Vertex &vert, const recob::Hit &hit,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  // project the vertex onto the hit Plane
  std::array<float, 2> v_plane = VertexVector(vert, hit.WireID(), geo, dprop); 

  // also get the 2d hit location
  std::array<float, 2> h_plane = HitVector(hit, geo, dprop);

  // Project the track direction onto the hit Plane
  std::array<float, 2> t_dir_plane = TrackDirection(trk, hit.WireID(), geo, dprop);

  float t_dir_perp_plane = TrackDirectionPerp(trk, hit.WireID(), geo, dprop);

  // Get the distance from the vertex to the hit on the plane projected along the track direction
  float plane_dist = (h_plane[0] - v_plane[0]) * t_dir_plane[0] + (h_plane[1] - v_plane[1]) * t_dir_plane[1];
 
  // Use this distance to place the hit along the track trajectory in 3D
  float dist_3d = plane_dist / sqrt(1 - t_dir_perp_plane * t_dir_perp_plane);

  // std::cout << "Placing Hit\n";
  // std::cout << "Hit wire: " << hit.WireID() << std::endl;
  // std::cout << "Hit plane coord: " << h_plane[0] << " " << h_plane[1] << std::endl;
  // std::cout << "Vtx plane coord: " << v_plane[0] << " " << v_plane[1] << std::endl;
  // std::cout << "Trk dir plane coord: " << t_dir_plane[0] << " " << t_dir_plane[1] << std::endl;
  // std::cout << "Trk dir plane proj: " << t_dir_plane[0]*(trk.StartDirection().x()/t_dir_plane[1]) << " " << t_dir_plane[1]*(trk.StartDirection().x()/t_dir_plane[1]) << std::endl;
  // std::cout << "Trk dir: " << trk.StartDirection().x() << " " << trk.StartDirection().y() << " " << trk.StartDirection().z() << std::endl;
  // std::cout << "Vtx pos: " << vert.position().x() << " " << vert.position().y() << " " << vert.position().z() << std::endl;

  // std::cout << "Trk perp dir: " << t_dir_perp_plane << std::endl;
  // std::cout << "Hit plane dist: " << plane_dist << std::endl;
  // std::cout << "Hit 3d dist: " << dist_3d << std::endl;

  TVector3 trk_dir(trk.StartDirection().x(), trk.StartDirection().y(), trk.StartDirection().z());
  TVector3 vert_v(vert.position().x(), vert.position().y(), vert.position().z()); 

  return vert_v + trk_dir * dist_3d;
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

  // get the PFParticle's and the associated data
  art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
  evt.getByLabel(fPFParticleLabel, pfparticle_handle);

  std::vector<art::Ptr<recob::PFParticle>> pfparticles;
  art::fill_ptr_vector(pfparticles, pfparticle_handle);

  art::FindManyP<recob::Vertex> pfparticleVertices(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::Cluster> pfparticleClusters(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::SpacePoint> pfparticleSpacePoints(pfparticles, evt, fPFParticleLabel);
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

    // we found a primary PFP! Get its vertex.
    const art::Ptr<recob::Vertex> &vtx_ptr = pfparticleVertices.at(i_pfp).at(0); 
    const recob::Vertex &vert = *vtx_ptr;
    TVector3 vert_v(vert.position().X(), vert.position().Y(), vert.position().Z()); 

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

    // also look up the space-points
    std::vector<int> daughterSPIDs;
    std::vector<TVector3> daughter3DDir;
    for (unsigned i_d = 0; i_d < daughterPFPs.size(); i_d++) {
      const art::Ptr<recob::PFParticle> &d = daughterPFPs[i_d];
      const std::vector<art::Ptr<recob::SpacePoint>> &dsp = pfparticleSpacePoints.at(d.key());
      TVector3 this_dir(0., 0., 0.);
      bool set = false;
      for (const art::Ptr<recob::SpacePoint> &sp: dsp) {
        TVector3 sp_v(sp->XYZ());
        if ((sp_v-vert_v).Mag() < fHitVacuumRadius) {
          this_dir += (sp_v - vert_v).Unit();
          set = true;
        }
      }
      if (set) {
        daughter3DDir.push_back(this_dir.Unit());
        daughterSPIDs.push_back(d->Self());
      }
    }

    // Get the Slice associated with the primary PFP
    art::Ptr<recob::Slice> thisSlc = pfparticleSlices.at(i_pfp).at(0);
    // look up the hits
    art::FindManyP<recob::Hit> thisSlcHits({thisSlc}, evt, fPFParticleLabel);
    const std::vector<art::Ptr<recob::Hit>> &hits = thisSlcHits.at(0);

    // work on each plane
    for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
      std::vector<std::array<float, 2>> daughterPlaneDirs;
      std::vector<int> daughterPFPIDs;
      // Get the PC-axis of each PFParticle daughter projected on this plane
      for (unsigned i_d = 0; i_d < daughterPFPs.size(); i_d++) {
        if (daughterPlaneHits[i_plane][i_d].size()) {
          std::array<float, 2> thisDir = HitDirection(daughterPlaneHits[i_plane][i_d], vert, geo, dprop);
          daughterPlaneDirs.push_back(thisDir);
          daughterPFPIDs.push_back(daughterPFPs[i_d]->Self());
        }
      }

      // vacuum up all the hits within the radius
      std::vector<art::Ptr<recob::Hit>> nearbyHits;
      for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
        const recob::Hit &hit = *hits[i_hit];
        if (hit.WireID().Plane == i_plane) {
          if (Vert2HitDistance(hit, vert, geo, dprop) < fHitVacuumRadius) {
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
        vhit.charge = fCaloAlg.ElectronsFromADCArea(hit.Integral(), hit.WireID().Plane) * fCaloAlg.LifetimeCorrection(clock_data, dprop, hit.PeakTime(), 0.);
	vhit.proj_dist_to_vertex = Vert2HitDistance(hit, vert, geo, dprop);
        vhit.vtxw = geo->WireCoordinate(vert.position(), hit.WireID());

        // lookup the spacepoint location
        const std::vector<art::Ptr<recob::SpacePoint>> &hit_sp = hitSPs.at(i_hit);

        bool has_xyz = false;
        int spID = -1;
        TVector3 spXYZ;
        // Space-Point!
        if (hit_sp.size()) {
          const recob::SpacePoint sp = *hit_sp.at(0);
          spID = sp.ID();
          spXYZ = TVector3(sp.XYZ());
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
            std::cout << "Found PFP: " << matchingPFP->Self() << std::endl;
            const std::vector<art::Ptr<recob::Track>> &pfptrack = pfparticleTracks.at(matchingPFP.key()); 
            std::cout << "PFP has track: " << pfptrack.size() << std::endl;
            if (pfptrack.size()) {
              const recob::Track &thisTrack = *pfptrack.at(0);
              spXYZ = PlaceHitAlongTrack(thisTrack, vert, hit, geo, dprop); 
              has_xyz = true;
            }
          }
        }

        if (has_xyz) {
          vhit.spID = spID;
          vhit.spXYZ = spXYZ;

          geo::Point_t pt (spXYZ);
          geo::Vector_t dir = (pt - vert.position()).Unit();
          vhit.pitch = sbn::GetPitch(geo, sce, pt, dir, hit.View(), hit.WireID(), true, true); 

          vhit.dqdx = fCaloAlg.ElectronsFromADCArea((hit.Integral() / vhit.pitch), hit.WireID().Plane) * fCaloAlg.LifetimeCorrection(clock_data, dprop, hit.PeakTime(), 0.);

          float EField = sbn::GetEfield(dprop, sce, pt, hit.WireID().TPC, true);

          vhit.dedx = fCaloAlg.dEdx_AREA(clock_data, dprop, vhit.dqdx, hit.PeakTime(), hit.WireID().Plane, 0., EField);

          std::vector<std::pair<float, int>> pfp_3d_dist; 
          for (unsigned i_dsp = 0; i_dsp < daughterSPIDs.size(); i_dsp++) {
            pfp_3d_dist.push_back({daughter3DDir[i_dsp].Dot((vhit.spXYZ - vert_v).Unit()), daughterSPIDs[i_dsp]});
          }
          std::sort(pfp_3d_dist.begin(), pfp_3d_dist.end(),
            [](auto const &lhs, auto const &rhs) { return lhs.second > rhs.second;});
          for (auto const &pair: pfp_3d_dist) {
            vhit.nearbyPFP3DDists.push_back(pair.first);
	    vhit.nearbyPFP3DIDs.push_back(pair.second);
          }
        }
        else {
          vhit.spID = -1;
          vhit.pitch = -1;
          vhit.dqdx = -1;
          vhit.dedx = -1;
        }

        // now for each hit compute the dot with each PFP PC
        std::vector<std::pair<float, int>> pfp_proj_dist; 
        for (unsigned i_pfp = 0; i_pfp < daughterPFPIDs.size(); i_pfp++) {
          std::array<float, 2> dir = daughterPlaneDirs[i_pfp];
          std::array<float, 2> vv = VertexVector(vert, hit.WireID(), geo, dprop);
          std::array<float, 2> hv = HitVector(hit, geo, dprop);
          hv[0] = hv[0] - vv[0];
          hv[1] = hv[1] - vv[1];

          float dot;

          if (hv[0] < 1e-4 && hv[1] < 1e-4) dot = 1;
          else dot = (dir[0]*hv[0] + dir[1]*hv[1]) / 
            (sqrt(hv[0]*hv[0]+hv[1]*hv[1]) *
             sqrt(dir[0]*dir[0]+dir[1]*dir[1]));

          pfp_proj_dist.push_back({dot, daughterPFPIDs[i_pfp]});
        }

        std::sort(pfp_proj_dist.begin(), pfp_proj_dist.end(),
          [](auto const &lhs, auto const &rhs) {return lhs.first > rhs.first;});

        for (auto const &pair: pfp_proj_dist) {
          vhit.nearbyPFPDists.push_back(pair.first);
	  vhit.nearbyPFPIDs.push_back(pair.second);
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
