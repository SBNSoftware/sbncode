////////////////////////////////////////////////////////////////////////
// Class:       PCAnglePlaneMakerFinder
// Plugin Type: producer (art v3_02_06)
// File:        PCAnglePlaneMaker_module.cc
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

#include "lardataobj/RecoBase/PFParticle.h"

#include <memory>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "Products/PCAnglePlane.h"

namespace sbn {
  class PCAnglePlaneMaker;

  // struct for organizing each angle info
  struct PCAngleInfo {
    sbn::PCAngle angle;
    unsigned generation;
    unsigned branch;
    unsigned hit_ind;
  };
}


class sbn::PCAnglePlaneMaker : public art::EDProducer {
public:
  explicit PCAnglePlaneMaker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PCAnglePlaneMaker(PCAnglePlaneMaker const&) = delete;
  PCAnglePlaneMaker(PCAnglePlaneMaker&&) = delete;
  PCAnglePlaneMaker& operator=(PCAnglePlaneMaker const&) = delete;
  PCAnglePlaneMaker& operator=(PCAnglePlaneMaker&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;


private:
  // private data 
  art::InputTag fPFParticleLabel;
  bool fFollowDaughters;
  float fHitGroupDistance;
  bool fOnlyPrimary;

};


// static helper functions
void SaveHits(std::map<unsigned, std::array<std::vector<unsigned>, 3>> &pfpToHits, const std::vector<art::Ptr<recob::Hit>> &hits, unsigned plane, const art::Ptr<recob::PFParticle> &pfp) {
  for (const art::Ptr<recob::Hit> &h: hits) {
    pfpToHits[pfp.key()][plane].push_back(h.key());
  }
}

float VecAngle(std::array<float, 2> A, std::array<float, 2> B) {
  float costh = (A[0] * B[0] + A[1] * B[1]) \
    / (sqrt(A[0]*A[0] + A[1] * A[1]) * sqrt(B[0]*B[0] + B[1] * B[1]));

  return acos(costh);
}

std::array<float, 2> HitVector(const recob::Hit &A, const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {
  // get the wire distance between A and B
  float wire_distance = A.WireID().Wire;
  // convert to cm
  float wire_distance_cm = wire_distance * geo->WirePitch();

  // and the time difference
  float time_distance = A.PeakTime();
  // convert to cm
  float time_distance_cm = dprop->ConvertTicksToX(time_distance, A.WireID());

  return {wire_distance_cm, time_distance_cm};
}

std::array<float, 2> HitVector(const recob::Hit &A, const recob::Hit &B, const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {
  // get each individual vec
  std::array<float, 2> vecA = HitVector(A, geo, dprop);
  std::array<float, 2> vecB = HitVector(B, geo, dprop);

  return {vecA[0] - vecB[0], vecA[1] - vecB[1]};
}

float HitDistance(const recob::Hit &A, const recob::Hit &B, const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {
  std::array<float, 2> vec = HitVector(A, B, geo, dprop);
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
}

/*
std::tuple<std::vector<art::Ptr<recob::Hit>>, std::vector<art::Ptr<recob::Hit>>, bool> GetNearestHits(
                                                 const std::vector<art::Ptr<recob::Hit>> &hits, int ihit, float distance,
                                                 const geo::GeometryCore *geo,
                                                 const detinfo::DetectorProperties *dprop) {

  std::vector<art::Ptr<recob::Hit>> retlo;
  std::vector<art::Ptr<recob::Hit>> rethi;


  std::vector<art::Ptr<recob::Hit>> nearby;
  // pull in all the hits
  for (unsigned i = 0; i < hits.size(); i++) {
    if (HitDistance(*hits[i], *hits[ihit], geo, dprop) < distance) {
      nearby.push_back(hits[i]);
    }
  } 

  art::Ptr<recob::Hit> center = hits[ihit];
  // run a PCA to find their general direction
  std::array<float, 2> dir = HitPCAVec(nearby, center, geo, dprop);

  // order the hits according to this direction
  std::sort(nearby.begin(), nearby.end(),
    [center, geo, dprop](auto const &lhs, auto const &rhs) {
    std::array<float, 2> lhs_vec = HitVector(*lhs, *center, geo, dprop);
    std::array<float, 2> rhs_vec = HitVector(*rhs, *center, geo, dprop);
    float proj_lhs = dir[0] * lhs_vec[0] + dir[1] + lhs_vec[1];
    float proj_rhs = dir[0] * rhs_vec[0] + dir[1] + rhs_vec[1];
    return proj_lhs < proj_rhs;
  });

  // everything up to the center hit is "lo"
  bool found_center = false;
  for (unsigned i = 0; i < nearby.size(); i++) {
    
  }

  return {retlo, rethi, true}; 
}*/

std::tuple<std::vector<art::Ptr<recob::Hit>>, std::vector<art::Ptr<recob::Hit>>, bool> GetNearestHits(
                                                 const std::vector<art::Ptr<recob::Hit>> &hits, int ihit, float distance,
                                                 const geo::GeometryCore *geo,
                                                 const detinfo::DetectorProperties *dprop) {
  std::vector<art::Ptr<recob::Hit>> retlo;
  std::vector<art::Ptr<recob::Hit>> rethi;

  bool lo_complete = false;
  // pull in smaller ones
  for (int j = ihit-1; j >= 0; j--) {
    if (HitDistance(*hits[j], *hits[ihit], geo, dprop) > distance) {
      lo_complete = true;
      break;
    }
    retlo.push_back(hits[j]);
  }

  bool hi_complete = false;
  // pull in larger ones
  for (unsigned j = ihit+1; j < hits.size(); j++) {
    if (HitDistance(*hits[j], *hits[ihit], geo, dprop) > distance) {
      hi_complete = true;
      break;
    }
    rethi.push_back(hits[j]);
  }

  return {retlo, rethi, lo_complete && hi_complete};
}

std::array<float, 2> HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Ptr<recob::Hit> &center,
               const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {

  if (hits.size() < 2) return {-100., -100.};

  std::array<float, 2> sum {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, *center, geo, dprop);
    sum[0] += vec[0];
    sum[1] += vec[1];
  }
  sum[0] = sum[0] / hits.size();
  sum[1] = sum[1] / hits.size();

  std::array<std::array<float, 2>, 2> scatter {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, *center, geo, dprop);
    vec[0] -= sum[0];
    vec[1] -= sum[1];

    scatter[0][0] += vec[0] * vec[0];
    scatter[0][1] += vec[0] * vec[1];
    scatter[1][0] += vec[1] * vec[0];
    scatter[1][1] += vec[1] * vec[1];
  }

  // first get the eigenvalues of the matrix
  float trace = scatter[0][0] + scatter[1][1];
  float det = scatter[0][0] * scatter[1][1] - scatter[0][1] * scatter[1][0];

  // this is always the max-eigenvalue
  float eigenP = (1. / 2.) * (trace + sqrt(trace*trace - 4 * det));
  // float eigenM = (1. / 2.) * (trace - sqrt(trace*trace - 4 * det));

  // and then the eigenvectors
  std::array<float, 2> ret {scatter[0][1], eigenP - scatter[0][0]};
  // std::array<float, 2> eigenVM {scatter[0][1], eigenM - scatter[0][0]};

  // make sure the sign is right
  if (sum[0] * ret[0] + sum[1] * ret[1] < 0.) {
    ret[0] = -ret[0];
    ret[1] = -ret[1];
  }

  return ret;

}

std::array<float, 2> HitPCAEigen(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Ptr<recob::Hit> &center,
               const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {

  std::array<float, 2> sum {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, *center, geo, dprop);
    sum[0] += vec[0];
    sum[1] += vec[1];
  }
  sum[0] = sum[0] / hits.size();
  sum[1] = sum[1] / hits.size();

  std::array<std::array<float, 2>, 2> scatter {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, *center, geo, dprop);
    vec[0] -= sum[0];
    vec[1] -= sum[1];

    scatter[0][0] += vec[0] * vec[0];
    scatter[0][1] += vec[0] * vec[1];
    scatter[1][0] += vec[1] * vec[0];
    scatter[1][1] += vec[1] * vec[1];
  }

  // first get the eigenvalues of the matrix
  float trace = scatter[0][0] + scatter[1][1];
  float det = scatter[0][0] * scatter[1][1] - scatter[0][1] * scatter[1][0];

  float eigenP = (1. / 2.) * (trace + sqrt(trace*trace - 4 * det));
  float eigenM = (1. / 2.) * (trace - sqrt(trace*trace - 4 * det));

  return {eigenP, eigenM};
}

float Vert2HitDistance(const recob::Hit &hit, const recob::Vertex &vert, const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {
  float vert_wire_coord = geo->WireCoordinate(vert.position().Y(), vert.position().Z(), hit.WireID()) * geo->WirePitch();
  float hit_wire_coord = hit.WireID().Wire * geo->WirePitch();

  float vert_time_coord = vert.position().X();
  float hit_time_coord = dprop->ConvertTicksToX(hit.PeakTime(), hit.WireID());

  return sqrt((vert_wire_coord - hit_wire_coord) * (vert_wire_coord - hit_wire_coord) +
    (vert_time_coord - hit_time_coord) * (vert_time_coord - hit_time_coord));
}

std::array<std::vector<art::Ptr<recob::Hit>>, 3> SortHits(const std::array<std::vector<art::Ptr<recob::Hit>>, 3> &hits, 
                                                          const recob::Vertex &start,
                                                          const geo::GeometryCore *geo, 
                                                          const detinfo::DetectorProperties *dprop) {
  // convert data to linked-list
  std::array<std::list<art::Ptr<recob::Hit>>, 3> planar {};
  for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
    planar[i_plane].insert(planar[i_plane].begin(), hits[i_plane].begin(), hits[i_plane].end());
  }

  // start from the zeroth hit and grow from the nearest hit
  std::array<std::vector<art::Ptr<recob::Hit>>, 3> ret {};
  for (unsigned i = 0; i < 3; i++) {
    if (planar[i].size()) {

      // find the zeroth hit
      std::list<art::Ptr<recob::Hit>>::iterator firstHit = std::min_element(planar[i].begin(), planar[i].end(),
          [start, geo, dprop](auto const &lhs, auto const &rhs) {
            return Vert2HitDistance(*lhs, start, geo, dprop) < Vert2HitDistance(*rhs, start, geo, dprop);
      });
      ret[i].push_back(*firstHit);
      planar[i].erase(firstHit);

      while (planar[i].size()) {
        const recob::Hit &lastHit = *ret[i].back();
        std::list<art::Ptr<recob::Hit>>::iterator closest = std::min_element(planar[i].begin(), planar[i].end(),
          [lastHit, geo, dprop](auto const &lhs, auto const &rhs) { 
            return HitDistance(*lhs, lastHit, geo, dprop) < HitDistance(*rhs, lastHit, geo, dprop);
          });

        ret[i].push_back(*closest);
        planar[i].erase(closest);
      }
    }
  }

  return ret;

}

// TODO: implement
//
// For now, always returning 'true' is ok -- it just means we do more work than necessary
bool DoBranch(art::Ptr<recob::PFParticle> particle, const std::map<unsigned, art::Ptr<recob::PFParticle>> &pfpMap) {
  (void) particle;
  (void) pfpMap;

  return true;
}

// remove duplicates induced by degenerate-branching and organize by branch ID
std::map<unsigned, std::vector<sbn::PCAngleInfo>> RemoveDupes(std::vector<sbn::PCAngleInfo> &angles) {
  // We will eliminate duplicates by converting to a set and then back to a vector
  //
  // Define a struct to implement the comparator
  struct angle_compare {
    bool operator() (const sbn::PCAngleInfo &lhs, const sbn::PCAngleInfo &rhs) const {
      // unique ID is hitID, hitIDLo, hitIDHi
      if (lhs.angle.hitID != rhs.angle.hitID) return lhs.angle.hitID < rhs.angle.hitID;
      if (lhs.angle.hitIDLo != rhs.angle.hitIDLo) return lhs.angle.hitIDLo < rhs.angle.hitIDLo;
      if (lhs.angle.hitIDHi != rhs.angle.hitIDHi) return lhs.angle.hitIDHi < rhs.angle.hitIDHi;
      // equal
      return false;
    }
  };

  // Given the choice, we want to keep angles with a smaller generation. 
  //
  // So first sort the input by generation (lo -> hi). So smaller generations
  // get inserted into the set first
  std::sort(angles.begin(), angles.end(),
     [](auto const &lhs, auto const &rhs) {return lhs.generation < rhs.generation;});

  // make the set to unique-ify the input
  std::set<sbn::PCAngleInfo, angle_compare> unique(angles.begin(), angles.end());

  // organize into output
  std::map<unsigned, std::vector<sbn::PCAngleInfo>> ret;
  for (const sbn::PCAngleInfo &angle: unique) {
    ret[angle.branch].push_back(angle);
  }

  // sort each ret by the hit order
  for (auto &pair: ret) {
    std::sort(pair.second.begin(), pair.second.end(),
      [](auto const &lhs, auto const &rhs) {
        return lhs.hit_ind < rhs.hit_ind;
   });
  }

  return ret;
}

sbn::PCAnglePlaneMaker::PCAnglePlaneMaker(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabel(p.get<std::string>("PFParticleLabel", "pandora")),
    fFollowDaughters(p.get<bool>("FollowDaughters", true)),
    fHitGroupDistance(p.get<float>("HitGroupDistance")),
    fOnlyPrimary(p.get<bool>("OnlyPrimary", true))
{
  produces<std::vector<sbn::PCAnglePlane>>();
  produces<art::Assns<recob::PFParticle, sbn::PCAnglePlane>>();
}

void sbn::PCAnglePlaneMaker::produce(art::Event& evt)
{
  // output stuff
  std::unique_ptr<std::vector<sbn::PCAnglePlane>> outAngles(new std::vector<sbn::PCAnglePlane>);
  std::unique_ptr<art::Assns<recob::PFParticle, sbn::PCAnglePlane>> assn(new art::Assns<recob::PFParticle, sbn::PCAnglePlane>);

  art::PtrMaker<sbn::PCAnglePlane> anglePtrMaker {evt};

  // collect services
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();
  detinfo::DetectorProperties const* dprop
          = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // get the PFParticle's and the associated data
  art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
  evt.getByLabel(fPFParticleLabel, pfparticle_handle);

  std::vector<art::Ptr<recob::PFParticle>> pfparticles;
  art::fill_ptr_vector(pfparticles, pfparticle_handle);

  // organize the PFPlist into a map
  std::map<unsigned, art::Ptr<recob::PFParticle>> id_to_pfp;
  for (unsigned i = 0; i < pfparticles.size(); i++) {
    id_to_pfp[pfparticles[i]->Self()] = pfparticles[i];
  }

  art::FindManyP<recob::Cluster> pfparticleClusters(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::Vertex> pfparticleVertices(pfparticles, evt, fPFParticleLabel);

  // iterate over each particle
  for (unsigned i_part = 0; i_part < pfparticles.size(); i_part++) {
    std::cout << "New particle! " << i_part << std::endl;
    art::Ptr<recob::PFParticle> thisPFP = pfparticles[i_part];

    std::cout << "Child of: " << thisPFP->Parent() << std::endl;

    // ignore non-primary or child of primary
    if (fOnlyPrimary &&
        thisPFP->Parent() != recob::PFParticle::kPFParticlePrimary &&
        id_to_pfp.at(thisPFP->Parent())->Parent() != recob::PFParticle::kPFParticlePrimary) {
      continue;
    }

    // keep track of the map between this PFParticle and all the hits
    std::map<unsigned, std::array<std::vector<unsigned>, 3>> pfpToHits;
    std::map<unsigned, std::array<std::vector<art::Ptr<recob::Hit>>,3>> allHits;
      
    const std::vector<art::Ptr<recob::Cluster>> &thisCluster = pfparticleClusters.at(i_part);
    art::FindManyP<recob::Hit> clusterHits(thisCluster, evt, fPFParticleLabel);

    // ignore PFP's w/out hits
    if (!thisCluster.size()) {
      std::cout << "No clusters :/\n";
      continue;
    }

    const std::vector<art::Ptr<recob::Vertex>> &pfpVerts = pfparticleVertices.at(i_part);

    // ignore PFP's w/out vertex
    if (!pfpVerts.size()) {
      std::cout << "No vertex :/\n";
      continue;
    }

    // grab the vertex
    const recob::Vertex &pfpVert = *pfpVerts.at(0);

    std::cout << "Valid!\n";

    for (unsigned i_clus = 0; i_clus < thisCluster.size(); i_clus++) {
      allHits[thisPFP->Self()][thisCluster[i_clus]->Plane().Plane].insert(allHits[thisPFP->Self()][thisCluster[i_clus]->Plane().Plane].begin(),
          clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
      SaveHits(pfpToHits, clusterHits.at(i_clus), thisCluster[i_clus]->Plane().Plane, thisPFP);
    }

    // also get the hits from each daughter PFP (and so on) if configured
    if (fFollowDaughters) {
      // define the data we keep in the stack
      struct daughter_stack_info {
        unsigned daughter_id;
        long unsigned hit_key;
      };

      std::stack<daughter_stack_info> d_pfps;
      for (unsigned d: thisPFP->Daughters()) {
        daughter_stack_info dsi {d, thisPFP->Self()};
        d_pfps.push(dsi);
      }
        
      // DFS
      while (!d_pfps.empty()) {
        daughter_stack_info daughter_info = d_pfps.top();
        unsigned daughter = daughter_info.daughter_id;
        unsigned hit_key = daughter_info.hit_key;
        d_pfps.pop();

        // get the PFP
        const art::Ptr<recob::PFParticle> &d_pfp = id_to_pfp.at(daughter);

        // determine whether to branch
        if (DoBranch(d_pfp, id_to_pfp)) {
          // branching -- make a new entry in the allHits-map with this particle
          //
          // First clone the hits we already have
          allHits[d_pfp->Self()][0] = allHits[hit_key][0];
          allHits[d_pfp->Self()][1] = allHits[hit_key][1];
          allHits[d_pfp->Self()][2] = allHits[hit_key][2];

          // and re-direct this hits to the new location
          hit_key = d_pfp->Self();
        }

        // add the daughter-daughters
        for (unsigned d: d_pfp->Daughters()) {
          daughter_stack_info dsi {d, hit_key};
          d_pfps.push(dsi);
        }

        // add the hits 
        const std::vector<art::Ptr<recob::Cluster>> &d_cluster = pfparticleClusters.at(d_pfp.key());

        // ignore PFP's w/out hits
        if (!d_cluster.size()) continue;

        art::FindManyP<recob::Hit> d_cluster_hits(d_cluster, evt, fPFParticleLabel);
        for (unsigned i_clus = 0; i_clus < d_cluster.size(); i_clus++) {
          allHits[hit_key][d_cluster[i_clus]->Plane().Plane].insert(allHits[hit_key][d_cluster[i_clus]->Plane().Plane].begin(),
            d_cluster_hits.at(i_clus).begin(), d_cluster_hits.at(i_clus).end());
          SaveHits(pfpToHits, d_cluster_hits.at(i_clus), d_cluster[i_clus]->Plane().Plane, d_pfp);
        }
      }
    }


    std::array<std::vector<sbn::PCAngleInfo>, 3> planeAngles;
    // Now we have all the hits we need!
    //
    // iterate through each branch
    for (auto const &pair: allHits) {
      unsigned branch_id = pair.first;
      const std::array<std::vector<art::Ptr<recob::Hit>>, 3> &branchHits = pair.second;

      // Now, for each hit, we'll consider the PCA angle and the vec angle
      //
      // First, order this hits
      std::array<std::vector<art::Ptr<recob::Hit>>, 3> sortedHits = SortHits(branchHits, pfpVert, geo, dprop);

      // iterate through all the Hits
      for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
        for (unsigned i_hit = 0; i_hit < sortedHits[i_plane].size(); i_hit++) {
          // get nearby hits
          auto [hitslo, hitshi, complete] = GetNearestHits(sortedHits[i_plane], i_hit, fHitGroupDistance, geo, dprop);

          //std::cout << "Plane: " << i_plane << " Branch: " << branch_id << " hit: " << i_hit << " " << sortedHits[i_plane][i_hit].key();
          //if (hitslo.size()) std::cout << " lo: " << hitslo.back().key();
          //else std::cout << " lo: none";
          //if (hitshi.size()) std::cout << " hi: " << hitshi.back().key();
          //else std::cout << " hi: none";
          //std::cout << std::endl;

          // invalid -- continue
          // if (hithi.size() < 2 || hitslo.size() < 2) continue;

          // Get the PCA dirs
          std::array<float, 2> pca_vec_lo = HitPCAVec(hitslo, sortedHits[i_plane][i_hit], geo, dprop);
          std::array<float, 2> pca_vec_hi = HitPCAVec(hitshi, sortedHits[i_plane][i_hit], geo, dprop);
          float angle = -100.;
          // check if valid vecs
          if (pca_vec_lo[0] > -99. && pca_vec_lo[1] > -99. && pca_vec_hi[0] > -99. && pca_vec_hi[1] > -99.) {
            // get angle between the two vecs and flip
            angle = M_PI - VecAngle(pca_vec_lo, pca_vec_hi);
          }

          std::array<float, 2> hv = HitVector(*sortedHits[i_plane][i_hit], geo, dprop);

          // and save
          sbn::PCAngleInfo thisAngle;
          thisAngle.angle.wire = sortedHits[i_plane][i_hit]->WireID();
          thisAngle.angle.angle = angle;
          thisAngle.angle.lo_vector = pca_vec_lo;
          thisAngle.angle.hi_vector = pca_vec_hi;
          thisAngle.angle.hit_wire_cm = hv[0];
          thisAngle.angle.hit_time_cm = hv[1];
          thisAngle.angle.complete = complete;
          thisAngle.angle.hitID = sortedHits[i_plane][i_hit].key();

          if (hitshi.size()) thisAngle.angle.dist_to_hi = HitDistance(*hitshi.back(), *sortedHits[i_plane][i_hit], geo, dprop); 
          else thisAngle.angle.dist_to_hi = -100.;
          if (hitslo.size()) thisAngle.angle.dist_to_lo = HitDistance(*hitslo.back(), *sortedHits[i_plane][i_hit], geo, dprop); 
          else thisAngle.angle.dist_to_lo = -100.;

          if (hitshi.size()) thisAngle.angle.hitIDHi = hitshi.back().key();
          else thisAngle.angle.hitIDHi = -1;
          if (hitslo.size()) thisAngle.angle.hitIDLo = hitslo.back().key();
          else thisAngle.angle.hitIDLo = -1;

          thisAngle.branch = branch_id;
          thisAngle.hit_ind = i_hit;

          // get the generation
          thisAngle.generation = 0;
          art::Ptr<recob::PFParticle> check = id_to_pfp.at(branch_id);
          while (check->Self() != thisPFP->Self()) {
            thisAngle.generation ++;
            check = id_to_pfp.at(check->Parent());
          }

          planeAngles[i_plane].push_back(thisAngle);
        } 
      }
    }

    // Save the angles on each plane
    for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
      std::map<unsigned, std::vector<sbn::PCAngleInfo>> selected = RemoveDupes(planeAngles[i_plane]);

      sbn::PCAnglePlane thisPlane;
      geo::PlaneID planeID(0, 0, i_plane); // TPC and cryo are arbitrary
      thisPlane.plane = planeID;
      for (auto const &pair: selected) {
        thisPlane.branchIDs.push_back(pair.first);

        // look up the hierarchy
        thisPlane.branchHierarchy.emplace_back();
        art::Ptr<recob::PFParticle> check = id_to_pfp.at(pair.first);
        thisPlane.branchHierarchy.back().push_back(check->Self());
        while (check->Self() != thisPFP->Self()) {
          check = id_to_pfp.at(check->Parent());
          thisPlane.branchHierarchy.back().push_back(check->Self());
        }

        if (pair.second.size()) thisPlane.generations.push_back(pair.second.front().generation);
        else thisPlane.generations.push_back(-1);

        if (pair.second.size()) thisPlane.plane = pair.second.front().angle.wire;

        thisPlane.angles.emplace_back();
        for (const sbn::PCAngleInfo &a: pair.second) thisPlane.angles.back().push_back(a.angle);
      }
      thisPlane.nBranches = thisPlane.angles.size();
      outAngles->push_back(thisPlane);
      art::Ptr<sbn::PCAnglePlane> thisAnglePtr = anglePtrMaker(outAngles->size()-1);
      assn->addSingle(thisPFP, thisAnglePtr);
    }

  }

  evt.put(std::move(outAngles));
  evt.put(std::move(assn));

}

DEFINE_ART_MODULE(sbn::PCAnglePlaneMaker)
