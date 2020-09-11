////////////////////////////////////////////////////////////////////////
// Class:       DThetaKinkFinder
// Plugin Type: producer (art v3_02_06)
// File:        DThetaKinkFinder_module.cc
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

#include "Products/DThetaKink.h"

namespace sbn {
  class DThetaKinkFinder;
}


class sbn::DThetaKinkFinder : public art::EDProducer {
public:
  explicit DThetaKinkFinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DThetaKinkFinder(DThetaKinkFinder const&) = delete;
  DThetaKinkFinder(DThetaKinkFinder&&) = delete;
  DThetaKinkFinder& operator=(DThetaKinkFinder const&) = delete;
  DThetaKinkFinder& operator=(DThetaKinkFinder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;


private:
  // private data 
  art::InputTag fPFParticleLabel;
  bool fFollowDaughters;
  float fHitGroupDistance;
  float fTruncateSigma;
  bool fSaveAllKinks;
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

std::tuple<std::vector<art::Ptr<recob::Hit>>, std::vector<art::Ptr<recob::Hit>>, bool> GetNearestHits(
                                                 const std::vector<art::Ptr<recob::Hit>> &hits, int ihit, float distance,
                                                 const geo::GeometryCore *geo,
                                                 const detinfo::DetectorProperties *dprop) {
  std::vector<art::Ptr<recob::Hit>> retlo;
  std::vector<art::Ptr<recob::Hit>> rethi;

  // pull in smaller ones
  for (int j = ihit-1; j >= 0; j--) {
    if (HitDistance(*hits[j], *hits[ihit], geo, dprop) > distance) {
      break;
    }
    retlo.push_back(hits[j]);
  }
  // complete if the distance to the last hit is at least half the circle
  bool lo_complete = false;
  if (retlo.size()) {
    lo_complete = HitDistance(*retlo.back(), *hits[ihit], geo, dprop) > distance/2.;
  }

  // pull in larger ones
  for (unsigned j = ihit+1; j < hits.size(); j++) {
    if (HitDistance(*hits[j], *hits[ihit], geo, dprop) > distance) {
      break;
    }
    rethi.push_back(hits[j]);
  }
  bool hi_complete = false;
  if (rethi.size()) {
    hi_complete = HitDistance(*rethi.back(), *hits[ihit], geo, dprop) > distance/2.;
  }

  return {retlo, rethi, lo_complete && hi_complete};
}

float HitVecAngle(const std::vector<art::Ptr<recob::Hit>> &hitslo, const std::vector<art::Ptr<recob::Hit>> &hitshi, const art::Ptr<recob::Hit> &center,
               const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {
  if (hitslo.size() == 0 || hitshi.size() == 0) return -100.;

  // get the furthest hit in both directions
  const art::Ptr<recob::Hit> &hlo = hitslo.back();  
  const art::Ptr<recob::Hit> &hhi = hitshi.back();  

  // get the vec to each
  std::array<float, 2> veclo = HitVector(*hlo, *center, geo, dprop);
  std::array<float, 2> vechi = HitVector(*hhi, *center, geo, dprop);

  return VecAngle(veclo, vechi);

}

std::array<float, 2> HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Ptr<recob::Hit> &center,
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

std::vector<art::Ptr<recob::Hit>> TruncateHitsPCA(const std::vector<art::Ptr<recob::Hit>> &allHits, 
                                        const std::vector<art::Ptr<recob::Hit>> &nearbyHits, const art::Ptr<recob::Hit> &center,
                                        const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop, 
                                        float hitDistance, float truncSigma) {

  std::array<float, 2> cvec = HitVector(*center, geo, dprop);
  std::cout << "Checking nearby hits to hit: " << center.key() << " at: " << cvec[0] << " " << cvec[1] << std::endl;

  // get the initial PCA eigen-system
  std::array<float, 2> PCA = HitPCAVec(nearbyHits, center, geo, dprop);
  std::array<float, 2> MCA = {PCA[1], -PCA[0]}; // orthogonal to PCA
  std::array<float, 2> eigen = HitPCAEigen(nearbyHits, center, geo, dprop);
  float eigenHi = eigen[0];
  float eigenLo = eigen[1];

  std::vector<art::Ptr<recob::Hit>> prunedNearbyHits;

  std::cout << "PCA: " << PCA[0] << " " << PCA[1] << std::endl;
  std::cout << "MCA: " << MCA[0] << " " << MCA[1] << std::endl;
  std::cout << "Eigen: " << eigenHi << " " << eigenLo << " " << (eigenLo/eigenHi) << std::endl;

  // now prune the hits if they are too anomalous
  for (unsigned i_hit = 0; i_hit < nearbyHits.size(); i_hit++) {
    std::array<float, 2> cvec = HitVector(*nearbyHits[i_hit], geo, dprop);
    std::cout << "Nearby hit: " << nearbyHits[i_hit].key() << " at: " << cvec[0] << " " << cvec[1] << std::endl;

    auto [hitslo, hitshi, complete] = GetNearestHits(allHits, i_hit, hitDistance, geo, dprop);
    // get the mean component along the PCA and MCA
    float mca_comp = 0.;
    float pca_comp = 0.;
    for (art::Ptr<recob::Hit> hit: hitslo) {
      std::array<float, 2> thisVec = HitVector(*hit, *nearbyHits[i_hit], geo, dprop);
      pca_comp += (thisVec[0] * PCA[0] + thisVec[1] * PCA[1]) * (thisVec[0] * PCA[0] + thisVec[1] * PCA[1]);
      mca_comp += (thisVec[0] * MCA[0] + thisVec[1] * MCA[1]) * (thisVec[0] * MCA[0] + thisVec[1] * MCA[1]);
    }
    for (art::Ptr<recob::Hit> hit: hitshi) {
      std::array<float, 2> thisVec = HitVector(*hit, *nearbyHits[i_hit], geo, dprop);
      pca_comp += (thisVec[0] * PCA[0] + thisVec[1] * PCA[1]) * (thisVec[0] * PCA[0] + thisVec[1] * PCA[1]);
      mca_comp += (thisVec[0] * MCA[0] + thisVec[1] * MCA[1]) * (thisVec[0] * MCA[0] + thisVec[1] * MCA[1]);
    }

    std::cout << "Ratio: " << pca_comp << " " << mca_comp << " " << (mca_comp / pca_comp) << std::endl; 

    // anomalous hits
    if (pca_comp < 1.e-4) continue;
    if (truncSigma * eigenLo / eigenHi < mca_comp / pca_comp) continue;

    prunedNearbyHits.push_back(nearbyHits[i_hit]);
  } 

  return prunedNearbyHits;
}

float HitPCAAngle(const std::vector<art::Ptr<recob::Hit>> &hitslo, const std::vector<art::Ptr<recob::Hit>> &hitshi, const art::Ptr<recob::Hit> &center,
               const geo::GeometryCore *geo, const detinfo::DetectorProperties *dprop) {
  if (hitslo.size() <= 1 || hitshi.size() <= 1) return -100.;

  // get the PCA up the "lo" and "hi" hits
  std::array<float, 2> loPCA = HitPCAVec(hitslo, center, geo, dprop);
  std::array<float, 2> hiPCA = HitPCAVec(hitshi, center, geo, dprop);

  // and then the angle between them
  return VecAngle(loPCA, hiPCA);

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

// remove duplicates induced by degenerate-branching
std::vector<sbn::DThetaKink> SelectKinks(std::vector<sbn::DThetaKink> &kinks) {
  // We will eliminate duplicates by converting to a set and then back to a vector
  //
  // Define a struct to implement the comparator
  struct kink_compare {
    bool operator() (const sbn::DThetaKink &lhs, const sbn::DThetaKink &rhs) const {
      // unique ID is hitID, hitIDLo, hitIDHi
      if (lhs.hitID != rhs.hitID) return lhs.hitID < rhs.hitID;
      if (lhs.hitIDLo != rhs.hitIDLo) return lhs.hitIDLo < rhs.hitIDLo;
      if (lhs.hitIDHi != rhs.hitIDHi) return lhs.hitIDHi < rhs.hitIDHi;
      // equal
      return false;
    }
  };

  // Given the choice, we want to keep kinks with a smaller generation. 
  //
  // So first sort the input by generation (lo -> hi). So smaller generations
  // get inserted into the set first
  std::sort(kinks.begin(), kinks.end(),
     [](auto const &lhs, auto const &rhs) {return lhs.generation < rhs.generation;});

  // make the set to unique-ify the input
  std::set<sbn::DThetaKink, kink_compare> unique(kinks.begin(), kinks.end());

  std::vector<sbn::DThetaKink> ret(unique.begin(), unique.end());

  return ret;
}

sbn::DThetaKinkFinder::DThetaKinkFinder(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabel(p.get<std::string>("PFParticleLabel", "pandora")),
    fFollowDaughters(p.get<bool>("FollowDaughters", true)),
    fHitGroupDistance(p.get<float>("HitGroupDistance")),
    fTruncateSigma(p.get<float>("TruncateSigma")),
    fSaveAllKinks(p.get<bool>("fSaveAllKinks", true)),
    fOnlyPrimary(p.get<bool>("OnlyPrimary", true))
{
  produces<std::vector<sbn::DThetaKink>>();
  produces<art::Assns<recob::PFParticle, sbn::DThetaKink>>();
}

void sbn::DThetaKinkFinder::produce(art::Event& evt)
{
  // output stuff
  std::unique_ptr<std::vector<sbn::DThetaKink>> outKinks(new std::vector<sbn::DThetaKink>);
  std::unique_ptr<art::Assns<recob::PFParticle, sbn::DThetaKink>> assn(new art::Assns<recob::PFParticle, sbn::DThetaKink>);

  art::PtrMaker<sbn::DThetaKink> kinkPtrMaker {evt};

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


    std::array<std::vector<sbn::DThetaKink>, 3> planeKinks;
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

          // Calculate the two angles
          float pca_angle = HitPCAAngle(hitslo, hitshi, sortedHits[i_plane][i_hit], geo, dprop);
          float vec_angle = HitVecAngle(hitslo, hitshi, sortedHits[i_plane][i_hit], geo, dprop);

          // also try the PCA after truncating stuff
          std::vector<art::Ptr<recob::Hit>> trunc_lo = TruncateHitsPCA(sortedHits[i_plane], hitslo, sortedHits[i_plane][i_hit], geo, dprop, fHitGroupDistance/2., fTruncateSigma);
          std::vector<art::Ptr<recob::Hit>> trunc_hi = TruncateHitsPCA(sortedHits[i_plane], hitshi, sortedHits[i_plane][i_hit], geo, dprop, fHitGroupDistance/2., fTruncateSigma);
          float truncated_pca_angle = HitPCAAngle(trunc_lo, trunc_hi, sortedHits[i_plane][i_hit], geo, dprop);

          // and save
          sbn::DThetaKink thisKink;
          thisKink.wire = sortedHits[i_plane][i_hit]->WireID();
          thisKink.pca = pca_angle;
          thisKink.trunc_pca = truncated_pca_angle;
          thisKink.vec = vec_angle;
          std::array<float, 2> hv = HitVector(*sortedHits[i_plane][i_hit], geo, dprop);
          thisKink.hit_wire_cm = hv[0];
          thisKink.hit_time_cm = hv[1];
          thisKink.complete = complete;
          thisKink.hitID = sortedHits[i_plane][i_hit].key();
          thisKink.hitOrder = i_hit;

          if (hitshi.size()) thisKink.hitIDHi = hitshi.back().key();
          else thisKink.hitIDHi = -1;
          if (hitslo.size()) thisKink.hitIDLo = hitslo.back().key();
          else thisKink.hitIDHi = -1;

          thisKink.branch = branch_id;

          // get the generation
          thisKink.generation = 0;
          art::Ptr<recob::PFParticle> check = id_to_pfp.at(branch_id);
          while (check->Self() != thisPFP->Self()) {
            thisKink.generation ++;
            check = id_to_pfp.at(check->Parent());
          }

          planeKinks[i_plane].push_back(thisKink);
        } 
      }
    }

    // select the kinks
    for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
      std::vector<sbn::DThetaKink> selected = SelectKinks(planeKinks[i_plane]);
      for (const sbn::DThetaKink &k: selected) {
        outKinks->push_back(k);
        art::Ptr<sbn::DThetaKink> thisKinkPtr = kinkPtrMaker(outKinks->size()-1);
        assn->addSingle(thisPFP, thisKinkPtr);
      }
    }

  }

  evt.put(std::move(outKinks));
  evt.put(std::move(assn));

}

DEFINE_ART_MODULE(sbn::DThetaKinkFinder)
