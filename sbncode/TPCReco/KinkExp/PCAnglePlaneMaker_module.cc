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

#include "sbncode/TPCReco/Products/PCAnglePlane.h"
#include "PCA.h"

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

std::array<std::vector<art::Ptr<recob::Hit>>, 3> SortHits(const std::array<std::vector<art::Ptr<recob::Hit>>, 3> &hits, 
                                                          const recob::Vertex &start,
                                                          const geo::GeometryCore *geo, 
                                                          const detinfo::DetectorPropertiesData &dprop) {
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
            return sbnpca::Vert2HitDistance(*lhs, start, geo, dprop) < sbnpca::Vert2HitDistance(*rhs, start, geo, dprop);
      });
      ret[i].push_back(*firstHit);
      planar[i].erase(firstHit);

      while (planar[i].size()) {
        const recob::Hit &lastHit = *ret[i].back();
        std::list<art::Ptr<recob::Hit>>::iterator closest = std::min_element(planar[i].begin(), planar[i].end(),
          [lastHit, geo, dprop](auto const &lhs, auto const &rhs) { 
            return sbnpca::HitDistance(*lhs, lastHit, geo, dprop) < sbnpca::HitDistance(*rhs, lastHit, geo, dprop);
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
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);


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
          auto [hitslo, hitshi, complete] = sbnpca::GetNearestHits(sortedHits[i_plane], i_hit, fHitGroupDistance, geo, dprop);

          //std::cout << "Plane: " << i_plane << " Branch: " << branch_id << " hit: " << i_hit << " " << sortedHits[i_plane][i_hit].key();
          //if (hitslo.size()) std::cout << " lo: " << hitslo.back().key();
          //else std::cout << " lo: none";
          //if (hitshi.size()) std::cout << " hi: " << hitshi.back().key();
          //else std::cout << " hi: none";
          //std::cout << std::endl;

          // invalid -- continue
          // if (hithi.size() < 2 || hitslo.size() < 2) continue;

          // Get the PCA dirs
          std::array<float, 2> pca_vec_lo = sbnpca::HitPCAVec(hitslo, *sortedHits[i_plane][i_hit], geo, dprop);
          std::array<float, 2> pca_vec_hi = sbnpca::HitPCAVec(hitshi, *sortedHits[i_plane][i_hit], geo, dprop);
          float angle = -100.;
          // check if valid vecs
          if (pca_vec_lo[0] > -99. && pca_vec_lo[1] > -99. && pca_vec_hi[0] > -99. && pca_vec_hi[1] > -99.) {
            // get angle between the two vecs and flip
            angle = M_PI - sbnpca::VecAngle(pca_vec_lo, pca_vec_hi);
          }

          std::array<float, 2> hv = sbnpca::HitVector(*sortedHits[i_plane][i_hit], geo, dprop);

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

          if (hitshi.size()) thisAngle.angle.dist_to_hi = sbnpca::HitDistance(*hitshi.back(), *sortedHits[i_plane][i_hit], geo, dprop); 
          else thisAngle.angle.dist_to_hi = -100.;
          if (hitslo.size()) thisAngle.angle.dist_to_lo = sbnpca::HitDistance(*hitslo.back(), *sortedHits[i_plane][i_hit], geo, dprop); 
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
