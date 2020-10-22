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

#include "Products/VertexHit.h"
#include "PCA.h"

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
  art::InputTag fHitLabel;
  float fHitVacuumRadius;

  // private data 

};


sbn::VertexChargeVacuum::VertexChargeVacuum(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabel(p.get<std::string>("PFParticleLabel", "pandora")),
    fHitLabel(p.get<std::string>("HitLabel", "gaushit")),
    fHitVacuumRadius(p.get<float>("HitVacuumRadius"))
{
}

void sbn::VertexChargeVacuum::produce(art::Event& evt)
{
  // output stuff
  std::unique_ptr<std::vector<sbn::VertexHit>> outVHit(new std::vector<sbn::VertexHit>);
  std::unique_ptr<art::Assns<recob::PFParticle, sbn::VertexHit>> assn(new art::Assns<recob::PFParticle, sbn::VertexHit>);

  art::PtrMaker<sbn::VertexHit> vhitPtrMaker {evt};

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

  art::FindManyP<recob::Vertex> pfparticleVertices(pfparticles, evt, fPFParticleLabel);
  art::FindManyP<recob::Cluster> pfparticleClusters(pfparticles, evt, fPFParticleLabel);

  // organize the PFPlist into a map
  std::map<unsigned, art::Ptr<recob::PFParticle>> id_to_pfp;
  for (unsigned i = 0; i < pfparticles.size(); i++) {
    id_to_pfp[pfparticles[i]->Self()] = pfparticles[i];
  }

  // get all the hits
  art::Handle<std::vector<recob::Hit>> hit_handle;
  evt.getByLabel(fHitLabel, hit_handle);

  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hit_handle);

  // iterate over the "primary" PFParticles
  for (unsigned i_pfp = 0; i_pfp < pfparticles.size(); i_pfp++) {
    const recob::PFParticle &pfp = *pfparticles[i_pfp];
    if (!pfp.IsPrimary()) continue;

    // we found a primary PFP! Get its vertex.
    const recob::Vertex &vert = *pfparticleVertices.at(i_pfp).at(0);

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

    // work on each plane
    for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
      std::vector<std::array<float, 2>> daughterPlanePCs;
      std::vector<int> daughterPFPIDs;
      // Get the PC-axis of each PFParticle daughter projected on this plane
      for (unsigned i_d = 0; i_d < daughterPFPs.size(); i_d++) {
        std::array<float, 2> thisPC = sbnpca::HitPCAVec(daughterPlaneHits[i_plane][i_d], vert, geo, dprop);
        if (thisPC[0] > -99. && thisPC[1] > -99.) {
          daughterPlanePCs.push_back(thisPC);
          daughterPFPIDs.push_back(daughterPFPs[i_d]->Self());
        }
       
      }

      // vacuum up all the hits within the radius
      std::vector<art::Ptr<recob::Hit>> nearbyHits;
      for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
        const recob::Hit &hit = *hits[i_hit];
        if (hit.WireID().Plane == i_plane) {
          if (sbnpca::Vert2HitDistance(hit, vert, geo, dprop) < fHitVacuumRadius) {
            nearbyHits.push_back(hits[i_hit]);
          }
        }
      }

      // now for each hit compute the dot with each PFP PC
      for (unsigned i_hit = 0; i_hit < nearbyHits.size(); i_hit++) {
        const recob::Hit &hit = *nearbyHits[i_hit];

	sbn::VertexHit vhit;
        vhit.wire = hit.WireID();
	vhit.charge = hit.Integral();
	vhit.proj_dist_to_vertex = sbnpca::Vert2HitDistance(hit, vert, geo, dprop);
	vhit.nearbyPFPIDs = daughterPFPIDs;

        for (unsigned i_pfp = 0; i_pfp < daughterPFPIDs.size(); i_pfp++) {
          std::array<float, 2> pc = daughterPlanePCs[i_pfp];
          std::array<float, 2> hv = sbnpca::HitVector(hit, geo, dprop);
          float dot = (pc[0]*hv[0] + pc[1]*hv[1]) / 
            (sqrt(hv[0]*hv[0]+hv[1]+hv[1]) *
             sqrt(pc[0]*pc[0]+pc[1]+pc[1]));
          vhit.PCproj.push_back(dot);
        }

        // Save!
        outVHit->push_back(vhit);
        art::Ptr<sbn::VertexHit> thisVHitPtr = vhitPtrMaker(outVHit->size()-1);
        assn->addSingle(pfparticles[i_pfp], thisVHitPtr);
      }
    }
  } 

  evt.put(std::move(outVHit));
  evt.put(std::move(assn));

}

DEFINE_ART_MODULE(sbn::VertexChargeVacuum)
