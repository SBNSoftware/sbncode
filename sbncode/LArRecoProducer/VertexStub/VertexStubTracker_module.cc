////////////////////////////////////////////////////////////////////////
// Class:       VertexStubTracker
// Plugin Type: producer (art v3_02_06)
// File:        VertexStubTracker_module.cc
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
#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcorealg/Geometry/Exceptions.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "sbncode/LArRecoProducer/Products/VertexHit.h"
#include "sbncode/LArRecoProducer/Products/Stub.h"

#include <memory>
#include <optional>

namespace sbn {
  class VertexStubTracker;
}


class sbn::VertexStubTracker : public art::EDProducer {
public:
  explicit VertexStubTracker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VertexStubTracker(VertexStubTracker const&) = delete;
  VertexStubTracker(VertexStubTracker&&) = delete;
  VertexStubTracker& operator=(VertexStubTracker const&) = delete;
  VertexStubTracker& operator=(VertexStubTracker&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // input labels
  art::InputTag fPFPLabel;
  art::InputTag fVertexChargeLabel;
  float fdQdxCut;
  calo::CalorimetryAlg fCaloAlg;

};

sbn::VertexStubTracker::VertexStubTracker(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFPLabel(p.get<art::InputTag>("PFPLabel", "pandora")),
    fVertexChargeLabel(p.get<art::InputTag>("VertexChargeLabel", "vhit")),
    fdQdxCut(p.get<float>("dQdxCut")),
    fCaloAlg(p.get<fhicl::ParameterSet >("CaloAlg"))
{
  produces<std::vector<sbn::Stub>>();
  produces<art::Assns<sbn::VertexHit, sbn::Stub>>();
  produces<art::Assns<sbn::Stub, recob::Hit>>();
  produces<art::Assns<sbn::Stub, recob::Slice>>();
  produces<art::Assns<sbn::Stub, recob::PFParticle>>();
  
}

art::Ptr<recob::PFParticle> FindPFP(
  const std::vector<art::Ptr<recob::Hit>> &hits, 
  const std::vector<art::Ptr<recob::PFParticle>> &pfps, 
  const std::vector<std::vector<art::Ptr<recob::Hit>>> &pfp_hits) {

  std::vector<unsigned> pfp_nhit(pfps.size(), 0);

  for (const art::Ptr<recob::Hit> h: hits) {
    for (unsigned i_pfp = 0; i_pfp < pfps.size(); i_pfp++) {
      const std::vector<art::Ptr<recob::Hit>> &thispfpHits = pfp_hits[i_pfp]; 
      bool found = false;
      for (const art::Ptr<recob::Hit> &pfp_h: thispfpHits) {
        if (h == pfp_h) {
          found = true;
          pfp_nhit[i_pfp] ++;
          break;
        }
      }
      if (found) break;
    }
  }

  for (unsigned i_pfp = 0; i_pfp < pfps.size(); i_pfp++) {
    if (pfp_nhit[i_pfp] >= hits.size() / 2) return pfps[i_pfp];
  }

  art::Ptr<recob::PFParticle> null;
  return null;
}

std::vector<art::Ptr<recob::Hit>> CollectHits(
    const std::vector<art::Ptr<recob::Hit>> &hits, 
    const sbn::VertexHit &vhit, 
    const recob::Hit &vhit_hit, 
    const geo::Point_t &vertex, 
    const geo::GeometryCore *geo, 
    const detinfo::DetectorPropertiesData &dprop) {

  // project the vertex onto the wireplane
  float vert_w = geo->WireCoordinate(vertex, vhit_hit.WireID());
  float vert_x = vertex.x();

  // get the vertex hit-coordinates
  float hit_w = vhit_hit.WireID().Wire;
  float hit_x = dprop.ConvertTicksToX(vhit_hit.PeakTime(), vhit_hit.WireID());
  geo::PlaneID vhit_plane = vhit_hit.WireID();

  // get the line-segment slope / intercept
  float slope = (hit_x - vert_x) / (hit_w - vert_w);
  float intercept = hit_x - slope * hit_w;

  // get all the hits that overlap between these two points
  std::vector<art::Ptr<recob::Hit>> ret;
  for (art::Ptr<recob::Hit> h: hits) {
    float this_hit_w = h->WireID().Wire;
    geo::PlaneID this_hit_plane = h->WireID();
    if (this_hit_plane == vhit_plane &&  // check plane
       ((this_hit_w <= hit_w && this_hit_w >= vert_w) || // check overlap
        (this_hit_w >= hit_w && this_hit_w <= vert_w))) {
      float this_proj_x = slope * this_hit_w + intercept;

      float h_time_lo = h->PeakTime() - h->RMS();
      float h_x_A = dprop.ConvertTicksToX(h_time_lo, h->WireID());
      float h_time_hi = h->PeakTime() + h->RMS();
      float h_x_B = dprop.ConvertTicksToX(h_time_hi, h->WireID());

      if ((h_x_A >= this_proj_x && h_x_B <= this_proj_x) ||
          (h_x_A <= this_proj_x && h_x_B >= this_proj_x)) {
        ret.push_back(h);
      }
    }
  }
  
  return ret;
}

void sbn::VertexStubTracker::produce(art::Event& e)
{
  // collect services
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  // output data products
  std::unique_ptr<std::vector<sbn::Stub>> outStubs(new std::vector<sbn::Stub>);
  std::unique_ptr<art::Assns<sbn::VertexHit, sbn::Stub>> assn(new art::Assns<sbn::VertexHit, sbn::Stub>);
  std::unique_ptr<art::Assns<sbn::Stub, recob::Hit>> hitAssn(new art::Assns<sbn::Stub, recob::Hit>);
  std::unique_ptr<art::Assns<sbn::Stub, recob::Slice>> slcAssn(new art::Assns<sbn::Stub, recob::Slice>);
  std::unique_ptr<art::Assns<sbn::Stub, recob::PFParticle>> pfpAssn(new art::Assns<sbn::Stub, recob::PFParticle>);

  art::PtrMaker<sbn::Stub> StubPtrMaker {e};

  // input data
  art::Handle<std::vector<sbn::VertexHit>> vhit_handle;
  e.getByLabel(fVertexChargeLabel, vhit_handle);

  std::vector<art::Ptr<sbn::VertexHit>> vhits;
  art::fill_ptr_vector(vhits, vhit_handle);

  art::FindManyP<recob::Slice> vhitSlcs(vhits, e, fVertexChargeLabel);
  art::FindManyP<recob::Vertex> vhitVtxs(vhits, e, fVertexChargeLabel);
  art::FindManyP<recob::Hit> vhitHits(vhits, e, fVertexChargeLabel);

  // Holders for saving data on slc info
  std::map<unsigned, std::vector<std::vector<art::Ptr<recob::Hit>>>> slicePFPHits;

  for (unsigned i_vhit = 0; i_vhit < vhits.size(); i_vhit++) {
    // Collect data on this Vertex-Hit
    const sbn::VertexHit &thisVHit = *vhits[i_vhit];
    const recob::Hit &thisVHitHit = *vhitHits.at(i_vhit).at(0);
    if (thisVHit.dqdx < fdQdxCut) continue;

    // Collect more data
    art::Ptr<recob::Slice> thisSlice = vhitSlcs.at(i_vhit).at(0);

    art::FindManyP<recob::PFParticle> findThisSlicePFP({thisSlice}, e, fPFPLabel);
    const std::vector<art::Ptr<recob::PFParticle>> &thisSlicePFPs = findThisSlicePFP.at(0);

    art::FindManyP<recob::Hit> findThisSliceHits({thisSlice}, e, fPFPLabel);
    const std::vector<art::Ptr<recob::Hit>> &hits = findThisSliceHits.at(0);

    // If we haven't seen this slice beore, collect the hits for its pfparticles
    if (!slicePFPHits.count(thisSlice.key())) {
      art::FindManyP<recob::Cluster> pfparticleClusters(thisSlicePFPs, e, fPFPLabel);

      std::vector<std::vector<art::Ptr<recob::Hit>>> thisSlicePFPHits;
      for (unsigned i_slc_pfp = 0; i_slc_pfp < thisSlicePFPs.size(); i_slc_pfp++) {
        thisSlicePFPHits.emplace_back();
        const std::vector<art::Ptr<recob::Cluster>> &thisPFPClusters = pfparticleClusters.at(i_slc_pfp);
        art::FindManyP<recob::Hit> clusterHits(thisPFPClusters, e, fPFPLabel);
        for (unsigned i_clus = 0; i_clus < thisPFPClusters.size(); i_clus++) {
          thisSlicePFPHits.back().insert(thisSlicePFPHits.back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
        }
      }
      slicePFPHits[thisSlice.key()] = thisSlicePFPHits;
    }

    const std::vector<std::vector<art::Ptr<recob::Hit>>> &thisSlicePFPHits = slicePFPHits.at(thisSlice.key());

    const recob::Vertex &vert = *vhitVtxs.at(i_vhit).at(0);
    TVector3 vert_v(vert.position().X(), vert.position().Y(), vert.position().Z());

    std::vector<art::Ptr<recob::Hit>> thisHits = CollectHits(hits, thisVHit, thisVHitHit, vert.position(), geo, dprop);

    // Find which pfparticle to group this VertexStub with 
    art::Ptr<recob::PFParticle> thisStubPFP = FindPFP(thisHits, thisSlicePFPs, thisSlicePFPHits);

    float charge = 0.;
    std::set<int> wires;
    for (unsigned i_hit = 0; i_hit < thisHits.size(); i_hit++) {
      const recob::Hit &hit = *thisHits[i_hit];
      charge += fCaloAlg.ElectronsFromADCArea(hit.Integral(), hit.WireID().Plane) * fCaloAlg.LifetimeCorrection(clock_data, dprop, hit.PeakTime(), 0.);
      wires.insert(hit.WireID().Wire);
    }

    float vert_w = geo->WireCoordinate(vert.position(), thisVHitHit.WireID());
    float hit_w = thisVHitHit.WireID().Wire;
    float wire_dist = abs(hit_w - vert_w);

    sbn::Stub stub;
    stub.charge = charge;
    stub.vtx = vert_v;
    stub.end = thisVHit.spXYZ;
    stub.pitch = thisVHit.charge / thisVHit.dqdx;
    stub.nwire = wires.size();
    stub.wire_dist = wire_dist;
    stub.plane = thisVHit.wire;

    // Save!
    outStubs->push_back(stub);
    art::Ptr<sbn::Stub> outStub = StubPtrMaker(outStubs->size() - 1);
    assn->addSingle(vhits[i_vhit], outStub);
    for (unsigned i_hit = 0; i_hit < thisHits.size(); i_hit++) { 
      hitAssn->addSingle(outStub, thisHits[i_hit]);
    }
    slcAssn->addSingle(outStub, thisSlice);
    if (thisStubPFP) {
      pfpAssn->addSingle(outStub, thisStubPFP);
    }

  } // end iterate over vertex hits

  // Save into event
  e.put(std::move(outStubs));
  e.put(std::move(assn));
  e.put(std::move(hitAssn));
  e.put(std::move(slcAssn));
  e.put(std::move(pfpAssn));
}

DEFINE_ART_MODULE(sbn::VertexStubTracker)
