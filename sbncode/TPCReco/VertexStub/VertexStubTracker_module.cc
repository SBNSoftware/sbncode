////////////////////////////////////////////////////////////////////////
// Class:       VertexStubTracker
// Plugin Type: producer (art v3_02_06)
// File:        VertexStubTracker_module.cc
// Author:      grayputnam@uchicago.edu
//
// Art module for building "Stub" reconstruction objects. Stubs are 
// correspond to low energy hadrons produced in neutrino interactions that
// produce short, highly-ionizing charge depositions near the vertex. 
//
// Module operates on sbn::VertexHit objects and connects the end-point Hit
// to the start-point Vertex to build a stub. Merging is then done between
// and across planes to consoldiate and refine output. Produces the sbn::Stub
// object.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
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
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"
#include "sbncode/TPCReco/VertexStub/StubBuilder.h"
#include "sbncode/TPCReco/VertexStub/IStubMerge.h"

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
  art::InputTag fTrackLabel;
  art::InputTag fVertexChargeLabel;
  float fdQdxCut;
  float fOneWiredQdxCut;
  sbn::StubBuilder fStubBuilder;
  std::vector<std::unique_ptr<sbn::IStubMerge>> fStubMergeTools;
};

sbn::VertexStubTracker::VertexStubTracker(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFPLabel(p.get<art::InputTag>("PFPLabel", "pandora")),
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrack")),
    fVertexChargeLabel(p.get<art::InputTag>("VertexChargeLabel", "vhit")),
    fdQdxCut(p.get<float>("dQdxCut")),
    fOneWiredQdxCut(p.get<float>("OneWiredQdxCut")),
    fStubBuilder(p.get<fhicl::ParameterSet >("CaloAlg"))
{
  // load the tools
  std::vector<fhicl::ParameterSet> merge_tool_configs(p.get<std::vector<fhicl::ParameterSet>>("MergeTools"));
  for (unsigned i = 0; i < merge_tool_configs.size(); i++) {
    fStubMergeTools.push_back(art::make_tool<IStubMerge>(merge_tool_configs[i]));
  }

  produces<std::vector<sbn::Stub>>();
  produces<art::Assns<sbn::VertexHit, sbn::Stub>>();
  produces<art::Assns<sbn::Stub, recob::Hit>>();
  produces<art::Assns<sbn::Stub, recob::Slice>>();
  produces<art::Assns<sbn::Stub, recob::PFParticle>>();
  
}

void sbn::VertexStubTracker::produce(art::Event& e)
{
  // collect services
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
  // TODO: fix -- for now, use a null space-charge service
  const spacecharge::SpaceCharge *sce = nullptr;

  // output data products
  std::unique_ptr<std::vector<sbn::Stub>> outStubs(new std::vector<sbn::Stub>);
  std::unique_ptr<art::Assns<sbn::VertexHit, sbn::Stub>> assn(new art::Assns<sbn::VertexHit, sbn::Stub>);
  std::unique_ptr<art::Assns<sbn::Stub, recob::Hit>> hitAssn(new art::Assns<sbn::Stub, recob::Hit>);
  std::unique_ptr<art::Assns<sbn::Stub, recob::Slice>> slcAssn(new art::Assns<sbn::Stub, recob::Slice>);
  std::unique_ptr<art::Assns<sbn::Stub, recob::PFParticle>> pfpAssn(new art::Assns<sbn::Stub, recob::PFParticle>);

  art::PtrMaker<sbn::Stub> StubPtrMaker {e};

  // input data
  art::Handle<std::vector<recob::Slice>> slice_handle;
  e.getByLabel(fPFPLabel, slice_handle);

  std::vector<art::Ptr<recob::Slice>> slices;
  art::fill_ptr_vector(slices, slice_handle);

  art::FindManyP<sbn::VertexHit> slcVHits(slices, e, fVertexChargeLabel);

  // Setup the stub builder
  fStubBuilder.Setup(e, fPFPLabel, fTrackLabel);

  for (unsigned i_slc = 0; i_slc < slices.size(); i_slc++) {
    const std::vector<art::Ptr<sbn::VertexHit>> &vhits = slcVHits.at(i_slc);
    
    // look up info per vertex hit
    art::FindManyP<recob::Vertex> vhitVtxs(vhits, e, fVertexChargeLabel);
    art::FindManyP<recob::Hit> vhitHits(vhits, e, fVertexChargeLabel);

    // stuff common to each vertex hit in a slice
    art::Ptr<recob::Slice> thisSlice = slices[i_slc];
    // each hit should have the same vertex
    recob::Vertex vertex;
    if (!vhits.empty()) vertex = *vhitVtxs.at(0).at(0);
  
    std::vector<sbn::StubInfo> stubs;
    for (unsigned i_vhit = 0; i_vhit < vhits.size(); i_vhit++) {
      // Collect data on this Vertex-Hit
      const sbn::VertexHit &thisVHit = *vhits[i_vhit];
      const recob::Hit &thisVHitHit = *vhitHits.at(i_vhit).at(0);

      bool passcut = (thisVHit.dqdx >= fdQdxCut) || ((abs(thisVHit.wire.Wire - thisVHit.vtxw) < 1.) && (thisVHit.dqdx >= fOneWiredQdxCut));

      if (!passcut) continue;

      sbn::StubInfo sinfo;
      sinfo.stub = fStubBuilder.FromVertexHit(thisSlice, thisVHit, thisVHitHit, vertex, geo, sce, clock_data, dprop, sinfo.hits, sinfo.pfp); 
      sinfo.vhit = vhits[i_vhit];
      sinfo.vhit_hit = vhitHits.at(i_vhit).at(0);

      stubs.push_back(sinfo);

    } // end iterate over vertex hits

    // Run all of the merging tools
    for (unsigned i_mrg = 0; i_mrg < fStubMergeTools.size(); i_mrg++) {
      stubs = fStubMergeTools[i_mrg]->Merge(stubs, geo, sce, clock_data, dprop);
    }

    // Save!
    for (unsigned i_stub = 0; i_stub < stubs.size(); i_stub++) {
      const sbn::StubInfo &sinfo = stubs[i_stub];

      outStubs->push_back(sinfo.stub);
      art::Ptr<sbn::Stub> outStub = StubPtrMaker(outStubs->size() - 1);
      assn->addSingle(sinfo.vhit, outStub);
      for (unsigned i_hit = 0; i_hit < sinfo.hits.size(); i_hit++) { 
        hitAssn->addSingle(outStub, sinfo.hits[i_hit]);
      }
      slcAssn->addSingle(outStub, thisSlice);
      if (sinfo.pfp) {
        pfpAssn->addSingle(outStub, sinfo.pfp);
      }
    } // end loop over stubs

  } // end loop over slices
  // Save into event
  e.put(std::move(outStubs));
  e.put(std::move(assn));
  e.put(std::move(hitAssn));
  e.put(std::move(slcAssn));
  e.put(std::move(pfpAssn));
}

DEFINE_ART_MODULE(sbn::VertexStubTracker)
