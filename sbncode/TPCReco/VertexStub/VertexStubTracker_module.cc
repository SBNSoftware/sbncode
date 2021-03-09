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
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"
#include "sbncode/TPCReco/VertexStub/StubBuilder.h"

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
  sbn::StubBuilder fStubBuilder;

};

sbn::VertexStubTracker::VertexStubTracker(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFPLabel(p.get<art::InputTag>("PFPLabel", "pandora")),
    fVertexChargeLabel(p.get<art::InputTag>("VertexChargeLabel", "vhit")),
    fdQdxCut(p.get<float>("dQdxCut")),
    fStubBuilder(p.get<fhicl::ParameterSet >("CaloAlg"))
{
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

  // Setup the stub builder
  fStubBuilder.Setup(e, fPFPLabel);

  // Holders for saving data on slc info
  std::map<unsigned, std::vector<std::vector<art::Ptr<recob::Hit>>>> slicePFPHits;

  for (unsigned i_vhit = 0; i_vhit < vhits.size(); i_vhit++) {
    // Collect data on this Vertex-Hit
    const sbn::VertexHit &thisVHit = *vhits[i_vhit];
    const recob::Hit &thisVHitHit = *vhitHits.at(i_vhit).at(0);
    if (thisVHit.dqdx < fdQdxCut) continue;

    // Collect more data
    art::Ptr<recob::Slice> thisSlice = vhitSlcs.at(i_vhit).at(0);
    const recob::Vertex &thisVertex = *vhitVtxs.at(i_vhit).at(0);

    std::vector<art::Ptr<recob::Hit>> thisHits;
    art::Ptr<recob::PFParticle> thisStubPFP;
    sbn::Stub stub = fStubBuilder.FromVertexHit(thisSlice, thisVHit, thisVHitHit, thisVertex, geo, clock_data, dprop, thisHits, thisStubPFP); 

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
