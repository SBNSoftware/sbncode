#ifndef StubBuilder_HH
#define StubBuilder_HH

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"

#include <map>
#include <memory>

namespace sbn {

class StubBuilder {
public:
  void Setup(const art::Event &e, const art::InputTag &pfplabel, const art::InputTag &trklabel);
  // Build a 1-Plane stub 
  sbn::Stub FromVertexHit(const art::Ptr<recob::Slice> &slice,
	    const sbn::VertexHit &vhit,
	    const recob::Hit &vhit_hit, 
	    const recob::Vertex &vertex,
	    const geo::GeometryCore *geo,
            const spacecharge::SpaceCharge *sce,
	    const detinfo::DetectorClocksData &dclock,
	    const detinfo::DetectorPropertiesData &dprop,
	    std::vector<art::Ptr<recob::Hit>> &stub_hits,
	    art::Ptr<recob::PFParticle> &stub_pfp);

  StubBuilder(fhicl::ParameterSet const& p): fCaloAlg(p) {}

private:
  calo::CalorimetryAlg fCaloAlg;
  std::map<unsigned, std::vector<art::Ptr<recob::Hit>>> fSliceHits;
  std::map<unsigned, std::vector<art::Ptr<recob::PFParticle>>> fSlicePFPs;
  std::map<unsigned, std::vector<art::Ptr<recob::Track>>> fSliceTrks;
  std::map<unsigned, std::vector<std::vector<art::Ptr<recob::Hit>>>> fSlicePFPHits;
  std::map<unsigned, std::vector<std::vector<art::Ptr<recob::Hit>>>> fSliceTrkHits;
  std::map<unsigned, std::vector<std::vector<const recob::TrackHitMeta *>>> fSliceTrkTHMs;

};

} // end namespace sbn
#endif
