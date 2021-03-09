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

#include <memory>

namespace sbn {

class StubBuilder {
public:
  void Setup(const art::Event &e, const art::InputTag &pfplabel);
  sbn::Stub FromVertexHit(const art::Ptr<recob::Slice> &slice,
	    const sbn::VertexHit &vhit,
	    const recob::Hit &vhit_hit, 
	    const recob::Vertex &vertex,
	    const geo::GeometryCore *geo,
	    const detinfo::DetectorClocksData &dclock,
	    const detinfo::DetectorPropertiesData &dprop,
	    std::vector<art::Ptr<recob::Hit>> &stub_hits,
	    art::Ptr<recob::PFParticle> &stub_pfp);

  StubBuilder(fhicl::ParameterSet const& p): fCaloAlg(p) {}

private:
  calo::CalorimetryAlg fCaloAlg;
  std::map<unsigned, std::vector<art::Ptr<recob::Hit>>> fSliceHits;
  std::map<unsigned, std::vector<art::Ptr<recob::PFParticle>>> fSlicePFPs;
  std::map<unsigned, std::vector<std::vector<art::Ptr<recob::Hit>>>> fSlicePFPHits;

};

} // end namespace sbn
#endif
