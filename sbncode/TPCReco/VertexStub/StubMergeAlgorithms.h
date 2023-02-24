#ifndef StubMergeAlgorithms_HH
#define StubMergeAlgorithms_HH

#include "canvas/Persistency/Common/Ptr.h"

#include "larcorealg/Geometry/fwd.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

#include "larevt/SpaceCharge/SpaceCharge.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"

// Helper algorithms/utilities useful for merging stubs

namespace sbn {

/// Internal struct: contains information on stub and other associated data products
struct StubInfo {
  sbn::Stub stub;
  art::Ptr<recob::PFParticle> pfp;
  std::vector<art::Ptr<recob::Hit>> hits;
  art::Ptr<sbn::VertexHit> vhit;
  art::Ptr<recob::Hit> vhit_hit;
};

/**
 * @brief Computes the track-pitch on a plane given an input direction and location.
 * @param loc location in world coordinate system [cm]
 * @param dir direction in world coordinate system (assumed to be unit vector)
 * @param view the view to project the direction onto
 * @param tpc ID of the TPC the selected `view` belongs to
 * @param correct_sce whether to apply space charge correction
 * @param track_is_sce_corrected whether location `loc` is already corrected for 
 *        spacecharge effects
 * @param xsign factor for the _x_ (drift) coordinate (`-1` flips it)
 * @return the pitch along direction `dir` [cm]
 * 
 * Is able to handle the presence of space charge (through the `SpaceCharge` service)
 * and whether or not the input location/direction are already corrected for SCE.
 */
double GetPitch(
    const geo::GeometryCore &geom,
    const geo::WireReadoutGeom &channelMap, const spacecharge::SpaceCharge *sce,
    geo::Point_t loc, geo::Vector_t dir, 
    geo::View_t view, geo::TPCID tpc, 
    bool correct_sce, bool track_is_sce_corrected, float xsign=1.);

/// Get the location in the presence of space charge
geo::Point_t GetLocation(const spacecharge::SpaceCharge *sce, geo::Point_t loc_w, geo::TPCID TPC, float xsign=1.);

/// Get the E-Field in the presence of space charge
double GetEfield(const detinfo::DetectorPropertiesData &dprop, const spacecharge::SpaceCharge *sce, geo::Point_t loc, geo::TPCID TPC, bool correct_loc_sce, float xsign=1.);

/// Get the SCE-distorted location (i.e. the location "seen" by the wireplanes)
geo::Point_t GetLocationAtWires(const spacecharge::SpaceCharge *sce, geo::Point_t loc, geo::Vector_t dir, float xsign=1.);

/// Returns whether stub `A` contains stub `B`.
bool StubContains(const sbn::StubInfo &A, const sbn::StubInfo &B);

/// Computes the dot product of two stubs
float StubDirectionDot(const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const geo::WireReadoutGeom &channelMap,
    const detinfo::DetectorPropertiesData &dprop);

// For matching stubs across planes
float StubTimeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const detinfo::DetectorClocksData &dclock,
    const detinfo::DetectorPropertiesData &dprop);

/// Returns an updated end position of a stub after merging across two planes
geo::Point_t TwoStubEndPosition(const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::WireReadoutGeom &channelMap,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop);

/// Difference of the total charge between two stubs
float StubChargeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B);

/// Difference of the endpoint charge between two stubs
float StubPeakChargeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B);

/// Difference of the endpoint dQ/dx between two stubs
float StubPeakdQdxOffset(const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::GeometryCore &geom,
    const geo::WireReadoutGeom &channelMap,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop);

} // end namespace sbn
#endif
