#include "StubBuilder.h"
#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

// Helper functions
bool HitOnTrack(const art::Ptr<recob::Hit> &hit,
    const art::Ptr<recob::Track> &trk,
    const std::vector<art::Ptr<recob::Hit>> &trk_hits,
    const std::vector<const recob::TrackHitMeta *> &trk_thms) {

  if (!trk) return false; // no track: hit can't be on a track
    
  for (unsigned i = 0; i < trk_hits.size(); i++) {
    // Found the hit on the track!
    if (hit == trk_hits[i]) { 
      // Is the hit part of the track Calo?
      if (trk_thms[i]->Index() != std::numeric_limits<int>::max() && trk->HasValidPoint(trk_thms[i]->Index())) { 
        return true;
      }
      break;
    }
  }
  return false;
}

int FindPFPInd(
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
    if (pfp_nhit[i_pfp] >= hits.size() / 2) return i_pfp;
  }

  return -1;
}

std::vector<art::Ptr<recob::Hit>> CollectHits(
    const std::vector<art::Ptr<recob::Hit>> &hits, 
    const sbn::VertexHit &vhit, 
    const recob::Hit &vhit_hit, 
    const geo::Point_t &vertex, 
    const geo::GeometryCore *geo, 
    const detinfo::DetectorPropertiesData &dprop) {

  // project the vertex onto the wireplane
  float vert_wf = geo->WireCoordinate(vertex, vhit_hit.WireID());
  float vert_x = vertex.x();

  // get the vertex hit-coordinates
  int hit_w = vhit_hit.WireID().Wire;
  float hit_x = dprop.ConvertTicksToX(vhit_hit.PeakTime(), vhit_hit.WireID());
  geo::PlaneID vhit_plane = vhit_hit.WireID();

  // get the line-segment slope / intercept
  float slope = (hit_x - vert_x) / (hit_w - vert_wf);
  float intercept = hit_x - slope * hit_w;

  // Include wires one after the hit and one before the vtx
  int vert_w = (int)((hit_w > vert_wf) ? std::floor(vert_wf) : std::ceil(vert_wf));
  hit_w = hit_w + ((hit_w > vert_wf) ? 1 : -1);

  // get all the hits that overlap between these two points
  std::vector<art::Ptr<recob::Hit>> ret;
  for (art::Ptr<recob::Hit> h: hits) {
    int this_hit_w = h->WireID().Wire;
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
  

void sbn::StubBuilder::Setup(const art::Event &e, const art::InputTag &pfplabel, const art::InputTag &trklabel) {
  // Clear out old data
  fSlicePFPHits.clear();
  fSlicePFPs.clear();
  fSliceHits.clear();
  fSliceTrkHits.clear();
  fSliceTrkTHMs.clear();
  fSliceTrks.clear();

  // Build the map of the slice key to the hits
  art::Handle<std::vector<recob::Slice>> slice_handle;
  e.getByLabel(pfplabel, slice_handle);

  std::vector<art::Ptr<recob::Slice>> slices;
  art::fill_ptr_vector(slices, slice_handle);

  art::FindManyP<recob::PFParticle> slicePFPs(slices, e, pfplabel);
  art::FindManyP<recob::Hit> sliceHits(slices, e, pfplabel);

  for (unsigned i_slc = 0; i_slc < slices.size(); i_slc++) {
    const std::vector<art::Ptr<recob::PFParticle>> thisSlicePFPs = slicePFPs.at(i_slc);
    art::FindManyP<recob::Cluster> pfparticleClusters(thisSlicePFPs, e, pfplabel);

    std::vector<art::Ptr<recob::Track>> thisSliceTracks;
    std::vector<std::vector<art::Ptr<recob::Hit>>> thisSlicePFPHits;
    std::vector<std::vector<art::Ptr<recob::Hit>>> thisSliceTrkHits;
    std::vector<std::vector<const recob::TrackHitMeta *>> thisSliceTrkTHMs;

    art::FindManyP<recob::Track> pfparticleTracks(thisSlicePFPs, e, trklabel);
    for (unsigned i_slc_pfp = 0; i_slc_pfp < thisSlicePFPs.size(); i_slc_pfp++) {
      thisSlicePFPHits.emplace_back();
      const std::vector<art::Ptr<recob::Cluster>> &thisPFPClusters = pfparticleClusters.at(i_slc_pfp);
      art::FindManyP<recob::Hit> clusterHits(thisPFPClusters, e, pfplabel);
      for (unsigned i_clus = 0; i_clus < thisPFPClusters.size(); i_clus++) {
        thisSlicePFPHits.back().insert(thisSlicePFPHits.back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
      }

      art::FindManyP<recob::Hit, recob::TrackHitMeta> FMthisTrackHits(pfparticleTracks.at(i_slc_pfp), e, trklabel);
      if (pfparticleTracks.at(i_slc_pfp).size()) {
        const std::vector<art::Ptr<recob::Hit>> &thisTrackHits = FMthisTrackHits.at(0);
        const std::vector<const recob::TrackHitMeta *> &thisTrackTHMs = FMthisTrackHits.data(0);
        thisSliceTrkHits.push_back(thisTrackHits);
        thisSliceTrkTHMs.push_back(thisTrackTHMs);
        thisSliceTracks.push_back(pfparticleTracks.at(i_slc_pfp).at(0));
      }
      else {
        thisSliceTrkHits.emplace_back();
        thisSliceTrkTHMs.emplace_back();
        thisSliceTracks.emplace_back(); // nullptr
      }

    }

    std::size_t const sliceKey = slices[i_slc].key();

    fSlicePFPHits[sliceKey] = thisSlicePFPHits;
    fSliceTrkHits[sliceKey] = thisSliceTrkHits;
    fSliceTrkTHMs[sliceKey] = thisSliceTrkTHMs;
    fSliceTrks[sliceKey] = thisSliceTracks;
    fSlicePFPs[sliceKey] = thisSlicePFPs;
    fSliceHits[sliceKey] = sliceHits.at(i_slc);
  }
}

sbn::Stub sbn::StubBuilder::FromVertexHit(const art::Ptr<recob::Slice> &slice,
                                  const sbn::VertexHit &vhit,
                                  const recob::Hit &vhit_hit, 
                                  const recob::Vertex &vertex,
                                  const geo::GeometryCore *geo,
                                  const spacecharge::SpaceCharge *sce,
                                  const detinfo::DetectorClocksData &dclock,
                                  const detinfo::DetectorPropertiesData &dprop,
                                  std::vector<art::Ptr<recob::Hit>> &stub_hits,
                                  art::Ptr<recob::PFParticle> &stub_pfp) {
    // look up stuff
    const std::vector<art::Ptr<recob::Hit>> &hits = fSliceHits.at(slice.key());
    const std::vector<art::Ptr<recob::PFParticle>> &pfps = fSlicePFPs.at(slice.key());
    const std::vector<art::Ptr<recob::Track>> &trks = fSliceTrks.at(slice.key());
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &pfp_hits = fSlicePFPHits.at(slice.key());
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &trk_hits = fSliceTrkHits.at(slice.key());
    const std::vector<std::vector<const recob::TrackHitMeta *>> &trk_thms = fSliceTrkTHMs.at(slice.key());
    
    TVector3 vertex_v(vertex.position().X(), vertex.position().Y(), vertex.position().Z());

    stub_hits = CollectHits(hits, vhit, vhit_hit, vertex.position(), geo, dprop);

    // Sort hits along the vtx -> end direction
    float vertex_w = geo->WireCoordinate(vertex.position(), vhit_hit.WireID());
    int stubdir = vertex_w <= vhit_hit.WireID().Wire ? 1 : -1;
    std::sort(stub_hits.begin(), stub_hits.end(), 
      [stubdir](auto const &lhs, auto const &rhs) {
        return lhs->WireID().Wire * stubdir < rhs->WireID().Wire * stubdir;});

    // Find which pfparticle to group this VertexStub with 
    int pfp_ind = FindPFPInd(stub_hits, pfps, pfp_hits);
    art::Ptr<recob::PFParticle> null; // nullptr
    stub_pfp = (pfp_ind >= 0) ? pfps[pfp_ind] : null;

    sbn::Stub stub;

    // Save all the charges
    stub.hits.emplace_back();
    for (unsigned i_hit = 0; i_hit < stub_hits.size(); i_hit++) {
      const recob::Hit &hit = *stub_hits[i_hit];

      sbn::StubHit stubhit; 

      stubhit.charge = fCaloAlg.ElectronsFromADCArea(hit.Integral(), hit.WireID().Plane) * fCaloAlg.LifetimeCorrection(dclock, dprop, hit.PeakTime(), 0.);
      stubhit.ontrack = (pfp_ind >= 0) ? HitOnTrack(stub_hits[i_hit], trks[pfp_ind], trk_hits[pfp_ind], trk_thms[pfp_ind]) : false;
      stubhit.wire = hit.WireID().Wire;

      stub.hits.back().push_back(stubhit);
    }

    // See if we can compute a track pitch
    if (pfp_ind >= 0 && trks[pfp_ind]) {
      stub.trkpitch.push_back(sbn::GetPitch(geo, sce, trks[pfp_ind]->Start(), trks[pfp_ind]->StartDirection(), geo->View(vhit_hit.WireID()), vhit_hit.WireID(), true, true));
    }
    else {
      stub.trkpitch.push_back(-1.);
    }

    stub.vtx = vertex_v;
    stub.end = vhit.spXYZ;
    stub.pitch.push_back(vhit.charge / vhit.dqdx);
    stub.plane.push_back(vhit.wire);
    stub.hit_w.push_back(vhit.wire.Wire);
    stub.vtx_w.push_back(vertex_w);

    return stub;
}

