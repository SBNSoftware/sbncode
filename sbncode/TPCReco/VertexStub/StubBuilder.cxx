#include "StubBuilder.h"

// Helper functions
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
  

void sbn::StubBuilder::Setup(const art::Event &e, const art::InputTag &pfplabel) {
  // Clear out old data
  fSlicePFPHits.clear();
  fSlicePFPs.clear();
  fSliceHits.clear();

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
    std::vector<std::vector<art::Ptr<recob::Hit>>> thisSlicePFPHits;
    for (unsigned i_slc_pfp = 0; i_slc_pfp < thisSlicePFPs.size(); i_slc_pfp++) {
      thisSlicePFPHits.emplace_back();
      const std::vector<art::Ptr<recob::Cluster>> &thisPFPClusters = pfparticleClusters.at(i_slc_pfp);
      art::FindManyP<recob::Hit> clusterHits(thisPFPClusters, e, pfplabel);
      for (unsigned i_clus = 0; i_clus < thisPFPClusters.size(); i_clus++) {
        thisSlicePFPHits.back().insert(thisSlicePFPHits.back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
      }
    }
    fSlicePFPHits[slices[i_slc].key()] = thisSlicePFPHits;
    fSlicePFPs[slices[i_slc].key()] = thisSlicePFPs;
    fSliceHits[slices[i_slc].key()] = sliceHits.at(i_slc);
  }
}

sbn::Stub sbn::StubBuilder::FromVertexHit(const art::Ptr<recob::Slice> &slice,
                                  const sbn::VertexHit &vhit,
                                  const recob::Hit &vhit_hit, 
                                  const recob::Vertex &vertex,
                                  const geo::GeometryCore *geo,
                                  const detinfo::DetectorClocksData &dclock,
                                  const detinfo::DetectorPropertiesData &dprop,
                                  std::vector<art::Ptr<recob::Hit>> &stub_hits,
                                  art::Ptr<recob::PFParticle> &stub_pfp) {
    // look up stuff
    const std::vector<art::Ptr<recob::Hit>> &hits = fSliceHits.at(slice.key());
    const std::vector<art::Ptr<recob::PFParticle>> &pfps = fSlicePFPs.at(slice.key());
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &pfp_hits = fSlicePFPHits.at(slice.key());
    
    TVector3 vertex_v(vertex.position().X(), vertex.position().Y(), vertex.position().Z());

    stub_hits = CollectHits(hits, vhit, vhit_hit, vertex.position(), geo, dprop);

    // Sort hits by their proximity to the vertex
    float vertex_w = geo->WireCoordinate(vertex.position(), vhit_hit.WireID());
    std::sort(stub_hits.begin(), stub_hits.end(), 
      [vertex_w](auto const &lhs, auto const &rhs) {
        return abs(lhs->WireID().Wire - vertex_w) < abs(rhs->WireID().Wire - vertex_w);});

    // Find which pfparticle to group this VertexStub with 
    stub_pfp = FindPFP(stub_hits, pfps, pfp_hits);

    sbn::Stub stub;

    // Save all the charges
    stub.charge.emplace_back();
    std::set<int> wires;
    for (unsigned i_hit = 0; i_hit < stub_hits.size(); i_hit++) {
      const recob::Hit &hit = *stub_hits[i_hit];
      float charge = fCaloAlg.ElectronsFromADCArea(hit.Integral(), hit.WireID().Plane) * fCaloAlg.LifetimeCorrection(dclock, dprop, hit.PeakTime(), 0.);
      stub.charge.back().push_back(charge);
      wires.insert(hit.WireID().Wire);
    }

    stub.vtx = vertex_v;
    stub.end = vhit.spXYZ;
    stub.pitch.push_back(vhit.charge / vhit.dqdx);
    stub.nwire.push_back(wires.size());
    stub.plane.push_back(vhit.wire);

    return stub;
}

