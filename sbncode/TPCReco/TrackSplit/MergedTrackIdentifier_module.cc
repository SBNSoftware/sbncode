////////////////////////////////////////////////////////////////////////
// Class:       MergedTrackIdentifier
// Plugin Type: producer (art v3_02_06)
// File:        MergedTrackIdentifier_module.cc
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
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <memory>

#include "sbnobj/Common/Reco/MergedTrackInfo.hh"

namespace sbn {
  class MergedTrackIdentifier;
}


class sbn::MergedTrackIdentifier : public art::EDProducer {
public:
  explicit MergedTrackIdentifier(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MergedTrackIdentifier(MergedTrackIdentifier const&) = delete;
  MergedTrackIdentifier(MergedTrackIdentifier&&) = delete;
  MergedTrackIdentifier& operator=(MergedTrackIdentifier const&) = delete;
  MergedTrackIdentifier& operator=(MergedTrackIdentifier&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fPFPLabel;
  art::InputTag fTrackLabel;

  // helper functions
  sbn::MergedTrackInfo BuildTrackInfo(
    const geo::GeometryCore *geo, const recob::Vertex &vtx, 
    const recob::PFParticle &A, const recob::Track &Atrk, 
    const recob::PFParticle &B, const recob::Track &Btrk, bool assume_a_is_trunk=false);

};

sbn::MergedTrackIdentifier::MergedTrackIdentifier(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFPLabel(p.get<art::InputTag>("PFPLabel", "pandora")),
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrack"))
{
  produces<std::vector<sbn::MergedTrackInfo>>();
  produces<art::Assns<sbn::MergedTrackInfo, recob::PFParticle>>();
}


// static helper
TVector3 MeanDirection(const recob::Vertex &vtx, const recob::Track &trk) {
  TVector3 start(vtx.position().x(), vtx.position().y(), vtx.position().z());

  TVector3 avg_dir(0., 0., 0.);

  for (unsigned i_tp = trk.FirstValidPoint(); i_tp < trk.NumberTrajectoryPoints(); i_tp = trk.NextValidPoint(i_tp+1)) {
    TVector3 p = trk.LocationAtPoint<TVector3>(i_tp);
    avg_dir += (p - start).Unit();
  }

  return avg_dir.Unit();

}

sbn::MergedTrackInfo sbn::MergedTrackIdentifier::BuildTrackInfo(
    const geo::GeometryCore *geo, const recob::Vertex &vtx, 
    const recob::PFParticle &A, const recob::Track &Atrk, 
    const recob::PFParticle &B, const recob::Track &Btrk, bool assume_a_is_trunk) {

  sbn::MergedTrackInfo ret; 

  ret.vertex = TVector3(vtx.position().x(), vtx.position().y(), vtx.position().z());

  ret.direction = (MeanDirection(vtx, Atrk) + MeanDirection(vtx, Btrk)).Unit();

  bool a_is_trunk = assume_a_is_trunk || ret.direction.Dot(Atrk.Start<TVector3>() - ret.vertex) < ret.direction.Dot(Btrk.Start<TVector3>() - ret.vertex);

  const recob::PFParticle &trunk =  (a_is_trunk) ? A : B;
  const recob::PFParticle &branch = (a_is_trunk) ? B : A;

  const recob::Track &trunk_trk =  (a_is_trunk) ? Atrk : Btrk;
  const recob::Track &branch_trk = (a_is_trunk) ? Btrk : Atrk;

  ret.trunk = trunk.Self();
  ret.branch = branch.Self();

  ret.branch_start = ret.direction.Dot(branch_trk.Start<TVector3>() - ret.vertex);
  ret.trunk_start =  ret.direction.Dot(trunk_trk.Start<TVector3>()  - ret.vertex);

  // TVector3 branch_start = branch_trk.LocationAtPoint<TVector3>(branch_trk.FirstValidPoint());
  // for (unsigned i = 0; i < 3; i++) {
  //   ret.branch_wire_start[i] = geo->WireCoordinate(branch_start.Y(), branch_start.Z(), i, 0, 0);

  //  TVector3 trunk_start = trunk_trk.LocationAtPoint<TVector3>(trunk_trk.FirstValidPoint());
  //  TVector3 trunk_start_p1 = trunk_trk.LocationAtPoint<TVector3>(trunk_trk.NextValidPoint(trunk_trk.FirstValidPoint()+1));
    
    // check if the tracks are ascending or descending
  //  ret.trunk_wire_direction_is_ascending[i] = 
  //    geo->WireCoordinate(trunk_start.Y(), trunk_start.Z(), i, 0, 0) <
  //    geo->WireCoordinate(trunk_start_p1.Y(), trunk_start_p1.Z(), i, 0, 0);
  // }

  bool set = false;
  float trunk_min = -1;
  float trunk_max = -1;
  
  for (unsigned i = trunk_trk.FirstValidPoint(); i < trunk_trk.NumberTrajectoryPoints(); i = trunk_trk.NextValidPoint(i+1)) {
    float proj = ret.direction.Dot(trunk_trk.DirectionAtPoint<TVector3>(i) - ret.vertex);
    if (!set || proj < trunk_min) trunk_min = proj;
    if (!set || proj > trunk_max) trunk_max = proj;
    set = true;
  }

  unsigned n_point = 0;
  unsigned n_overlap_point = 0;
  for (unsigned i = branch_trk.FirstValidPoint(); i < branch_trk.NumberTrajectoryPoints(); i = branch_trk.NextValidPoint(i+1)) {
    float proj = ret.direction.Dot(branch_trk.DirectionAtPoint<TVector3>(i) - ret.vertex);
    if (proj > trunk_min && proj < trunk_max) {
      n_overlap_point ++;
    }
    n_point ++;
  }

  ret.branch_overlap = ((float)n_overlap_point) / n_point;
  

  return ret;
}

void sbn::MergedTrackIdentifier::produce(art::Event& e)
{
  // output data products
  std::unique_ptr<std::vector<sbn::MergedTrackInfo>> infos(new std::vector<sbn::MergedTrackInfo>);
  std::unique_ptr<art::Assns<sbn::MergedTrackInfo, recob::PFParticle>> assn(new art::Assns<sbn::MergedTrackInfo, recob::PFParticle>);

  art::PtrMaker<sbn::MergedTrackInfo> infoPtrMaker {e};

  // input data
  art::Handle<std::vector<recob::Slice>> slice_handle;
  e.getByLabel(fPFPLabel, slice_handle);

  std::vector<art::Ptr<recob::Slice>> slices;
  art::fill_ptr_vector(slices, slice_handle);

  art::FindManyP<recob::PFParticle> slicePFPs(slices, e, fPFPLabel);

  // services
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();

  for (unsigned i_slc = 0; i_slc < slices.size(); i_slc++) {
    // collect the PFPs
    const std::vector<art::Ptr<recob::PFParticle>> &this_slice_pfps = slicePFPs.at(i_slc);

    art::FindManyP<recob::Vertex> slicePFPVtxs(this_slice_pfps, e, fPFPLabel);

    art::Ptr<recob::Vertex> vtx;
    art::Ptr<recob::PFParticle> primary;
    // get the primary particle in each Slice
    for (unsigned i_pfp = 0; i_pfp < this_slice_pfps.size(); i_pfp++) {
      if (this_slice_pfps[i_pfp]->IsPrimary() && slicePFPVtxs.at(i_pfp).size()) {
        primary = this_slice_pfps[i_pfp];
        vtx = slicePFPVtxs.at(i_pfp).at(0);
        break;
      }
    }

    // no primary PFParticle -- ignore 
    if (!primary) continue;

    std::vector<art::Ptr<recob::PFParticle>> topTracks;

    // primary is a Neutrino -- will have a number of daughters -- which we call primary 
    if (lar_pandora::LArPandoraHelper::IsNeutrino(primary)) {
      for (unsigned d: primary->Daughters()) {
        for (unsigned i_pfp = 0; i_pfp < this_slice_pfps.size(); i_pfp++) {
          if (this_slice_pfps[i_pfp]->Self() == d) {
            if (lar_pandora::LArPandoraHelper::IsTrack(this_slice_pfps[i_pfp])) topTracks.push_back(this_slice_pfps[i_pfp]);
            break;
          }
        }
      }
    }
    // primary is a track -- this is the only top-level-track
    else if (lar_pandora::LArPandoraHelper::IsTrack(primary)) {
      topTracks.push_back(primary);
    }

    art::FindManyP<recob::Track> topTrackTracks(topTracks, e, fTrackLabel);

    std::map<art::Ptr<recob::PFParticle>, std::vector<std::pair<art::Ptr<recob::PFParticle>, sbn::MergedTrackInfo>>> possible_merge_parents;
    std::map<art::Ptr<recob::PFParticle>, std::vector<std::pair<art::Ptr<recob::PFParticle>, sbn::MergedTrackInfo>>> possible_merge_children;

    // First iterate over all pairs of top-level tracks
    for (unsigned i = 0; i < topTracks.size(); i++) {
      for (unsigned j = i+1; j < topTracks.size(); j++) {
        if (!topTrackTracks.at(i).size() || !topTrackTracks.at(j).size()) continue;

        sbn::MergedTrackInfo info = BuildTrackInfo(geo, *vtx, *topTracks[i], *topTrackTracks.at(i).at(0), *topTracks[j], *topTrackTracks.at(j).at(0));

        // std::cout << "CONSIDERING MATCH: " << topTracks[i]->Self() << " paired " << topTracks[j]->Self() << " overlap: " << info.branch_overlap << " trunk: " << info.trunk << std::endl;

        // only take cases with 50% overlap 
        if (info.branch_overlap < 0.5) continue;

        // std::cout << "Overlaps!\n";

        if ((unsigned)info.trunk == topTracks[i]->Self()) {
          possible_merge_parents[topTracks[i]].push_back({topTracks[j], info});
          possible_merge_children[topTracks[j]].push_back({topTracks[i], info});

        }
        else {
          possible_merge_parents[topTracks[j]].push_back({topTracks[i], info});
          possible_merge_children[topTracks[i]].push_back({topTracks[j], info});
        }

      }
    }

    // also get cases with a top-track w/ a daughter track
    for (unsigned i = 0; i < topTracks.size(); i++) {
      std::vector<art::Ptr<recob::PFParticle>> daughters;
      for (unsigned d: topTracks[i]->Daughters()) {
        for (unsigned i_pfp = 0; i_pfp < this_slice_pfps.size(); i_pfp++) {
          if (this_slice_pfps[i_pfp]->Self() == d) {
            if (lar_pandora::LArPandoraHelper::IsTrack(this_slice_pfps[i_pfp])) daughters.push_back(this_slice_pfps[i_pfp]);
            break;
          }
        }
      }

      art::FindManyP<recob::Track> daughterTracks(daughters, e, fTrackLabel);
      for (unsigned j = 0; j < daughters.size(); j++) {
        if (!topTrackTracks.at(i).size() || !daughterTracks.at(j).size()) continue;

        sbn::MergedTrackInfo info = BuildTrackInfo(geo, *vtx, *topTracks[i], *topTrackTracks.at(i).at(0), *daughters[j], *daughterTracks.at(j).at(0), true);

        // std::cout << "CONSIDERING MATCH: " << topTracks[i]->Self() << " daughter " << daughters[j]->Self() << " overlap: " << info.branch_overlap << " trunk: " << info.trunk << std::endl;

        if (info.branch_overlap < 0.5) continue;

        possible_merge_parents[topTracks[i]].push_back({daughters[j], info});
        possible_merge_children[daughters[j]].push_back({topTracks[j], info});
      }
    }

    // Evaluate all the possible merges
    // Each track should only have one other possible merge
    for (auto const &pair: possible_merge_parents) {
      // std::cout << "Possible Merged: " << pair.first->Self()  << " N parents: " << possible_merge_children[pair.first].size() << std::endl;
      // for (auto const &pair_v: pair.second) {
        // std::cout << "Child: " << pair_v.first->Self() << " with N parents: " << possible_merge_children[pair_v.first].size() << std::endl;
        // std::cout << "Info: trunk: " << pair_v.second.trunk << " overlap: " << pair_v.second.branch_overlap << std::endl;
      // }

      if (pair.second.size() == 1 && possible_merge_children[pair.first].size() == 0) {
        art::Ptr<recob::PFParticle> branch = pair.second[0].first;
        if (possible_merge_parents[branch].size() == 0 && possible_merge_children[branch].size() == 1) {
          infos->push_back(pair.second[0].second);
          art::Ptr<sbn::MergedTrackInfo> infoPtr = infoPtrMaker(infos->size()-1);
          assn->addSingle(infoPtr, pair.first);
          assn->addSingle(infoPtr, branch);
        }
      }
    }
  } // end iterate over slices

  e.put(std::move(infos));
  e.put(std::move(assn));

}

DEFINE_ART_MODULE(sbn::MergedTrackIdentifier)
