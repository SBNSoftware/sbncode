////////////////////////////////////////////////////////////////////////
// Class:       TrackSplitter
// Plugin Type: producer (art v3_02_06)
// File:        TrackSplitter_module.cc
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

#include <memory>
#include <optional>
#include "sbnobj/Common/Reco/MergedTrackInfo.hh"

namespace sbn {
  class TrackSplitter;
}


class sbn::TrackSplitter : public art::EDProducer {
public:
  explicit TrackSplitter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackSplitter(TrackSplitter const&) = delete;
  TrackSplitter(TrackSplitter&&) = delete;
  TrackSplitter& operator=(TrackSplitter const&) = delete;
  TrackSplitter& operator=(TrackSplitter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // input labels
  art::InputTag fTrackLabel;
  art::InputTag fMergedPFPLabel;

  // service holders
  const geo::GeometryCore *fGeo;
  std::optional<detinfo::DetectorPropertiesData> fDetProp;
  std::optional<detinfo::DetectorClocksData> fDetClock;

  // helper functions
  std::pair<int, float> ClosestTrajectoryPoint(const detinfo::DetectorPropertiesData &dprop, const recob::Track &trk, const recob::Hit &hit);
  std::vector<art::Ptr<recob::Hit>> SplitTrunkHits(
    const detinfo::DetectorPropertiesData &dprop, 
    const std::vector<art::Ptr<recob::Hit>> &trunk_hits,
    unsigned main_hit_ind,
    const std::vector<const recob::TrackHitMeta *> &trunk_hmetas,
    const recob::Track &trunk_trk,
    const recob::Track &branch_trk);

  std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<const recob::TrackHitMeta *>> OrganizeHits(
    const recob::Track &trk, const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<const recob::TrackHitMeta *> &thm, unsigned plane);

  // Setup struct for ret for the DeMergeHits thingy
  struct PairedHits {
    std::vector<recob::Hit> trunk_hits;
    std::vector<recob::Hit> branch_hits;
    std::vector<recob::TrackHitMeta> trunk_hmetas;
    std::vector<recob::TrackHitMeta> branch_hmetas;
  };

  PairedHits DeMergePlaneHits(
    const detinfo::DetectorPropertiesData &dprop,
    const recob::Track &trunk_trk,
    const recob::Track &new_branch_trk,
    const recob::Track &old_branch_trk,
    const std::vector<art::Ptr<recob::Hit>> &trunk_hits,
    const std::vector<art::Ptr<recob::Hit>> &branch_hits,
    const std::vector<const recob::TrackHitMeta *> &trunk_hmeta,
    const std::vector<const recob::TrackHitMeta *> &branch_hmeta,
    const sbn::MergedTrackInfo &info,
    unsigned plane);

  PairedHits DeMergeHits(
    const detinfo::DetectorPropertiesData &dprop,
    const recob::Track &trunk_trk,
    const recob::Track &new_branch_trk,
    const recob::Track &old_branch_trk,
    const std::vector<art::Ptr<recob::Hit>> &trunk_hits,
    const std::vector<art::Ptr<recob::Hit>> &branch_hits,
    const std::vector<const recob::TrackHitMeta *> &trunk_hmetas,
    const std::vector<const recob::TrackHitMeta *> &branch_hmetas,
    const sbn::MergedTrackInfo &info);

  recob::Track DeMergeTrack(const recob::Track &trunk, const recob::Track &branch, const sbn::MergedTrackInfo &info);
};

sbn::TrackSplitter::TrackSplitter(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrack")),
    fMergedPFPLabel(p.get<art::InputTag>("MergedPFPLabel", "mergeIdent"))
{
  produces<std::vector<recob::Track>>();
  produces<std::vector<recob::Hit>>();
  produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  produces<art::Assns<recob::PFParticle, recob::Track>>();
}

std::pair<int, float> sbn::TrackSplitter::ClosestTrajectoryPoint(const detinfo::DetectorPropertiesData &dprop, const recob::Track &trk, const recob::Hit &hit) {
  float hit_x = dprop.ConvertTicksToX(hit.PeakTime(), hit.WireID()); 
  float hit_w = hit.WireID().Wire * fGeo->WirePitch();

  // std::cout << "Hit time: " << hit.PeakTime() << std::endl;
  // std::cout << "Hit X: " << hit_x << std::endl;
  // std::cout << "Hit W: " << hit_w << std::endl;

  float dist = -1000.;
  int ret = -1;

  for (unsigned i_tp = trk.FirstValidPoint(); i_tp < trk.NumberTrajectoryPoints(); i_tp = trk.NextValidPoint(i_tp+1)) {
    TVector3 pos = trk.LocationAtPoint<TVector3>(i_tp);
    float pos_x = pos.X();           
    float pos_w = fGeo->WireCoordinate(trk.LocationAtPoint<geo::Point_t>(i_tp), hit.WireID()) * fGeo->WirePitch();
    float this_dist = sqrt((pos_x - hit_x) *(pos_x - hit_x) + (pos_w - hit_w) *(pos_w - hit_w));
    if (ret < 0 || this_dist < dist) {
      dist = this_dist;
      ret = i_tp;
    }
  }

  // if (ret == -1) std::cout << "No SP :(\n";
  // else {
  //   TVector3 pos = trk.LocationAtPoint<TVector3>(ret);
  //   std::cout << "Closest TP X: " << pos.X() << std::endl;
  //   std::cout << "Closest TP W: " << fGeo->WireCoordinate(pos.Y(), pos.Z(), hit.WireID()) * fGeo->WirePitch() << std::endl;
  // }

  return {ret, dist};
}

std::vector<art::Ptr<recob::Hit>> sbn::TrackSplitter::SplitTrunkHits(
    const detinfo::DetectorPropertiesData &dprop, 
    const std::vector<art::Ptr<recob::Hit>> &trunk_hits,
    unsigned main_hit_ind,
    const std::vector<const recob::TrackHitMeta *> &trunk_hmetas,
    const recob::Track &trunk_trk,
    const recob::Track &branch_trk) {

  std::vector<art::Ptr<recob::Hit>> split_trunk_hits;

  geo::WireID wire = trunk_hits[main_hit_ind]->WireID();

  for (unsigned i = 0; i < trunk_hits.size(); i++) {
    art::Ptr<recob::Hit> h = trunk_hits[i];
    const recob::TrackHitMeta *hmeta = trunk_hmetas[i];
    if (i != main_hit_ind // different hit
      && (hmeta->Index() == std::numeric_limits<unsigned int>::max() || !trunk_trk.HasValidPoint(hmeta->Index())) // not used in trunk traj. 
      && h->WireID() == wire) { // on same wire

      // Found a hit! So there are two hits on the trunk on this wire
      // and one may belong on the branch
      //
      // We'll say this hit "belongs" to the branch if it closer to the projection
      // of the branch on the plane

      // std::cout << "TRUNK\n";
      std::pair<int, float> trunk_dist = ClosestTrajectoryPoint(dprop, trunk_trk, *h);
      // std::cout << "BRANCH\n";
      std::pair<int, float> branch_dist = ClosestTrajectoryPoint(dprop, branch_trk, *h);

      // std::cout << "Possible split trunk hit on wire: " << h->WireID().Wire << ". Trunk dist: " << trunk_dist.second << " branch dist: " << branch_dist.second << std::endl;

      if (trunk_dist.first >= 0 && branch_dist.first > 0 && branch_dist.second < trunk_dist.second) {
        split_trunk_hits.push_back(h);
      }
    }
  }

  return split_trunk_hits;
}

sbn::TrackSplitter::PairedHits sbn::TrackSplitter::DeMergePlaneHits(
    const detinfo::DetectorPropertiesData &dprop, 
    const recob::Track &trunk_trk,
    const recob::Track &new_branch_trk,
    const recob::Track &old_branch_trk,
    const std::vector<art::Ptr<recob::Hit>> &trunk_hits,
    const std::vector<art::Ptr<recob::Hit>> &branch_hits,
    const std::vector<const recob::TrackHitMeta *> &trunk_hmeta,
    const std::vector<const recob::TrackHitMeta *> &branch_hmeta,
    const sbn::MergedTrackInfo &info,
    unsigned plane) {

  sbn::TrackSplitter::PairedHits ret;

  unsigned n_new_branch_trajpoints = new_branch_trk.NumberTrajectoryPoints();
  unsigned n_old_branch_trajpoints = old_branch_trk.NumberTrajectoryPoints();
  unsigned n_new_points = n_new_branch_trajpoints - n_old_branch_trajpoints;

  unsigned n_hits = trunk_hits.size();

  // std::cout << "De-merging plane: " << plane << std::endl;

  // std::cout << "Branch start: " << info.branch_start << std::endl;

  // consider the energy depositions in the trunk before the branch starts
  // Consider the hits in the trunk before the branch starts
  unsigned i_hit_trunk = 0;
  for (; i_hit_trunk < n_hits; i_hit_trunk++) {
    const recob::Hit &hit = *trunk_hits[i_hit_trunk];
    const recob::TrackHitMeta *thm = trunk_hmeta[i_hit_trunk];

    // ignore hits not used in the trunk trajectory
    if (thm->Index() == std::numeric_limits<unsigned int>::max() || !trunk_trk.HasValidPoint(thm->Index())) continue;

    // Are we past the track start?
    TVector3 loc = trunk_trk.LocationAtPoint<TVector3>(thm->Index());
    bool this_point_past_branch_start = (loc - info.vertex).Dot(info.direction) > info.branch_start; 
    // bool this_point_past_branch_start = info.trunk_wire_direction_is_ascending[plane] ? 
    //    (int)hit.WireID().Wire >= info.branch_wire_start[plane] :
    //    (int)hit.WireID().Wire <= info.branch_wire_start[plane];

    // std::cout << "This wire: " << hit.WireID().Wire << std::endl;
    // std::cout << "This point: " << (loc - info.vertex).Dot(info.direction) << std::endl;

    // lookup if the branch has a hit on this wire 
    std::vector<art::Ptr<recob::Hit>> branch_wire_hits;
    std::vector<const recob::TrackHitMeta *> branch_wire_hmetas;
    for (unsigned i_hit_branch = 0; i_hit_branch < branch_hits.size(); i_hit_branch++) {
      art::Ptr<recob::Hit> b_hit = branch_hits[i_hit_branch];
      if (b_hit->WireID() == hit.WireID()) {
        branch_wire_hits.push_back(b_hit);
        branch_wire_hmetas.push_back(branch_hmeta[i_hit_branch]);
      }
    } 

    // std::cout << "Number of branch hits on this wire: " << branch_wire_hits.size() << std::endl;

    bool valid_tp = false;
    unsigned i_hit_branch = 0;
    if (branch_wire_hits.size()) {
      // Check any hits on the branch already have a valid trajectory point.
      for (; i_hit_branch < branch_wire_hits.size(); i_hit_branch++) {
        if (branch_wire_hmetas[i_hit_branch]->Index() != std::numeric_limits<unsigned int>::max() && 
            old_branch_trk.HasValidPoint(branch_wire_hmetas[i_hit_branch]->Index())) {
          valid_tp = true;
          break;
        }
      }
    }

    // std::cout << "Valid traj-point: " << valid_tp << std::endl;

    // Past the start of the branch track! We should be all set now
    if (this_point_past_branch_start && valid_tp) {
      ret.trunk_hits.push_back(hit);
      ret.trunk_hmetas.push_back(*thm);
      break;
    }

    // No hit lurking on the branch -- we'll have to split the trunk
    // First see if we can find a stray hit on the trunk
    if (!branch_wire_hits.size()) {
      branch_wire_hits = SplitTrunkHits(dprop, trunk_hits, i_hit_trunk, trunk_hmeta, trunk_trk, new_branch_trk);
    }

    // std::cout << "Number split trunk hits: " << branch_wire_hits.size() << std::endl;

    // Ok! Now we have all the information we need to determine what hit to assign
    // to the branch on this wire. There are three possiblities:
    //
    // 1. The branch has a valid hit/traj-point on this wire.
    // 2. There is a spare hit on the branch or trunk we can associate to the branch trajectory
    // 3. There are no spare hits on the branch and we need to split the trunk hit.
    //
    // These three possibilties are handled in order below: 
    if (valid_tp) {
      // Valid hit on this wire! Save it and go onto the next trunk hit
      //
      // Make a new TrackHitMeta that updates the traj point to the new correct place
      recob::TrackHitMeta new_branch_meta(branch_wire_hmetas[i_hit_branch]->Index() + n_new_points);
      ret.branch_hits.push_back(*branch_wire_hits[i_hit_branch]);
      ret.branch_hmetas.push_back(new_branch_meta);

      ret.trunk_hits.push_back(hit);
      ret.trunk_hmetas.push_back(*thm);
    }
    // Found a stray hit! We can work with that one
    else if (branch_wire_hits.size()) {
      art::Ptr<recob::Hit> new_branch_hit = *std::max_element(branch_wire_hits.begin(), branch_wire_hits.end(),
        [this, dprop, new_branch_trk](auto const &lhs, auto const &rhs) { 
          return ClosestTrajectoryPoint(dprop, new_branch_trk, *lhs).second < ClosestTrajectoryPoint(dprop, new_branch_trk, *rhs).second;
        });
      unsigned i_tp = (unsigned)ClosestTrajectoryPoint(dprop, new_branch_trk, *new_branch_hit).first;

      recob::TrackHitMeta new_meta(i_tp);

      // save this hit
      ret.branch_hits.push_back(*new_branch_hit);
      ret.branch_hmetas.push_back(new_meta);

      // save the main hit for the trunk
      ret.trunk_hits.push_back(hit);
      ret.trunk_hmetas.push_back(*thm);
    }
    // No stray hit -- we'll have to split the tracks single hit on this wire in 2
    else {
      recob::Hit halfhit(hit.Channel(), hit.StartTick(), hit.EndTick(), hit.PeakTime(),
        hit.SigmaPeakTime(), hit.RMS(), 
        hit.PeakAmplitude() / 2., hit.SigmaPeakAmplitude() / 2., hit.SummedADC() / 2.,
        hit.Integral() / 2., hit.SigmaIntegral() / 2., 
        hit.Multiplicity(), hit.LocalIndex(), hit.GoodnessOfFit(), hit.DegreesOfFreedom(), hit.View(),
        hit.SignalType(), hit.WireID());

      unsigned branch_tp = ClosestTrajectoryPoint(dprop, new_branch_trk, halfhit).first;
      recob::TrackHitMeta branch_meta(branch_tp);

      // save the half-hit for the branch and the trunk
      ret.branch_hits.push_back(halfhit);
      ret.branch_hmetas.push_back(branch_meta);

      ret.trunk_hits.push_back(halfhit);
      ret.trunk_hmetas.push_back(*thm);
    }
  } 

  // fill in the rest of the trunk hits
  for (; i_hit_trunk < n_hits; i_hit_trunk++) {
    const recob::Hit &hit = *trunk_hits[i_hit_trunk];
    const recob::TrackHitMeta *thm = trunk_hmeta[i_hit_trunk];

    // ignore hits not used in the trunk trajectory
    if (thm->Index() == std::numeric_limits<unsigned int>::max() || !trunk_trk.HasValidPoint(thm->Index())) continue;

    ret.trunk_hits.push_back(hit);
    ret.trunk_hmetas.push_back(*thm);
  }

  // fill in the rest of the branch hits
  unsigned n_split_hits = ret.branch_hits.size();
  for (unsigned i_hit_branch = 0; i_hit_branch < branch_hits.size(); i_hit_branch++) {

    // don't save invalid hits
    if (branch_hmeta[i_hit_branch]->Index() == std::numeric_limits<unsigned int>::max() ||
        !old_branch_trk.HasValidPoint(branch_hmeta[i_hit_branch]->Index())) {
      continue;
    }

    // don't save hits that were already added to the track
    bool found = false;
    for (unsigned j_hit_branch = 0; j_hit_branch < n_split_hits; j_hit_branch++) {
      if (ret.branch_hits[j_hit_branch].WireID() == branch_hits[i_hit_branch]->WireID()) {
        found = true;
        break;
      }
    }
    if (found) continue;

    // Make a new TrackHitMeta that updates the traj point to the new correct place
    recob::TrackHitMeta new_branch_meta(branch_hmeta[i_hit_branch]->Index() + n_new_points);
    ret.branch_hits.push_back(*branch_hits[i_hit_branch]);
    ret.branch_hmetas.push_back(new_branch_meta);
  }

  return ret;
}

std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<const recob::TrackHitMeta *>> sbn::TrackSplitter::OrganizeHits(
    const recob::Track &trk, const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<const recob::TrackHitMeta *> &thm, unsigned plane) {

  std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<const recob::TrackHitMeta *>> ret;

  // Collect the hits on the provided plane
  for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
    if (hits[i_hit]->WireID().Plane == plane) {
      ret.first.push_back(hits[i_hit]);
      ret.second.push_back(thm[i_hit]);
    }
  }

  // Now decide if the hit order needs to be reversed
  //
  // DeMergePlaneHits will assume hits are ordered start -> end

  int first_valid_hit = -1;
  int last_valid_hit = -1;

  for (unsigned i_hit = 0; i_hit < ret.first.size(); i_hit++) {
    if (ret.second[i_hit]->Index() != std::numeric_limits<unsigned int>::max() && trk.HasValidPoint(ret.second[i_hit]->Index())) {
      first_valid_hit = i_hit;
      break;
    }
  }

  for (int i_hit = ret.first.size()-1; i_hit >= 0; i_hit--) {
    if (ret.second[i_hit]->Index() != std::numeric_limits<unsigned int>::max() && trk.HasValidPoint(ret.second[i_hit]->Index())) {
      last_valid_hit = i_hit;
      break;
    }
  }

  bool reverse_hits = false;

  // If there are two valid hits/traj-point's, then we can decide to re-order using which
  // is closer to the start
  if (first_valid_hit >= 0 && last_valid_hit >= 0 && first_valid_hit != last_valid_hit) {
    reverse_hits = \
      (trk.LocationAtPoint<TVector3>(ret.second[first_valid_hit]->Index()) - trk.Start<TVector3>()).Mag() >
      (trk.LocationAtPoint<TVector3>(ret.second[last_valid_hit]->Index()) -  trk.Start<TVector3>()).Mag();
  }
  // Otherwise, try to see if we can at least determine whether wires and traj-points are ascending or descending
  //
  // If not enough info, assume the order is correct
  else if (ret.first.size() && trk.NumberTrajectoryPoints() >= 2) {
    TVector3 start = trk.Start<TVector3>();
    TVector3 start_p1 = trk.LocationAtPoint<TVector3>(trk.NextValidPoint(trk.FirstValidPoint()+1));

    geo::PlaneID thisplane(0, 0, plane);
    bool trk_is_ascending = fGeo->WireCoordinate(trk.Start<geo::Point_t>(), thisplane) < 
                           fGeo->WireCoordinate(trk.LocationAtPoint<geo::Point_t>(trk.NextValidPoint(trk.FirstValidPoint()+1)), thisplane);
    int wire0 = ret.first[0]->WireID().Wire;
    int wire1 = -1;
    for (unsigned i_hit = 1; i_hit < ret.first.size(); i_hit++) {
      if ((int)ret.first[i_hit]->WireID().Wire != wire0) {
        wire1 = ret.first[i_hit]->WireID().Wire;
        break;
      }
    }
    bool wire_is_ascending = (wire1 >= 0) ? wire1 > wire0 : trk_is_ascending;

    reverse_hits = trk_is_ascending != wire_is_ascending;
  }

  if (reverse_hits) {
    std::reverse(ret.first.begin(), ret.first.end());
    std::reverse(ret.second.begin(), ret.second.end());
  }

  return ret;
}

sbn::TrackSplitter::PairedHits sbn::TrackSplitter::DeMergeHits(
    const detinfo::DetectorPropertiesData &dprop,
    const recob::Track &trunk_trk,
    const recob::Track &new_branch_trk,
    const recob::Track &old_branch_trk,
    const std::vector<art::Ptr<recob::Hit>> &trunk_hits,
    const std::vector<art::Ptr<recob::Hit>> &branch_hits,
    const std::vector<const recob::TrackHitMeta *> &trunk_hmetas,
    const std::vector<const recob::TrackHitMeta *> &branch_hmetas,
    const sbn::MergedTrackInfo &info) {

  sbn::TrackSplitter::PairedHits ret;

  for (unsigned plane = 0; plane < 3; plane++) {
    // collect the valid hits on each plane
    std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<const recob::TrackHitMeta *>> trunkPlaneHits = \
      OrganizeHits(trunk_trk, trunk_hits, trunk_hmetas, plane);

    std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<const recob::TrackHitMeta *>> branchPlaneHits = \
      OrganizeHits(old_branch_trk, branch_hits, branch_hmetas, plane);

    // Split
    sbn::TrackSplitter::PairedHits paired_hits = DeMergePlaneHits(dprop, trunk_trk, new_branch_trk, old_branch_trk, trunkPlaneHits.first, branchPlaneHits.first, 
        trunkPlaneHits.second, branchPlaneHits.second, info, plane);

    // Save
    ret.trunk_hits.insert(ret.trunk_hits.end(), paired_hits.trunk_hits.begin(), paired_hits.trunk_hits.end());
    ret.trunk_hmetas.insert(ret.trunk_hmetas.end(), paired_hits.trunk_hmetas.begin(), paired_hits.trunk_hmetas.end());
    ret.branch_hits.insert(ret.branch_hits.end(), paired_hits.branch_hits.begin(), paired_hits.branch_hits.end());
    ret.branch_hmetas.insert(ret.branch_hmetas.end(), paired_hits.branch_hmetas.begin(), paired_hits.branch_hmetas.end());
  }

  return ret;

}

recob::Track sbn::TrackSplitter::DeMergeTrack(const recob::Track &trunk, const recob::Track &branch, const sbn::MergedTrackInfo &info) {
  // add points from the trunk-track
  std::vector<TVector3> add_p;
  std::vector<TVector3> add_m;
  std::vector<recob::TrajectoryPointFlags> add_f;

  for (unsigned i = 0; i < trunk.NumberTrajectoryPoints(); i++) {
    if (!trunk.HasValidPoint(i)) continue;

    TVector3 p = trunk.LocationAtPoint<TVector3>(i);
    float p_proj = (p - info.vertex).Dot(info.direction);

    if (p_proj >= info.branch_start) break;

    add_p.push_back(p);
    add_m.push_back(trunk.DirectionAtPoint<TVector3>(i));
    add_f.push_back(trunk.FlagsAtPoint(i));
  }

  // smooth the addded points to avoid big jump between trunk and branch
  if (add_p.size() >= 2) {
    std::vector<TVector3> add_p_old = add_p;

    float dist = info.branch_start;  
    TVector3 diff = branch.Start<TVector3>() - add_p_old[add_p_old.size()-1];

    for (unsigned i = 1; i < add_p_old.size(); i++) {
      float this_dist = (add_p_old[i] - info.vertex).Dot(info.direction);
      add_p[i] = add_p_old[i] + diff * (this_dist / dist);
    }

    for (int i = 0; i < (int)add_m.size() - 1; i++) {
      add_m[i] = (add_p[i+1] - add_p[i]).Unit();
      if (i > 0) add_m[i] = (add_m[i] + (add_p[i] - add_p[i-1]).Unit()).Unit();
    }
    add_m[add_p.size()-1] = ((branch.Start<TVector3>() - add_p[add_p.size()-1]).Unit() +
      (add_p[add_p.size()-1] - add_p[add_p.size()-2])).Unit();

    // TODO: fix what the initial direction gets set to
    add_m[0] = (branch.Start<TVector3>() - add_p[0]).Unit();
  }

  /*
  TVector3 p = trunk.Start<TVector3>();
  float p_proj = (p - info.vertex).Dot(info.direction);

  if (p_proj < info.branch_start) {
    add_p.push_back(p);
    add_m.push_back((branch.Start<TVector3>() - p).Unit());
    add_f.push_back(trunk.FlagsAtPoint(trunk.FirstValidPoint()));
  }
  */

  // Save the new points
  recob::Track::Positions_t positions; 
  recob::Track::Momenta_t momenta;
  recob::Track::Flags_t flags = std::move(add_f);

  for (unsigned i = 0; i < add_p.size(); i++) {
    positions.emplace_back(add_p[i].X(), add_p[i].Y(), add_p[i].Z());
    momenta.emplace_back(add_m[i].X(), add_m[i].Y(), add_m[i].Z());
  }

  // copy the starting info from the branch-track
  for (unsigned i_tp = 0; i_tp < branch.NumberTrajectoryPoints(); i_tp += 1) {
    positions.push_back(branch.LocationAtPoint(i_tp));
    momenta.push_back(branch.DirectionAtPoint(i_tp));
    flags.push_back(branch.FlagsAtPoint(i_tp));
  }

  recob::Track::SMatrixSym55 cov_start = trunk.StartCovariance();
  recob::Track::SMatrixSym55 cov_end = branch.EndCovariance();
  return recob::Track(std::move(positions), std::move(momenta), std::move(flags),
    branch.HasMomentum(), branch.ParticleId(), branch.Chi2(), branch.Ndof(), std::move(cov_start), std::move(cov_end), branch.ID());
}

void sbn::TrackSplitter::produce(art::Event& e)
{
  // update services
  fGeo = lar::providerFrom<geo::Geometry>();

  auto const &dclock = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e); 
  auto const &dprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, dclock);

  // output data products
  std::unique_ptr<art::Assns<recob::PFParticle, recob::Track>> trkAssn(new art::Assns<recob::PFParticle, recob::Track>);
  std::unique_ptr<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>> hitAssn(new art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>);
  std::unique_ptr<std::vector<recob::Track>> outTracks(new std::vector<recob::Track>);
  std::unique_ptr<std::vector<recob::Hit>> outHits(new std::vector<recob::Hit>);

  art::PtrMaker<recob::Track> trkPtrMaker{e};
  art::PtrMaker<recob::Hit> hitPtrMaker{e};

  // input data
  art::Handle<std::vector<sbn::MergedTrackInfo>> mergedtrack_handle;
  e.getByLabel(fMergedPFPLabel, mergedtrack_handle);

  std::vector<art::Ptr<sbn::MergedTrackInfo>> mergedtracks;
  art::fill_ptr_vector(mergedtracks, mergedtrack_handle);

  art::FindManyP<recob::PFParticle> mergedPFPs(mergedtracks, e, fMergedPFPLabel);

  for (unsigned i_mrg = 0; i_mrg < mergedtracks.size(); i_mrg++) {
    const sbn::MergedTrackInfo &merge_info = *mergedtracks[i_mrg];
    const std::vector<art::Ptr<recob::PFParticle>> pfps = mergedPFPs.at(i_mrg);
    art::FindManyP<recob::Track> fmTracks(pfps, e, fTrackLabel);

    assert(pfps.size() == 2);

    art::Ptr<recob::PFParticle> trunk_pfp;
    art::Ptr<recob::PFParticle> branch_pfp;
    art::Ptr<recob::Track> trunk_trk;
    art::Ptr<recob::Track> branch_trk;
    if ((unsigned)merge_info.trunk == pfps[0]->Self()) {
      trunk_pfp = pfps[0];
      trunk_trk = fmTracks.at(0).at(0);
      branch_pfp = pfps[1];
      branch_trk = fmTracks.at(1).at(0);
    }
    else {
      trunk_pfp = pfps[1];
      trunk_trk = fmTracks.at(1).at(0);
      branch_pfp = pfps[0];
      branch_trk = fmTracks.at(0).at(0);
    }

    // std::cout << "Splitting branch: " << branch_pfp->Self() << " from trunk: " << trunk_pfp->Self() << std::endl;

    // De-merge with the provided match
    recob::Track demerge = DeMergeTrack(*trunk_trk, *branch_trk, merge_info);

    // look-up the track hits
    art::FindManyP<recob::Hit, recob::TrackHitMeta> hits({trunk_trk, branch_trk}, e, fTrackLabel);

    sbn::TrackSplitter::PairedHits demerged_hits = DeMergeHits(dprop, 
      *trunk_trk, demerge, *branch_trk, 
      hits.at(0), hits.at(1), hits.data(0), hits.data(1), 
      merge_info);

    // Save
    outTracks->push_back(demerge);
    art::Ptr<recob::Track> branchTrackPtr = trkPtrMaker(outTracks->size()-1);   
    trkAssn->addSingle(branch_pfp, branchTrackPtr);

    for (unsigned i_hit = 0; i_hit < demerged_hits.branch_hits.size(); i_hit++) {
      outHits->push_back(demerged_hits.branch_hits[i_hit]);
      art::Ptr<recob::Hit> thisHitPtr = hitPtrMaker(outHits->size()-1);
      hitAssn->addSingle(branchTrackPtr, thisHitPtr, demerged_hits.branch_hmetas[i_hit]);
    }

    // Also make a new track for the trunk
    recob::Track trunk_copy = *trunk_trk;
    outTracks->push_back(trunk_copy);
    art::Ptr<recob::Track> trunkTrackPtr = trkPtrMaker(outTracks->size()-1);
    trkAssn->addSingle(trunk_pfp, trunkTrackPtr);

    for (unsigned i_hit = 0; i_hit < demerged_hits.trunk_hits.size(); i_hit++) {
      outHits->push_back(demerged_hits.trunk_hits[i_hit]);
      art::Ptr<recob::Hit> thisHitPtr = hitPtrMaker(outHits->size()-1);
      hitAssn->addSingle(trunkTrackPtr, thisHitPtr, demerged_hits.trunk_hmetas[i_hit]);
    }

  } // end iterate over merged PFParticles
  
  // Save into event
  e.put(std::move(outTracks));
  e.put(std::move(trkAssn));
  e.put(std::move(outHits));
  e.put(std::move(hitAssn));
}

DEFINE_ART_MODULE(sbn::TrackSplitter)
