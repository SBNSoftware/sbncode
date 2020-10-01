////////////////////////////////////////////////////////////////////////
// Class:       CRTHitMatch
// Plugin Type: producer (art v3_02_06)
// File:        CRTHitMatch_module.cc
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
#include "lardataobj/AnalysisBase/T0.h"

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"

#include "Products/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

#include <memory>

namespace sbn {
  class CRTHitMatch;
}


class sbn::CRTHitMatch : public art::EDProducer {
public:
  explicit CRTHitMatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTHitMatch(CRTHitMatch const&) = delete;
  CRTHitMatch(CRTHitMatch&&) = delete;
  CRTHitMatch& operator=(CRTHitMatch const&) = delete;
  CRTHitMatch& operator=(CRTHitMatch&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  sbnd::CRTT0MatchAlg fMatchAlg;
  art::InputTag fCRTHitLabel;
  art::InputTag fTrackLabel;
  int fTSMode;
  float fTimeCorrection;
  float fMinTrackLength;

};

sbnd::crt::CRTHit SBN2SBNDCrtHit(const sbn::crt::CRTHit &inp) {
  sbnd::crt::CRTHit ret;

  ret.peshit = inp.peshit;
  ret.ts0_s = inp.ts0_s;
  ret.ts0_s_corr = inp.ts0_s_corr;
  ret.ts0_ns = inp.ts0_ns;
  ret.ts0_ns_corr = inp.ts0_ns_corr; 
  ret.ts1_ns = inp.ts1_ns;
  ret.plane = inp.plane;
  ret.x_pos = inp.x_pos;
  ret.x_err = inp.x_err;
  ret.y_pos = inp.y_pos;
  ret.y_err = inp.y_err;
  ret.z_pos = inp.z_pos;
  ret.z_err = inp.z_err;
  ret.tagger = inp.tagger;

  return ret;
} 

sbn::CRTHitMatch::CRTHitMatch(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fMatchAlg(fhicl::Table<sbnd::CRTT0MatchAlg::Config>(p.get<fhicl::ParameterSet>("Alg"))()),
    fCRTHitLabel(p.get<std::string>("CRTHitLabel", "crthit")),
    fTrackLabel(p.get<std::string>("TrackLabel", "pandoraTrack")),
    fTSMode(p.get<int>("Alg.TSMode", 1)),
    fTimeCorrection(p.get<float>("Alg.TimeCorrection", 0.)),
    fMinTrackLength(p.get<float>("MinTrackLength", 10.))
{
  produces< art::Assns<recob::Track, sbn::crt::CRTHit, anab::T0> >();
}

void sbn::CRTHitMatch::produce(art::Event& e)
{
  std::unique_ptr<art::Assns<recob::Track, sbn::crt::CRTHit, anab::T0>> assn(new art::Assns<recob::Track, sbn::crt::CRTHit, anab::T0>);

  art::Handle<std::vector<recob::Track>> track_handle;
  e.getByLabel(fTrackLabel, track_handle);

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, track_handle);

  art::FindManyP<recob::Hit> fmHits(tracks, e, fTrackLabel); 
  // TODO: use hits directly in CRT hit matching
  (void) fmHits;

  art::Handle<std::vector<sbn::crt::CRTHit>> crthit_handle;
  e.getByLabel(fCRTHitLabel, crthit_handle);

  std::vector<art::Ptr<sbn::crt::CRTHit>> crthits;
  art::fill_ptr_vector(crthits, crthit_handle);
  
  std::vector<sbnd::crt::CRTHit> sbnd_crthits;
  for (art::Ptr<sbn::crt::CRTHit> hit: crthits) {
     sbnd_crthits.push_back(SBN2SBNDCrtHit(*hit));
  }

  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  for (unsigned i = 0; i < tracks.size(); i++) {
    if (fMinTrackLength > 0. && tracks[i]->Length() < fMinTrackLength) continue;

    std::pair<sbnd::crt::CRTHit, double> hit_pair = fMatchAlg.ClosestCRTHit(det_prop, *tracks[i], sbnd_crthits, e);
    if (hit_pair.second >= 0) {
      // TODO: fix hacky BS
      // figure out which hit was matched
      int match_ind = -1;
      for (unsigned j = 0; j < crthits.size(); j++) {
        if (hit_pair.first.ts0_ns == crthits[j]->ts0_ns &&
            hit_pair.first.ts1_ns == crthits[j]->ts1_ns &&
            hit_pair.first.x_pos == crthits[j]->x_pos &&
            hit_pair.first.y_pos == crthits[j]->y_pos &&
            hit_pair.first.z_pos == crthits[j]->z_pos) {
          match_ind = j;
          break;
        }
      }
      anab::T0 t0;
      if (fTSMode == 1) {
        t0.fTime = ((double)(int)crthits[match_ind]->ts1_ns) * 1e-3 + fTimeCorrection;
      }
      else {
        t0.fTime = ((double)(int)crthits[match_ind]->ts0_ns) * 1e-3 + fTimeCorrection;
      }
      t0.fTriggerConfidence = hit_pair.second;
      assn->addSingle(tracks[i], crthits[match_ind], t0);
    }
  }
  e.put(std::move(assn));
}

DEFINE_ART_MODULE(sbn::CRTHitMatch)
