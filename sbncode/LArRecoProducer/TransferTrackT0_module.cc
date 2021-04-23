////////////////////////////////////////////////////////////////////////
// Class:       TransferTrackT0
// Plugin Type: producer (art v3_02_06)
// File:        TransferTrackT0_module.cc
//
// Generated at Wed Feb 19 17:38:15 2020 by Gray Putnam using cetskelgen
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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <memory>

namespace sbn {
  class TransferTrackT0;
}


class sbn::TransferTrackT0 : public art::EDProducer {
public:
  explicit TransferTrackT0(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TransferTrackT0(TransferTrackT0 const&) = delete;
  TransferTrackT0(TransferTrackT0&&) = delete;
  TransferTrackT0& operator=(TransferTrackT0 const&) = delete;
  TransferTrackT0& operator=(TransferTrackT0&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fTrackLabelIn;
  art::InputTag fT0Label;

  art::InputTag fTrackLabelOut;
};


sbn::TransferTrackT0::TransferTrackT0(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fTrackLabelIn(p.get<art::InputTag>("TrackLabelIn")),
    fT0Label(p.get<art::InputTag>("T0Label")),
    fTrackLabelOut(p.get<art::InputTag>("TrackLabelOut"))
{
  produces<art::Assns<recob::Track, anab::T0>>();
}

void sbn::TransferTrackT0::produce(art::Event& e)
{
  // define output data
  std::unique_ptr<art::Assns<recob::Track, anab::T0>> assn(new art::Assns<recob::Track, anab::T0>);

  // Collect input data
  art::Handle<std::vector<recob::Track>> intrack_handle;
  e.getByLabel(fTrackLabelIn, intrack_handle);

  std::vector<art::Ptr<recob::Track>> intracks;
  art::fill_ptr_vector(intracks, intrack_handle);

  art::FindManyP<anab::T0> fmT0s(intracks, e, fT0Label);
  art::FindManyP<recob::PFParticle> fmInParticles(intracks, e, fTrackLabelIn);

  art::Handle<std::vector<recob::Track>> outtrack_handle;
  e.getByLabel(fTrackLabelOut, outtrack_handle);

  std::vector<art::Ptr<recob::Track>> outtracks;
  art::fill_ptr_vector(outtracks, outtrack_handle);

  art::FindManyP<recob::PFParticle> fmOutParticles(outtracks, e, fTrackLabelOut);

  // process input data
  for (unsigned i_in = 0; i_in < intracks.size(); i_in++) {
    const recob::PFParticle &pfp = *fmInParticles.at(i_in).at(0);
    const std::vector<art::Ptr<anab::T0>> &t0s = fmT0s.at(i_in);
    if (!t0s.size()) continue; // ignore tracks without T0's 

    for (unsigned i_out = 0; i_out < outtracks.size(); i_out++) {
      const recob::PFParticle &out_pfp = *fmOutParticles.at(i_out).at(0);
      if (pfp.Self() == out_pfp.Self()) { // The tracks match!
        // Transfer over the T0 association
        assn->addMany(outtracks.at(i_out), t0s);
        break; // done with this input track
      }
    }
  }

  e.put(std::move(assn));
}

DEFINE_ART_MODULE(sbn::TransferTrackT0)
