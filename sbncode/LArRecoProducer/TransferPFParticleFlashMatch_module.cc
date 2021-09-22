////////////////////////////////////////////////////////////////////////
// Class:       TransferPFParticleFlashMatch
// Plugin Type: producer (art v3_02_06)
// File:        TransferPFParticleFlashMatch.cc
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
#include "lardataobj/RecoBase/PFParticle.h"

#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"

#include <memory>

namespace sbn {
  class TransferPFParticleFlashMatch;
}


class sbn::TransferPFParticleFlashMatch : public art::EDProducer {
public:
  explicit TransferPFParticleFlashMatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TransferPFParticleFlashMatch(TransferPFParticleFlashMatch const&) = delete;
  TransferPFParticleFlashMatch(TransferPFParticleFlashMatch&&) = delete;
  TransferPFParticleFlashMatch& operator=(TransferPFParticleFlashMatch const&) = delete;
  TransferPFParticleFlashMatch& operator=(TransferPFParticleFlashMatch&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fPFParticleLabelIn;
  art::InputTag fFlashLabel;

  art::InputTag fPFParticleLabelOut;
};


sbn::TransferPFParticleFlashMatch::TransferPFParticleFlashMatch(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabelIn(p.get<art::InputTag>("PFParticleLabelIn")),
    fFlashLabel(p.get<art::InputTag>("FlashLabel")),
    fPFParticleLabelOut(p.get<art::InputTag>("PFParticleLabelOut"))
{
  produces<art::Assns<recob::PFParticle, sbn::SimpleFlashMatch>>();
}

void sbn::TransferPFParticleFlashMatch::produce(art::Event& e)
{
  // define output data
  std::unique_ptr<art::Assns<recob::PFParticle, sbn::SimpleFlashMatch>> assn(new art::Assns<recob::PFParticle, sbn::SimpleFlashMatch>);

  // Collect input data
  art::Handle<std::vector<recob::PFParticle>> inpfp_handle;
  e.getByLabel(fPFParticleLabelIn, inpfp_handle);

  std::vector<art::Ptr<recob::PFParticle>> inpfps;
  art::fill_ptr_vector(inpfps, inpfp_handle);

  art::FindManyP<sbn::SimpleFlashMatch> fmFlash(inpfps, e, fFlashLabel);

  art::Handle<std::vector<recob::PFParticle>> outpfp_handle;
  e.getByLabel(fPFParticleLabelOut, outpfp_handle);

  std::vector<art::Ptr<recob::PFParticle>> outpfps;
  art::fill_ptr_vector(outpfps, outpfp_handle);

  // process input data
  for (unsigned i_in = 0; i_in < inpfps.size(); i_in++) {
    const recob::PFParticle &pfp = *inpfps[i_in];
    const std::vector<art::Ptr<sbn::SimpleFlashMatch>> &flashVec = fmFlash.at(i_in);
    if (!flashVec.size()) continue; // ignore pfps without flashes 

    for (unsigned i_out = 0; i_out < outpfps.size(); i_out++) {
      const recob::PFParticle &out_pfp = *outpfps[i_out];
      if (pfp.Self() == out_pfp.Self()) { // The pfps match!
        // Transfer over the flash association
        assn->addMany(outpfps.at(i_out), flashVec);
        break; // done with this input pfp
      }
    }
  }

  e.put(std::move(assn));
}

DEFINE_ART_MODULE(sbn::TransferPFParticleFlashMatch)
