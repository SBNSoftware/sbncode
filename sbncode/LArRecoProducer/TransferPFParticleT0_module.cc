////////////////////////////////////////////////////////////////////////
// Class:       TransferPFParticleT0
// Plugin Type: producer (art v3_02_06)
// File:        TransferPFParticleT0_module.cc
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
#include "lardataobj/AnalysisBase/T0.h"

#include <memory>

namespace sbn {
  class TransferPFParticleT0;
}


class sbn::TransferPFParticleT0 : public art::EDProducer {
public:
  explicit TransferPFParticleT0(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TransferPFParticleT0(TransferPFParticleT0 const&) = delete;
  TransferPFParticleT0(TransferPFParticleT0&&) = delete;
  TransferPFParticleT0& operator=(TransferPFParticleT0 const&) = delete;
  TransferPFParticleT0& operator=(TransferPFParticleT0&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fPFParticleLabelIn;
  art::InputTag fT0Label;

  art::InputTag fPFParticleLabelOut;
};


sbn::TransferPFParticleT0::TransferPFParticleT0(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabelIn(p.get<art::InputTag>("PFParticleLabelIn")),
    fT0Label(p.get<art::InputTag>("T0Label")),
    fPFParticleLabelOut(p.get<art::InputTag>("PFParticleLabelOut"))
{
  produces<art::Assns<recob::PFParticle, anab::T0>>();
}

void sbn::TransferPFParticleT0::produce(art::Event& e)
{
  // define output data
  std::unique_ptr<art::Assns<recob::PFParticle, anab::T0>> assn(new art::Assns<recob::PFParticle, anab::T0>);

  // Collect input data
  art::Handle<std::vector<recob::PFParticle>> inpfp_handle;
  e.getByLabel(fPFParticleLabelIn, inpfp_handle);

  std::vector<art::Ptr<recob::PFParticle>> inpfps;
  art::fill_ptr_vector(inpfps, inpfp_handle);

  art::FindManyP<anab::T0> fmT0s(inpfps, e, fT0Label);

  art::Handle<std::vector<recob::PFParticle>> outpfp_handle;
  e.getByLabel(fPFParticleLabelOut, outpfp_handle);

  std::vector<art::Ptr<recob::PFParticle>> outpfps;
  art::fill_ptr_vector(outpfps, outpfp_handle);

  // process input data
  for (unsigned i_in = 0; i_in < inpfps.size(); i_in++) {
    const recob::PFParticle &pfp = *inpfps[i_in];
    const std::vector<art::Ptr<anab::T0>> &t0s = fmT0s.at(i_in);
    if (!t0s.size()) continue; // ignore pfps without T0's 

    for (unsigned i_out = 0; i_out < outpfps.size(); i_out++) {
      const recob::PFParticle &out_pfp = *outpfps[i_out];
      if (pfp.Self() == out_pfp.Self()) { // The pfps match!
        // Transfer over the T0 association
        assn->addMany(outpfps.at(i_out), t0s);
        break; // done with this input pfp
      }
    }
  }

  e.put(std::move(assn));
}

DEFINE_ART_MODULE(sbn::TransferPFParticleT0)
