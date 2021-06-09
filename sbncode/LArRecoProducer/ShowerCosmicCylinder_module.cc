////////////////////////////////////////////////////////////////////////
// Class:       ShowerCosmicCylinder
// Plugin Type: producer (art v3_03_01)
// File:        ShowerCosmicCylinder_module.cc
//
// Generated at Tue Sep 15 08:54:39 2020 by Edward Tyley using cetskelgen
// from cetlib version v3_08_00.
//
// Producer to find the distance of closest approach from a shower to
// a unambiguous cosmic ray
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

#include <memory>

#include "TF1.h"
#include "TGraph.h"

namespace sbn {
class ShowerCosmicCylinder;
}

class sbn::ShowerCosmicCylinder : public art::EDProducer {
  public:
  explicit ShowerCosmicCylinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerCosmicCylinder(ShowerCosmicCylinder const&) = delete;
  ShowerCosmicCylinder(ShowerCosmicCylinder&&) = delete;
  ShowerCosmicCylinder& operator=(ShowerCosmicCylinder const&) = delete;
  ShowerCosmicCylinder& operator=(ShowerCosmicCylinder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  private:
  // FCL params
  const art::InputTag fPandoraLabel;
  const art::InputTag fShowerLabel;

  const float fMinShowerEnergy;
};

sbn::ShowerCosmicCylinder::ShowerCosmicCylinder(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fPandoraLabel(p.get<art::InputTag>("PandoraLabel"))
    , fShowerLabel(p.get<art::InputTag>("ShowerLabel"))
    , fMinShowerEnergy(p.get<float>("MinShowerEnergy"))
{
  produces<std::vector<float>>();
  produces<art::Assns<recob::Shower, float>>();
}

void sbn::ShowerCosmicCylinder::produce(art::Event& e)
{
  //Get the showers
  auto const showerHandle = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);
  auto const sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fPandoraLabel);
  auto const pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
  // auto const spHandle = e.getValidHandle<std::vector<recob::SpacePoint> >(fPandoraLabel);

  std::vector<art::Ptr<recob::Shower>> showers;
  art::fill_ptr_vector(showers, showerHandle);

  std::vector<art::Ptr<recob::Slice>> slices;
  art::fill_ptr_vector(slices, sliceHandle);

  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  const art::FindManyP<recob::PFParticle> fmSlicePFP(sliceHandle, e, fPandoraLabel);
  if (!fmSlicePFP.isValid()) {
    throw cet::exception("ShowerCosmicCylinder") << "Slice-PFP association is somehow not valid. Stopping";
    return;
  }
  const art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta(pfpHandle, e, fPandoraLabel);
  if (!fmPFPMeta.isValid()) {
    throw cet::exception("ShowerCosmicCylinder") << "PFP-Meta association is somehow not valid. Stopping";
    return;
  }
  const art::FindManyP<recob::SpacePoint> fmPFPSP(pfpHandle, e, fPandoraLabel);
  if (!fmPFPSP.isValid()) {
    throw cet::exception("ShowerCosmicCylinder") << "PFP-SP association is somehow not valid. Stopping";
    return;
  }

  // Find the primary cosmic PFPs
  std::vector<art::Ptr<recob::PFParticle>> cosmicPFPs;
  for (auto const& slice : slices) {
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(fmSlicePFP.at(slice.key()));

    if (slicePFPs.empty())
      continue;

    // Find the Primary PFP
    std::vector<art::Ptr<recob::PFParticle>> pfpPrimaries;
    for (auto const& pfp : slicePFPs) {
      if (pfp->IsPrimary())
        pfpPrimaries.push_back(pfp);
    }

    // We should always have 1 primary
    if (pfpPrimaries.size() != 1)
      throw cet::exception("ShowerCosmicCylinder") << "Wrong number of primaries in slice: " << pfpPrimaries.size() << ". Stopping";

    // Ignore anything that is not a clear cosmic
    // TODO: Should expand to check for T0 tags from CRT/Flash matching
    auto const& pfpMetaVec(fmPFPMeta.at(pfpPrimaries.front().key()));
    if (pfpMetaVec.size() != 1)
      throw cet::exception("ShowerCosmicCylinder") << "Wrong metadata entries for PFP: " << pfpMetaVec.size() << ". Stopping";
    if (!pfpMetaVec.front()->GetPropertiesMap().count("IsClearCosmic"))
      continue;

    cosmicPFPs.push_back(pfpPrimaries.front());
  } // end pfp: slicePFPs

  std::unique_ptr<std::vector<float>> residualCol(std::make_unique<std::vector<float>>());
  std::unique_ptr<art::Assns<recob::Shower, float>> residualAssns(std::make_unique<art::Assns<recob::Shower, float>>());

  for (auto const& shower : showers) {

    // To reduce the combinatorics, put a minumum energy cut on the showers we consider
    if (shower->best_plane() < 0 || shower->Energy().at(shower->best_plane()) < fMinShowerEnergy)
      continue;

    float res(std::numeric_limits<float>::max());
    const TVector3 showerStart3(shower->ShowerStart());
    const geo::Point_t showerStart(showerStart3.X(), showerStart3.Y(), showerStart3.Z());

    // Loop over the cosmic PFPs and find the closest
    for (auto const& cosmicPFP : cosmicPFPs) {
      auto const& pfpSPs(fmPFPSP.at(cosmicPFP.key()));
      // Loop over all of the spacepoints in the PFP and find the closest
      for (auto const& sp : pfpSPs) {
        const float dist((showerStart - sp->position()).r());
        res = std::min(res, dist);
      }
    }

    residualCol->push_back(res);
    util::CreateAssn(*this, e, *residualCol, shower, *residualAssns);
  }

  e.put(std::move(residualCol));
  e.put(std::move(residualAssns));
}

DEFINE_ART_MODULE(sbn::ShowerCosmicCylinder)
