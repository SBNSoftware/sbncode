////////////////////////////////////////////////////////////////////////
// Class:       ShowerCosmicDistance
// Plugin Type: producer (art v3_03_01)
// File:        ShowerCosmicDistance_module.cc
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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

#include <memory>

#include "TF1.h"
#include "TGraph.h"

namespace sbn {
class ShowerCosmicDistance : public art::EDProducer {
  public:
  explicit ShowerCosmicDistance(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerCosmicDistance(ShowerCosmicDistance const&) = delete;
  ShowerCosmicDistance(ShowerCosmicDistance&&) = delete;
  ShowerCosmicDistance& operator=(ShowerCosmicDistance const&) = delete;
  ShowerCosmicDistance& operator=(ShowerCosmicDistance&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  private:
  // FCL params
  const art::InputTag fPandoraLabel;
  const art::InputTag fShowerLabel;

  const float fMinShowerEnergy;

  const std::vector<art::Ptr<recob::PFParticle>> GetCosmicPFPs(const std::vector<art::Ptr<recob::PFParticle>>& pfps,
      const art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta) const;

  const float FindShowerResidual(const recob::Shower& shower,
      const std::vector<art::Ptr<recob::PFParticle>>& cosmicPFPs,
      const art::FindManyP<recob::SpacePoint>& fmPFPSP) const;
};

ShowerCosmicDistance::ShowerCosmicDistance(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fPandoraLabel(p.get<art::InputTag>("PandoraLabel"))
    , fShowerLabel(p.get<art::InputTag>("ShowerLabel"))
    , fMinShowerEnergy(p.get<float>("MinShowerEnergy"))
{
  produces<std::vector<float>>();
  produces<art::Assns<recob::Shower, float>>();
}

void ShowerCosmicDistance::produce(art::Event& e)
{
  //Get the showers
  auto const showerHandle = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);
  auto const pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);

  std::vector<art::Ptr<recob::Shower>> showers;
  art::fill_ptr_vector(showers, showerHandle);

  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  const art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta(pfpHandle, e, fPandoraLabel);
  if (!fmPFPMeta.isValid()) {
    throw cet::exception("ShowerCosmicDistance") << "PFP-Meta association is somehow not valid. Stopping";
    return;
  }
  const art::FindManyP<recob::SpacePoint> fmPFPSP(pfpHandle, e, fPandoraLabel);
  if (!fmPFPSP.isValid()) {
    throw cet::exception("ShowerCosmicDistance") << "PFP-SP association is somehow not valid. Stopping";
    return;
  }

  const std::vector<art::Ptr<recob::PFParticle>> cosmicPFPs(GetCosmicPFPs(pfps, fmPFPMeta));

  std::unique_ptr<std::vector<float>> residualCol(std::make_unique<std::vector<float>>());
  std::unique_ptr<art::Assns<recob::Shower, float>> residualAssns(std::make_unique<art::Assns<recob::Shower, float>>());

  for (auto const& shower : showers) {

    // To reduce the combinatorics, put a minumum energy cut on the showers we consider
    if (shower->best_plane() < 0 || shower->Energy().at(shower->best_plane()) < fMinShowerEnergy)
      continue;

    const float res(FindShowerResidual(*shower, cosmicPFPs, fmPFPSP));

    residualCol->push_back(res);
    util::CreateAssn(*this, e, *residualCol, shower, *residualAssns);
  }

  e.put(std::move(residualCol));
  e.put(std::move(residualAssns));
}

const std::vector<art::Ptr<recob::PFParticle>> ShowerCosmicDistance::GetCosmicPFPs(
    const std::vector<art::Ptr<recob::PFParticle>>& pfps,
    const art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta) const
{
  std::vector<art::Ptr<recob::PFParticle>> cosmicPFPs;
  for (auto const& pfp : pfps) {

    if (!pfp->IsPrimary())
      continue;

    auto const& pfpMetaVec(fmPFPMeta.at(pfp.key()));
    if (pfpMetaVec.size() != 1)
      throw cet::exception("ShowerCosmicDistance") << "Wrong metadata entries for PFP: " << pfpMetaVec.size() << ". Stopping";

    // Ignore anything that is not a clear cosmic
    if (!pfpMetaVec.front()->GetPropertiesMap().count("IsClearCosmic"))
      continue;

    // TODO: Should expand to check for T0 tags from CRT/Flash matching

    cosmicPFPs.push_back(pfp);
  }
  return cosmicPFPs;
}

const float ShowerCosmicDistance::FindShowerResidual(const recob::Shower& shower,
    const std::vector<art::Ptr<recob::PFParticle>>& cosmicPFPs,
    const art::FindManyP<recob::SpacePoint>& fmPFPSP) const
{
  float res(std::numeric_limits<float>::max());
  const TVector3 showerStart3(shower.ShowerStart());
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
  return res;
}
}

DEFINE_ART_MODULE(sbn::ShowerCosmicDistance)
