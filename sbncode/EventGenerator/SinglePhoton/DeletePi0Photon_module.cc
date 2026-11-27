////////////////////////////////////////////////////////////////////////
/// \file   DeletePi0Photon_module.cc
/// \brief  Builds single-photon training samples from GENIE events by
///         decaying pi0s at generator level and deleting one photon.
///
/// Intended to run in the gen stage directly after GENIEGen, which is
/// relabelled (e.g. "rawgenerator") so that this module can take the
/// "generator" label. It re-publishes every product GENIEGen makes,
/// so downstream stages (G4, detsim, reco, CAFMaker, event weights) run
/// unchanged on the modified truth:
///
///   event  : std::vector<simb::MCTruth>        (modified)
///            std::vector<simb::MCFlux>         (copied, index-aligned with MCTruth)
///            std::vector<simb::GTruth>         (copied, index-aligned with MCTruth)
///            art::Assns<MCTruth, MCFlux/GTruth> (rebuilt)
///            std::vector<sim::BeamGateInfo>    (copied)
///            Dk2Nu / NuChoice + Assns          (copied, if ForwardDk2Nu)
///   subrun : sumdata::POTSummary               (copied)
///   run    : sumdata::RunData                  (rebuilt from Geometry, as GENIEGen does)
///
/// The input label's products MUST be dropped on output, otherwise
/// larg4Main (which reads every MCTruth collection when inputCollections
/// is empty) simulates both the original and the modified interaction.
///
/// Note that GTruth, MCFlux and the MCNeutrino describe the *unmodified*
/// interaction; the modified events do not conserve 4-momentum.
///
/// As an EDFilter it can also require an exact number of final-state
/// photons after modification (RequiredNPhotons; < 0 disables). With
/// OnlyDeletePi0Photons, both the deletion and this count consider only
/// photons whose mother is a pi0, so other photons in the event (e.g.
/// GENIE nuclear de-excitation gammas) are left alone and not counted.
////////////////////////////////////////////////////////////////////////

#include "sbncode/EventGenerator/SinglePhoton/PhotonTruthTools.h"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dk2nu/tree/NuChoice.h"
#include "dk2nu/tree/dk2nu.h"

#include "CLHEP/Random/RandomEngine.h"

#include <memory>
#include <vector>

namespace sbn {

  class DeletePi0Photon : public art::EDFilter {
  public:
    explicit DeletePi0Photon(fhicl::ParameterSet const& p);

    bool filter(art::Event& e) override;
    bool beginRun(art::Run& r) override;
    bool endSubRun(art::SubRun& sr) override;
    void endJob() override;

  private:
    art::InputTag const fInputTag;
    bool const fDecayPi0s;
    bool const fDeletePhoton;
    bool const fOnlyDeletePi0Photons;
    int const fRequiredNPhotons;
    bool const fForwardDk2Nu;
    bool const fVerbose;

    CLHEP::HepRandomEngine& fEngine;

    // counters for the end-of-job summary
    unsigned long fNEvents = 0;
    unsigned long fNPassed = 0;
    unsigned long fNPi0sDecayed = 0;
    unsigned long fNPhotonsDeleted = 0;
  };

} // namespace sbn

//------------------------------------------------------------------------------
sbn::DeletePi0Photon::DeletePi0Photon(fhicl::ParameterSet const& p)
  : EDFilter{p}
  , fInputTag{p.get<art::InputTag>("InputLabel")}
  , fDecayPi0s{p.get<bool>("DecayPi0s", true)}
  , fDeletePhoton{p.get<bool>("DeletePhoton", true)}
  , fOnlyDeletePi0Photons{p.get<bool>("OnlyDeletePi0Photons", false)}
  , fRequiredNPhotons{p.get<int>("RequiredNPhotons", -1)}
  , fForwardDk2Nu{p.get<bool>("ForwardDk2Nu", false)}
  , fVerbose{p.get<bool>("Verbose", false)}
  , fEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
      createEngine(0, "HepJamesRandom", "deletepi0photon"),
      "HepJamesRandom",
      "deletepi0photon",
      p,
      "Seed"))
{
  produces<std::vector<simb::MCTruth>>();
  produces<std::vector<simb::MCFlux>>();
  produces<std::vector<simb::GTruth>>();
  produces<art::Assns<simb::MCTruth, simb::MCFlux>>();
  produces<art::Assns<simb::MCTruth, simb::GTruth>>();
  produces<std::vector<sim::BeamGateInfo>>();
  produces<sumdata::POTSummary, art::InSubRun>();
  produces<sumdata::RunData, art::InRun>();

  if (fForwardDk2Nu) {
    produces<std::vector<bsim::Dk2Nu>>();
    produces<std::vector<bsim::NuChoice>>();
    produces<art::Assns<simb::MCTruth, bsim::Dk2Nu>>();
    produces<art::Assns<simb::MCTruth, bsim::NuChoice>>();
  }
}

//------------------------------------------------------------------------------
namespace {

  /// Copy the T associated (one-to-one) with each input MCTruth into a new
  /// collection, in MCTruth order, and associate it with the new MCTruth.
  /// CAFMaker indexes MCFlux and Dk2Nu by MCTruth position, so the order matters.
  template <typename T, typename H>
  void CopyAssociated(art::Event& e,
                      H const& inTruthH,
                      art::InputTag const& tag,
                      art::PtrMaker<simb::MCTruth> const& truthPtr,
                      std::vector<T>& outCol,
                      art::Assns<simb::MCTruth, T>& outAssn,
                      char const* what)
  {
    art::FindOneP<T> const fop(inTruthH, e, tag);
    art::PtrMaker<T> const makePtr{e};
    for (std::size_t i = 0; i < inTruthH->size(); ++i) {
      art::Ptr<T> const& in = fop.at(i);
      if (in.isNull())
        throw cet::exception("DeletePi0Photon")
          << "MCTruth #" << i << " from '" << tag.encode() << "' has no associated " << what
          << ".\n";
      outCol.push_back(*in);
      outAssn.addSingle(truthPtr(i), makePtr(outCol.size() - 1));
    }
  }

} // namespace

//------------------------------------------------------------------------------
bool sbn::DeletePi0Photon::filter(art::Event& e)
{
  ++fNEvents;

  auto const inTruthH = e.getValidHandle<std::vector<simb::MCTruth>>(fInputTag);

  auto outTruth = std::make_unique<std::vector<simb::MCTruth>>();
  outTruth->reserve(inTruthH->size());

  unsigned int nPhotons = 0;
  for (std::size_t i = 0; i < inTruthH->size(); ++i) {
    simb::MCTruth truth = (*inTruthH)[i];

    if (fDecayPi0s) {
      simb::MCTruth decayed;
      fNPi0sDecayed += singlephoton::DecayPi0sToTwoPhotons(truth, decayed, fEngine);
      truth = std::move(decayed);
    }

    if (fDeletePhoton) {
      simb::MCTruth deleted;
      int const idx =
        singlephoton::DeleteOneRandomPhoton(truth, deleted, fEngine, fOnlyDeletePi0Photons);
      if (idx >= 0) {
        ++fNPhotonsDeleted;
        if (fVerbose) {
          auto const& g = truth.GetParticle(idx);
          mf::LogInfo("DeletePi0Photon")
            << "Event " << e.id() << " MCTruth #" << i << ": deleted photon track ID "
            << g.TrackId() << " (mother " << g.Mother() << ", E = " << g.E() << " GeV)";
        }
      }
      truth = std::move(deleted);
    }

    nPhotons += singlephoton::CountFinalStatePhotons(truth, fOnlyDeletePi0Photons);
    outTruth->push_back(std::move(truth));
  }

  // Everything else is copied unchanged, index-aligned with the MCTruth.
  art::PtrMaker<simb::MCTruth> const truthPtr{e};

  auto outFlux = std::make_unique<std::vector<simb::MCFlux>>();
  auto outFluxAssn = std::make_unique<art::Assns<simb::MCTruth, simb::MCFlux>>();
  CopyAssociated(e, inTruthH, fInputTag, truthPtr, *outFlux, *outFluxAssn, "MCFlux");

  auto outGTruth = std::make_unique<std::vector<simb::GTruth>>();
  auto outGTruthAssn = std::make_unique<art::Assns<simb::MCTruth, simb::GTruth>>();
  CopyAssociated(e, inTruthH, fInputTag, truthPtr, *outGTruth, *outGTruthAssn, "GTruth");

  auto outGates = std::make_unique<std::vector<sim::BeamGateInfo>>(
    *e.getValidHandle<std::vector<sim::BeamGateInfo>>(fInputTag));

  e.put(std::move(outTruth));
  e.put(std::move(outFlux));
  e.put(std::move(outGTruth));
  e.put(std::move(outFluxAssn));
  e.put(std::move(outGTruthAssn));
  e.put(std::move(outGates));

  // dk2nu fluxes (e.g. SBND's G4BNB) give Dk2Nu products that CAFMaker requires;
  // refuse to silently lose them.
  if (!fForwardDk2Nu && e.getHandle<std::vector<bsim::Dk2Nu>>(fInputTag))
    throw cet::exception("DeletePi0Photon")
      << "Input '" << fInputTag.encode() << "' has bsim::Dk2Nu products (dk2nu flux) but "
      << "ForwardDk2Nu is false; set ForwardDk2Nu: true.\n";

  if (fForwardDk2Nu) {
    auto outDk2Nu = std::make_unique<std::vector<bsim::Dk2Nu>>();
    auto outDk2NuAssn = std::make_unique<art::Assns<simb::MCTruth, bsim::Dk2Nu>>();
    CopyAssociated(e, inTruthH, fInputTag, truthPtr, *outDk2Nu, *outDk2NuAssn, "Dk2Nu");

    auto outNuChoice = std::make_unique<std::vector<bsim::NuChoice>>();
    auto outNuChoiceAssn = std::make_unique<art::Assns<simb::MCTruth, bsim::NuChoice>>();
    CopyAssociated(
      e, inTruthH, fInputTag, truthPtr, *outNuChoice, *outNuChoiceAssn, "NuChoice");

    e.put(std::move(outDk2Nu));
    e.put(std::move(outNuChoice));
    e.put(std::move(outDk2NuAssn));
    e.put(std::move(outNuChoiceAssn));
  }

  bool const pass =
    fRequiredNPhotons < 0 || nPhotons == static_cast<unsigned int>(fRequiredNPhotons);
  if (pass) ++fNPassed;
  return pass;
}

//------------------------------------------------------------------------------
bool sbn::DeletePi0Photon::beginRun(art::Run& r)
{
  // Same as GENIEGen::beginRun.
  art::ServiceHandle<geo::Geometry const> geo;
  r.put(std::make_unique<sumdata::RunData>(geo->DetectorName()), art::fullRun());
  return true;
}

//------------------------------------------------------------------------------
bool sbn::DeletePi0Photon::endSubRun(art::SubRun& sr)
{
  // GENIEGen puts its POTSummary in endSubRun; this module must therefore be
  // scheduled after it in the same path.
  auto const pot = sr.getHandle<sumdata::POTSummary>(fInputTag);
  if (!pot)
    throw cet::exception("DeletePi0Photon")
      << "No sumdata::POTSummary found under '" << fInputTag.encode()
      << "' in endSubRun. Is this module after the generator in the path?\n";
  sr.put(std::make_unique<sumdata::POTSummary>(*pot), art::subRunFragment());
  return true;
}

//------------------------------------------------------------------------------
void sbn::DeletePi0Photon::endJob()
{
  mf::LogInfo("DeletePi0Photon") << "Summary: " << fNEvents << " events read, "
                                 << fNPi0sDecayed << " pi0s decayed, " << fNPhotonsDeleted
                                 << " photons deleted, " << fNPassed
                                 << " events passed (RequiredNPhotons = " << fRequiredNPhotons
                                 << ").";
}

DEFINE_ART_MODULE(sbn::DeletePi0Photon)
