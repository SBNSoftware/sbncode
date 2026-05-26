////////////////////////////////////////////////////////////////////////
// Class:       BadChannelHitCheck
// Plugin Type: analyzer
// File:        BadChannelHitCheck_module.cc
//
// Validation analyzer for the GaussHitFinderSBN "ExcludeBadChannels"
// feature.  For each configured recob::Hit collection it checks every
// hit's channel against the channel-status database service; if any hit
// sits on a channel the database flags as bad, the art job is failed
// (by throwing).  This lets a ctest verify end-to-end that bad channels
// are excluded from the hit finder output.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include <vector>

namespace sbn {

  class BadChannelHitCheck : public art::EDAnalyzer {
  public:
    struct Config {
      fhicl::Sequence<art::InputTag> HitLabels{
        fhicl::Name("HitLabels"),
        fhicl::Comment("recob::Hit collections to validate against the bad-channel database")};
      fhicl::Atom<bool> RequirePresent{
        fhicl::Name("RequirePresent"),
        fhicl::Comment("Also fail when a hit lands on a channel the database does not list as present"),
        false};
      fhicl::Atom<bool> ThrowOnFailure{
        fhicl::Name("ThrowOnFailure"),
        fhicl::Comment("Throw on the first bad-channel hit; if false, collect and throw in endJob()"),
        true};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit BadChannelHitCheck(Parameters const& p);

    void analyze(art::Event const& e) override;
    void endJob() override;

  private:
    std::vector<art::InputTag> const fHitLabels;
    bool const fRequirePresent;
    bool const fThrowOnFailure;

    unsigned long fNHits{0};      ///< total hits checked
    unsigned long fNBad{0};       ///< hits found on bad/not-present channels
    bool fReportedStatus{false};  ///< whether the DB bad-channel count was logged
  };

  //-------------------------------------------------------------------
  BadChannelHitCheck::BadChannelHitCheck(Parameters const& p)
    : art::EDAnalyzer{p}
    , fHitLabels{p().HitLabels()}
    , fRequirePresent{p().RequirePresent()}
    , fThrowOnFailure{p().ThrowOnFailure()}
  {}

  //-------------------------------------------------------------------
  void BadChannelHitCheck::analyze(art::Event const& e)
  {
    lariov::ChannelStatusProvider const& chanStatus =
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    if (!fReportedStatus) {
      fReportedStatus = true;
      mf::LogInfo("BadChannelHitCheck")
        << "channel-status database reports " << chanStatus.BadChannels().size()
        << " bad and " << chanStatus.NoisyChannels().size() << " noisy channel(s).";
    }

    for (art::InputTag const& tag : fHitLabels) {
      auto const& hits = e.getProduct<std::vector<recob::Hit>>(tag);

      for (recob::Hit const& hit : hits) {
        raw::ChannelID_t const channel = hit.Channel();
        ++fNHits;

        bool const isBad = chanStatus.IsBad(channel);
        bool const isAbsent = fRequirePresent && !chanStatus.IsPresent(channel);
        if (!isBad && !isAbsent) continue;

        ++fNBad;
        mf::LogError("BadChannelHitCheck")
          << "recob::Hit on " << (isBad ? "BAD" : "NOT-PRESENT") << " channel " << channel
          << " in product '" << tag.encode() << "', event " << e.id();

        if (fThrowOnFailure) {
          throw art::Exception(art::errors::LogicError)
            << "BadChannelHitCheck: recob::Hit found on "
            << (isBad ? "bad" : "not-present") << " channel " << channel
            << " (product '" << tag.encode() << "', event " << e.id() << ")\n";
        }
      }
    }
  }

  //-------------------------------------------------------------------
  void BadChannelHitCheck::endJob()
  {
    mf::LogInfo("BadChannelHitCheck")
      << "Checked " << fNHits << " hit(s); " << fNBad << " on bad/not-present channels.";

    if (fNBad > 0) {
      throw art::Exception(art::errors::LogicError)
        << "BadChannelHitCheck FAILED: " << fNBad << " of " << fNHits
        << " hit(s) are on bad/not-present channels.\n";
    }
  }

} // namespace sbn

DEFINE_ART_MODULE(sbn::BadChannelHitCheck)
