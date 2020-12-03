#ifndef CAF_CAFMAKERPARAMS_H
#define CAF_CAFMAKERPARAMS_H

#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDAnalyzer.h"

namespace caf
{
  struct CAFMakerParams
  {
    template<class T> using Atom = fhicl::Atom<T>;
    template<class T> using Table = fhicl::Table<T>;
    using Comment  = fhicl::Comment;
    using Name     = fhicl::Name;
    using string   = std::string;
    using InputTag = art::InputTag;

    /* Atom<bool> EnableBlindness */
    /* { */
    /*   Name("EnableBlindness"), */
    /*   Comment("true = hide sensitive info, false = include full record") */
    /* }; */

    Atom<std::string> CAFFilename { Name("CAFFilename"),
      Comment("Provide a string to override the automatic filename.")
    };

    Atom<string> DataTier        { Name("DataTier") };
    Atom<string> FileExtension   { Name("FileExtension") };
    Atom<string> GeneratorLabel  { Name("GeneratorInput") };

    Atom<bool> StrictMode        { Name("StrictMode"),
	Comment("Abort if any required product not found, unless label is empty")
    };

    Atom<bool> CutClearCosmic {
      Name("CutClearCosmic"),
      Comment("Cut slices which are marked as a 'clear-cosmic' by pandora"),
      false
    };

    Atom<bool> SelectOneSlice {
      Name("SelectOneSlice"),
      Comment("Only select one slice per spill (ranked by nu_score) [TODO: implement]."),
      false
    };

    fhicl::OptionalSequence<std::string> PandoraTagSuffixes {
      Name("PandoraTagSuffixes"),
      Comment("List of suffixes to add to TPC reco tag names (e.g. cryo0 cryo1)")
    };

    Atom<string> G4Label {
      Name("G4Label"),
      Comment("Label of G4 module."),
      "largeant"
    };

    Atom<string> GenLabel {
      Name("GenLabel"),
      Comment("Label of neutrino gen module."),
      "generator"
    };

    Atom<string> CosmicGenLabel {
      Name("CosmicGenLabel"),
      Comment("Label of cosmic gen module."),
      "cosmgen"
    };

    Atom<string> ParticleGunGenLabel {
      Name("ParticleGunGenLabel"),
      Comment("Label of particle gun gen module."),
      "particlegun"
    };

    Atom<string> PFParticleLabel {
      Name("PFParticleLabel"),
      Comment("Base label of PFParticle producer."),
      "pandora"
    };

    Atom<string> FlashMatchLabel {
      Name("FlashMatchLabel"),
      Comment("Base label of flash match producer."),
      "fmatch" // same for icarus and sbnd
    };

    Atom<string> RecoTrackLabel {
      Name("RecoTrackLabel"),
      Comment("Base label of reco-base track producer."),
      "pandoraTrack"
    };

    Atom<string> RecoShowerLabel {
      Name("RecoShowerLabel"),
      Comment("Base label of reco-base shower producer."),
      "pandoraShowerSBN"
    };

    Atom<string> RecoShowerSelectionLabel {
      Name("RecoShowerSelectionLabel"),
      Comment("Base label of shower selection vars producer."),
      "pandoraShowerSelectionVars"
    };

    Atom<string> TrackCaloLabel {
      Name("TrackCaloLabel"),
      Comment("Base label of track calorimetry producer."),
      "pandoraCalo"
    };

    Atom<string> TrackPidLabel {
      Name("TrackPidLabel"),
      Comment("Base label of track particle-id producer."),
      "pandoraPid"
    };

    Atom<string> CRTHitMatchLabel {
      Name("CRTHitMatchLabel"),
      Comment("Base label of track to CRT hit matching producer."),
      "pandoraTrackCRTHit"
    };

    Atom<string> CRTTrackMatchLabel {
      Name("CRTTrackMatchLabel"),
      Comment("Base label of track to CRT track matching producer."),
      "pandoraTrackCRTTrack"
    };

    Atom<string> TrackMCSLabel {
      Name("TrackMCSLabel"),
      Comment("Base label of track MCS momentum calculation producer."),
      "pandoraTrackMCS"
    };

    Atom<string> TrackRangeLabel {
      Name("TrackRangeLabel"),
      Comment("Base label of track range momentum calculation producer."),
      "pandoraTrackRange"
    };

    Atom<string> CRTHitLabel {
      Name("CRTHitLabel"),
      Comment("Label of sbn CRT hits."),
      "crthit" // same for icarus and sbnd
    };

    Atom<string> CRTTrackLabel {
      Name("CRTTrackLabel"),
      Comment("Label of sbn CRT tracks."),
      "crttrack" // same for icarus and sbnd
    };
    
    Atom<string> FlashTrigLabel {
      Name("FlashTrigLabel"),
      Comment("Label of bool of passing flash trigger."),
      "flashtrigfilter"
    };

    Atom<bool> CRTUseTS0 {
      Name("CRTUseTS0"),
      Comment("Whether to use ts0 or ts1 to fill the time of the SRCRTHit and SRCRTTrack"),
      false
    };

    Atom<string> SimChannelLabel {
      Name("SimChannelLabel"),
      Comment("Label of input sim::SimChannel objects."),
      "largeant"
    };

    fhicl::Sequence<float, 3u> CalorimetryConstants {
      Name("CalorimetryConstants"),
      Comment("Constants to convert ADC*tick charge measurement to electrons."
              "Ordered 1st Induction, 2nd Induction, Collection."
              "In units of ADC*tick / electrons")
    };

  };
}

#endif
