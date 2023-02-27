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
    template<class T> using Sequence = fhicl::Sequence<T>;
    template<class T> using Table = fhicl::Table<T>;
    using Comment  = fhicl::Comment;
    using Name     = fhicl::Name;
    using string   = std::string;
    using InputTag = art::InputTag;

    Atom<bool> CreateCAF { Name("CreateCAF"),
      Comment("Whether to produce an output file in CAF format"), true
    };

    Atom<bool> CreateFlatCAF { Name("CreateFlatCAF"),
      Comment("Whether to produce an output file in FlatCAF format"), true
    };

    Atom<bool> CreateBlindedCAF { Name("CreateBlindedCAF"),
      Comment("Whether to produce output files with one consisting of a fraction of events and the other consisting of the remainder of the events with critical information obscured"), true
    };

    Atom<std::string> CAFFilename { Name("CAFFilename"),
      Comment("Provide a string to override the automatic filename."), ""
    };

    Atom<std::string> FlatCAFFilename { Name("FlatCAFFilename"),
      Comment("Provide a string to override the automatic filename."), ""
    };

    Atom<float> PrescaleFactor { Name("PrescaleFactor"),
	Comment("Factor by which to prescale unblind events"), 10
    };

    Atom<int> POTBlindSeed { Name("POTBlindNum"),
	Comment("Integer used to derive POT scaling factor for blind events"), 655277
    };

    Atom<std::string> DetectorOverride { Name("DetectorOverride"),
      Comment("Override the automatically detectected detector using 'sbnd' or 'icarus'. This parameter should usually be unset - ''"),
      ""
    };

    Atom<string> DataTier        { Name("DataTier") };
    Atom<string> FileExtension   { Name("FileExtension"), ".caf.root" };
    Atom<string> FlatCAFFileExtension { Name("FlatCAFFileExtension"), ".flat.caf.root" };
    Atom<string> UnblindFileExtension   { Name("UnblindFileExtension"), ".Unblind.DONOTLOOK.dum" };
    Atom<string> BlindFileExtension { Name("BlindFileExtension"), ".Blind.OKTOLOOK.dum" };
    Atom<string> PrescaleFileExtension { Name("PrescaleFileExtension"), ".Prescaled.OKTOLOOK.dum" };

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

    Atom<string> BNBPOTDataLabel {
      Name("BNBPOTDataLabel"),
      Comment("Label of BNBRetriever module"),
      "bnbinfo"
    };

    Atom<string> NuMIPOTDataLabel {
      Name("NuMIPOTDataLabel"),
      Comment("Label of NuMIRetriever module"),
      "numiinfo"
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

    Atom<string> StubLabel {
      Name("StubLabel"),
      Comment("Base label of Stub producer."),
      "vertexStub"
    };

    Atom<string> FlashMatchLabel {
      Name("FlashMatchLabel"),
      Comment("Base label of flash match producer."),
      "fmatch" // same for icarus and sbnd
    };

    Atom<string> CRUMBSLabel {
      Name("CRUMBSLabel"),
      Comment("Base label of CRUMBS ID producer."),
      "crumbs"
    };

    Atom<string> CNNScoreLabel {
      Name("CNNScoreLabel"),
      Comment("Base label of CNNScore producer."),
      "cnnid"
    };

    Atom<bool> FillHits {
      Name("FillHits"),
      Comment("Label deciding if you want to fill SRHits"),
      false // default is false; set to true if you want SRHits filled
    };

    Atom<string> HitLabel {
      Name("HitLabel"),
      Comment("Base label of the TPC Hit producer."),
      "gaushit"
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

    Atom<string> ShowerRazzleLabel {
      Name("ShowerRazzleLabel"),
      Comment("Base label of shower mva particle-id producer."),
      "pandoraShowerRazzle"
    };

    Atom<string> RecoShowerSelectionLabel {
      Name("RecoShowerSelectionLabel"),
      Comment("Base label of shower selection vars producer."),
      "pandoraShowerSelectionVars"
    };

    Atom<string> ShowerCosmicDistLabel {
      Name("ShowerCosmicDistLabel"),
      Comment("Base label of shower selection vars producer."),
      "pandoraShowerCosmicDist"
    };

    Atom<string> TrackCaloLabel {
      Name("TrackCaloLabel"),
      Comment("Base label of track calorimetry producer."),
      "pandoraCalo"
    };

    Atom<string> TrackChi2PidLabel {
      Name("TrackChi2PidLabel"),
      Comment("Base label of track chi2 particle-id producer."),
      "pandoraPid"
    };

    Atom<string> TrackScatterClosestApproachLabel {
      Name("TrackScatterClosestApproachLabel"),
      Comment("Base label of track track scatter closestapproach producer."),
      "pandoraTrackClosestApproach"
    };

    Atom<string> TrackStoppingChi2FitLabel {
      Name("TrackStoppingChi2FitLabel"),
      Comment("Base label of track stopping chi2 fit producer."),
      "pandoraTrackStoppingChi2"
    };

    Atom<string> TrackDazzleLabel {
      Name("TrackDazzleLabel"),
      Comment("Base label of track mva particle-id producer."),
      "pandoraTrackDazzle"
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

    Atom<string> OpFlashLabel {
      Name("OpFlashLabel"),
      Comment("Label of PMT flash."),
      "OpFlash"
    };

    Atom<long long> CRTSimT0Offset {
      Name("CRTSimT0Offset"),
      Comment("start of beam gate/simulation time in the simulated CRT clock"),
      0,
    };

    Atom<art::InputTag> TriggerLabel {
      Name("TriggerLabel"),
      Comment("Label of trigger."),
      "daqTrigger"
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
      "simdrift"
    };

    Atom<bool> FillTrueParticles {
      Name("FillTrueParticles"),
      Comment("Whether to fill the rec.true_particles branch. The information on true particles"
              " will still be stored for the neutirno primaries and for trk/shw truth matching."),
      true
    };

    Sequence<std::string> SystWeightLabels {
      Name("SystWeightLabels"),
      Comment("Labels for EventWeightMap objects for mc.nu.wgt")
    };

    Atom<bool> FillHitsAllSlices {
      Name("FillHitsAllSlices"),
      Comment("Fill per-hit information in all reconstructed slices."),
      false
    };

    Atom<bool> FillHitsNeutrinoSlices {
      Name("FillHitsNeutrinoSlices"),
      Comment("Fill per-hit information in neutrino ID-d reconstructed slices."),
      true
    };

    Atom<float> TrackHitFillRRStartCut {
      Name("TrackHitFillRRStartCut"),
      Comment("How long from the start of a track to save calo-point information. Set to -1 to save nothing"),
      5.
    };

    Atom<float> TrackHitFillRREndCut {
      Name("TrackHitFillRREndCut"),
      Comment("How long from the end of a track to save calo-point information. Set to -1 to save nothing"),
      25.
    };

    Atom<bool> ReferencePMTFromTriggerToBeam {
      Name("ReferencePMTFromTriggerToBeam"),
      Comment("Whether to switch the reference time of PMT reco from 'trigger' to 'beam spill' time."),
      true
    };

    Atom<bool> ReferenceCRTT0ToBeam {
      Name("ReferenceCRTT0ToBeam"),
      Comment("Whether to switch the reference time of CRT T0 reco to the 'beam spill' time."),
      true
    };

    Atom<bool> ReferenceCRTT1FromTriggerToBeam {
      Name("ReferenceCRTT1FromTriggerToBeam"),
      Comment("Whether to switch the reference time of CRT T1 reco from 'trigger' to the 'beam spill' time."),
      true
    };

  };
}

#endif
