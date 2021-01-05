#pragma once

#include "helper.h"

using namespace ana;

const SpillVar kLeptonEnergy([](const caf::SRSpillProxy* sr) {
  if (sr->mc.nnu != 1)
    return -5.f;

  for (auto const& prim : sr->mc.nu[0].prim) {
    if (std::abs(prim.pdg) != 11 && std::abs(prim.pdg) != 13)
      continue;

    return (float)prim.startE;
  }
  return -5.f;
});

const SpillCut kCosmicSpill([](const caf::SRSpillProxy* sr) {
  return (sr->mc.nnu == 0);
});

const SpillCut kSpillTrueNuFV([](const caf::SRSpillProxy* sr) {
  if (kCosmicSpill(sr))
    return true;

  if (!kSpillSingleNu(sr))
    return false;

  return PtInVolAbsX(sr->mc.nu[0].position, fvndAbs);
});

const SpillCut kSpillNuEEnergyCut([](const caf::SRSpillProxy* sr) {
  if (kCosmicSpill(sr) || !(std::abs(sr->mc.nu[0].pdg) == 12 && sr->mc.nu[0].iscc))
    return true;

  for (auto const& prim : sr->mc.nu[0].prim) {
    if (std::abs(prim.pdg) != 11)
      continue;

    if (prim.startE > 0.2)
      return true;
  }

  return false;
});

const SpillCut kNuECCSpill([](const caf::SRSpillProxy* sr) {
  if (kCosmicSpill(sr) || !kSpillSingleNu(sr))
    return false;

  return (std::abs(sr->mc.nu[0].pdg) == 12 && sr->mc.nu[0].iscc && PtInVolAbsX(sr->mc.nu[0].position, fvndAbs) && kSpillNuEEnergyCut(sr));
});

const SpillCut kNuMuCCSpill([](const caf::SRSpillProxy* sr) {
  if (kCosmicSpill(sr) || !kSpillSingleNu(sr))
    return false;
  return (std::abs(sr->mc.nu[0].pdg) == 14 && sr->mc.nu[0].iscc && PtInVolAbsX(sr->mc.nu[0].position, fvndAbs));
});

const SpillCut kNCSpill([](const caf::SRSpillProxy* sr) {
  if (kCosmicSpill(sr) || !kSpillSingleNu(sr))
    return false;
  return ((bool)sr->mc.nu[0].isnc && PtInVolAbsX(sr->mc.nu[0].position, fvndAbs));
});

const SpillCut kOtherNuSpill([](const caf::SRSpillProxy* sr) {
  return (kSpillSingleNu(sr) && !(kCosmicSpill(sr) || kNuECCSpill(sr) || kNuMuCCSpill(sr) || kNCSpill(sr)));
});

struct PlotDef {
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3, 0, 3);
  SpillVar var = kNuEnergy;
};

const Binning kBDTBinning = Binning::Simple(40, 0, 1.0);
const Binning kEnergyBinning = Binning::Simple(40, 0, 3000);
const Binning kEnergyBinningGeV = Binning::Simple(30, 0, 3);
const Binning kdEdxBinning = Binning::Simple(40, 0, 10);
const Binning kGapBinning = Binning::Simple(20, 0, 10);
const Binning kDensityBinning = Binning::Simple(20, 0., 20);
const Binning kOpenAngleBinning = Binning::Simple(20, 0., 1);
const Binning kLengthBinning = Binning::Simple(42, -5.01, 100);
const Binning kPEBinning = Binning::Simple(70, 0., 1400);
const Binning kTimeBinning = Binning::Simple(70, -3500, 3500);
const Binning kFlashBinning = Binning::Simple(40, -6.f, 34.f);

std::vector<PlotDef> plots = {
  { "count", "", Binning::Simple(3, 0, 3), kSpillCounting },
  { "nuEnergy", "True Neutrino Energy [GeV]", kEnergyBinningGeV, kNuEnergy },
  { "leptonEnergy", "True Lepton Energy [GeV]", kEnergyBinningGeV, kLeptonEnergy }
};

struct SelDef {
  std::string suffix = "";
  std::string label = "";
  SpillCut cut = kNoSpillCut;
  int color = kBlack;
};

const SpillCut kSpillVtxDistCut([](const caf::SRSpillProxy* sr) {
  if (kCosmicSpill(sr))
    return true;

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0 || !kSlcIsRecoNu(&slc))
      continue;
    if (kVtxDistMagCut(&slc))
      ++counter;
  }

  return (bool)counter;
});

const Cut kSpillPreSelSlc = kSlcIsRecoNu && kPreNueSelND && kVtxDistMagCut && kSlcCompletenessCut;
// const Cut kSpillPreSelSlc = kSlcIsRecoNu && kPreNueSelND && kRecoComp;

const SpillCut kSpillPreSel([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kPreNueSelND(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kSpillRecoSel([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kPreNueSelND(&slc) && kRecoNueSel(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kSpillFullSel([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kPreNueSelND(&slc) && kRecoNueSel(&slc) && kFullNueSel(&slc))
      ++counter;
  }

  return counter;
});

std::vector<SelDef> types = {
  { "nuecc", "NuE CC", kNuECCSpill, kBlue + 1 },
  { "numucc", "NuMu CC", kNuMuCCSpill, kRed + 1 },
  { "nunc", "NC", kNCSpill, kOrange },
  { "othernu", "Other Nu", kOtherNuSpill, kOrange },
  // { "cosmic", "Cosmic", kCosmicSpill, kBlack },
};

// const SpillCut kNuESpillTruthCuts = kSpillSingleNu; //&& kSpillNuEEnergyCut
const SpillCut kNuESpillTruthCuts = kSpillVtxDistCut && kSpillSingleNu; //&& kSpillNuEEnergyCut

std::vector<SelDef> sels = {
  { "spill_trueshowerenergy", "True Shower Energy", kNuESpillTruthCuts, kBlue },
  { "spill_crtVeto", "CRT Veto", kNuESpillTruthCuts&& kCRTHitVetoND, kBlue },
  { "spill_presel", "Pre Selection", kNuESpillTruthCuts&& kCRTHitVetoND&& kSpillPreSel, kBlue },
  { "spill_recosel", "Reco Selection", kNuESpillTruthCuts&& kCRTHitVetoND&& kSpillRecoSel, kBlue },
  { "spill_fullsel", "Full Selection", kNuESpillTruthCuts&& kCRTHitVetoND&& kSpillFullSel, kBlue },
};
