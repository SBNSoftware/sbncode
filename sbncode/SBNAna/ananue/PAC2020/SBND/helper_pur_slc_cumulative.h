#pragma once

#include "helper.h"

using namespace ana;

// These are examples of useful structs to
// use when making a bunch of Spectra
struct PlotDef {
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3, 0, 3);
  Var var = kCounting;
};

const Binning kBDTBinning = Binning::Simple(40, 0, 1.0);
const Binning kEnergyBinning = Binning::Simple(40, 0, 3000);
const Binning kEnergyBinningGeV = Binning::Simple(40, 0, 3);
const Binning kdEdxBinning = Binning::Simple(40, 0, 10);
const Binning kGapBinning = Binning::Simple(20, 0, 20);
const Binning kDensityBinning = Binning::Simple(20, 0., 20);
const Binning kOpenAngleBinning = Binning::Simple(20, 0., 1);
const Binning kLengthBinning = Binning::Simple(42, -5.01, 200);
const Binning kPEBinning = Binning::Simple(70, 0., 1400);
const Binning kTimeBinning = Binning::Simple(70, -3500, 3500);
const Binning kFlashBinning = Binning::Simple(40, -6.f, 34.f);
const Binning kPdgBinning = Binning::Simple(5000, -2500, 2500);

// // In this example, we are making the following Spectra
std::vector<PlotDef> plots = {
  { "count", "", Binning::Simple(3, 0, 3), kCounting },
  { "bestenergy", "Best plane energy (MeV)", kEnergyBinning, kRecoShower_BestEnergy },
  { "conversion", "Conversion gap (cm)", kGapBinning, kRecoShower_ConversionGap },
  { "dedx", "dEdx (MeV/cm)", kdEdxBinning, kRecoShower_BestdEdx },
  { "slength", "Shower Length (cm)", kLengthBinning, kRecoShower_Length },
  { "tlength", "Longest Track Length (cm)", kLengthBinning, kLongestTrackLength },
  { "openangle", "Opening angle (rad)", kOpenAngleBinning, kRecoShower_OpenAngle },
  { "density", "Shower density (MeV/cm)", kDensityBinning, kRecoShower_Density },
  { "startx", "Shower start position X (cm)", kPositionXNDBinning, kRecoShower_StartX },
  { "starty", "Shower start position Y (cm)", kPositionYNDBinning, kRecoShower_StartY },
  { "startz", "Shower start position Z (cm)", kPositionZNDBinning, kRecoShower_StartZ },
  { "vtxx", "Slice vertex X (cm)", kPositionXNDBinning, kSlcVtxX },
  { "vtxy", "Slice vertes Y (cm)", kPositionYNDBinning, kSlcVtxY },
  { "vtxz", "Slice vertex Z (cm)", kPositionZNDBinning, kSlcVtxZ },
  { "nuScore", "Slice Nu Score", kBDTBinning, kSlcNuScore },
  { "flashScore", "Slice Flash Score", kFlashBinning, kSlcFlashScore },
  { "vtxdist", "Vertex Distance (cm)", kGapBinning, kTruthVtxDistMag },
  { "trueshwpdg", "True Shower Pdg", kPdgBinning, kRecoShower_TruePdg },
  { "truetrkpdg", "True Track Pdg", kPdgBinning, kLongestTrackTruePdg },
};

// Selection Struc
struct SelDef {
  std::string suffix = "";
  std::string label = "";
  Cut cut = kNoCut;
  int color = kBlack;
};

std::vector<SelDef> types = {
  { "nuecc", "NuE CC", kSlcIsRecoNu&& kNuECC, kBlue + 1 },
  { "numucc", "NuMu CC", kSlcIsRecoNu&& kNuMuCC, kRed + 1 },
  { "nunc", "NC", kSlcIsRecoNu&& kNC, kOrange },
  { "othernu", "Other Nu", kSlcIsRecoNu&& kOtherNu, kOrange },
  { "cosmic", "Cosmic", kSlcIsRecoNu&& kCosmicRay, kBlack },
};

std::vector<SelDef> sels = {
  { "fvcut", "FV Cut", kFiducialVolumeND&& kVtxDistMagCut, kBlue },
  { "nuscorecut", "Nu Score Cut", kSlcNuScoreCut, kBlue },
  { "fmatchcut", "Flah Match Cut", kSlcFlashMatchCut, kBlue },
  { "recosel", "Reco Selection", kRecoNueSel, kBlue },
  { "trklencut", "Track Length Cut", kNueTrackLenCut, kBlue },
  { "shwdensitycut", "Shower Density Cut", kShowerDensityCut, kBlue },
  { "shwconvgapcut", "Conversion Gap Cut", kShowerConvGapCut, kBlue },
  { "shwdedxcut", "Shower dEdx Cut", kShowerdEdxCut, kBlue },
  { "fullsel", "Full Selection", kFullNueSel, kBlue },
};

// std::vector<SelDef> sels = {
//   { "fvcut", "FV Cut", kFiducialVolumeND&& kVtxDistMagCut, kBlue },
//   { "presel", "Pre Selection", kPreNueSel, kBlue },
//   { "recosel", "Reco Selection", kRecoNueSel, kBlue },
//   { "fullsel", "Full Selection", kFullNueSel, kBlue },
// };
