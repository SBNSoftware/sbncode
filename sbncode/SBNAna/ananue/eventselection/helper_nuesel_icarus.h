#pragma once 

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include "SBNAna/Vars/Binnings.h"
#include "SBNAna/Vars/NueVars.h"
#include "SBNAna/Vars/TruthVars.h"
#include "SBNAna/Vars/Vars.h"
#include "SBNAna/Cuts/Cuts.h"
#include "SBNAna/Cuts/NueCuts.h"
#include "SBNAna/Cuts/TruthCuts.h"

using namespace ana;

Color_t color_nue   = kBlue-7;
Color_t color_numu  = kGreen+1;
Color_t color_nc    = kMagenta;
Color_t color_other = kOrange+8;
Color_t color_cos   = kGray+2;
Style_t line_nue    = kSolid;
Style_t line_numu   = kSolid;
Style_t line_nc     = kDashed;
Style_t line_other  = kDotted;
Style_t line_cos    = kSolid;

// These are examples of useful structs to
// use when making a bunch of Spectra
struct PlotDef
{
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  Var var = kCounting;
};

struct PlotDefSpill
{
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  SpillMultiVar var = kCRTHitX;
};

const Binning kEnergyBinning    = Binning::Simple(40,0.,3000.); // to define
const Binning kDedxBinning      = Binning::Simple(20,0.,10); // to define
const Binning kGapBinning       = Binning::Simple(20,0.,10);
const Binning kDensityBinning 	= Binning::Simple(25,0.,10);
const Binning kOpenAngleBinning = Binning::Simple(20,0.,2);
const Binning kLengthBinning  	= Binning::Simple(40,0.,200);
const Binning kPEBinning        = Binning::Simple(60,0.,600);
const Binning kTimeBinning      = Binning::Simple(155,-1550.,1550.);
const Binning kFlashBinning     = Binning::Simple(40,-6.f,34.f);

// // In this example, we are making the following Spectra
std::vector<PlotDef> plots = 
  {{"count",      "Number of slices",             Binning::Simple(3,0,3), kCounting},
   {"openangle",  "Opening angle",                kOpenAngleBinning,      kRecoShower_OpenAngle},
   {"startx",     "Shower start position X (cm)", kPositionXFDBinning,    kRecoShower_StartX},
   {"starty",     "Shower start position Y (cm)", kPositionYFDBinning,    kRecoShower_StartY},
   {"startz",     "Shower start position Z (cm)", kPositionZFDBinning,    kRecoShower_StartZ},
   {"endx",       "Shower end position X (cm)",   kPositionXFDBinning,    kRecoShower_EndX},
   {"endy",       "Shower end position Y (cm)",   kPositionYFDBinning,    kRecoShower_EndY},
   {"endz",       "Shower end position Z (cm)",   kPositionZFDBinning,    kRecoShower_EndZ},
   {"vtxx",       "Slice vertex X (cm)",          kPositionXFDBinning,    kSlcVtxX},
   {"vtxy",       "Slice vertes Y (cm)",          kPositionYFDBinning,    kSlcVtxY},
   {"vtxz",       "Slice vertex Z (cm)",          kPositionZFDBinning,    kSlcVtxZ},
   {"conversion", "Conversion gap (cm)",          kGapBinning,            kRecoShower_ConversionGap},  
   {"bestdedx",   "Best plane dEdx (MeV)",        kDedxBinning,           kRecoShower_BestdEdx},  
   {"bestenergy", "Best plane energy (MeV)",      kEnergyBinning,         kRecoShower_BestEnergy},  
   {"density",    "Shower density (MeV/cm)",      kDensityBinning,        kRecoShower_Density},
   {"energy",     "Shower energy (MeV)",          kNueEnergyBinning,      kRecoShower_Energy},
   {"length",     "Length (cm)",                  kLengthBinning,         kRecoShower_Length},
   {"nuscore",    "Pandora #nu score",            kBDTBinning,            kSlcNuScore},
   {"flashscore", "Flash score",                  kFlashBinning,          kSlcFlashScore},
   {"truthenergy","True #nu energy (GeV)",        kTruthLowEnergyBinning, kTruthEnergy},
 };

std::vector<PlotDefSpill> plots_spill =
  {{"crtx",   "CRT Hit Position X (cm)", kCRTXFDBinning, kCRTHitX},
  {"crty",    "CRT Hit Position Y (cm)", kCRTYFDBinning, kCRTHitY},
  {"crtz",    "CRT Hit Position Z (cm)", kCRTZFDBinning, kCRTHitZ},
  {"crttime", "CRT Hit Time (#mus)",     kTimeBinning,   kCRTHitTimeFD},
  {"crtpe",   "CRT PE",                  kPEBinning,     kCRTHitPE}
};

// Selection Struc
struct SelDef
{
  std::string suffix = "";
  std::string label = "";
  Cut cut = kNoCut;
  int color = kBlack;
};

struct SelDefSpill
{
  std::string suffix = "";
  std::string label = "";
  SpillCut cut = kNoSpillCut;
  int color = kBlack;
};

const Cut kNuECC    = kIsNue && !kIsNC;
const Cut kNuMuCC   = kIsNumu && !kIsNC;
const Cut kNC       = kIsNC;
const Cut kTotal    = kNoCut; // Get the total and substract everything else when plotting
const Cut kIsCosmic = !kHasNu;
const Cut kOtherNu  = kHasNu && !(kNuECC || kNuMuCC || kNC);

// Step by step cuts
const Cut kContained = kContainedFD;
const Cut kRecoCut   = kRecoShowerFD;
const Cut kFullCut   = kNueFD;

// N-1 cuts: apply all except one cut
// // const Cut kAllCuts      = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kN1Contained  = kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShowerFD;
const Cut kN1Flash      = kContained && kNuePandoraScoreFDCut && kRecoShowerFD;
const Cut kN1Pandora    = kContained && kNueFlashScoreFDCut && kRecoShowerFD;
const Cut kN1Reco       = kContained && kNueFlashScoreFDCut && kNuePandoraScoreFDCut;
const Cut kN1RecoShower = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kN1NumShowers = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kN1Dedx       = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kN1ConvGap    = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kN1TrkLen     = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kN1Density    = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerEnergyCut;
const Cut kN1Energy     = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut;

std::vector<SelDef> types =
{
  {"nue",   "NuE CC",  kSlcIsRecoNu&&kNuECC,     color_nue},
  {"numu",  "NuMu CC", kSlcIsRecoNu&&kNuMuCC,    color_numu},
  {"nunc",  "NC",      kSlcIsRecoNu&&kNC,        color_nc},
  {"other", "Other",   kSlcIsRecoNu&&kOtherNu,   color_other},
  {"total", "Total",   kSlcIsRecoNu&&kTotal,     color_other},
  {"cosmic", "Cosmic", kSlcIsRecoNu&&kIsCosmic,  color_cos}
};

std::vector<SelDef> sels ={
  {"nocut",      "No cut",               kNoCut,            kBlack},
  {"cont",       "Containment",          kContained,        kBlack},
  {"flash",      "Flash score",          kSlcFlashMatchCut, kBlack},
  {"pandnu",     "Neutrino score",       kSlcNuScoreCut,    kBlack},
  // Subcuts that go into the fd reco cut
  {"recoshw",    "Reconstructed shower", kRecoShower,       kBlack},
  {"nshws",      "Number of showers",    kNueNumShowersCut, kBlack},
  {"dedx",       "Shower dE/dx",         kShowerdEdxCut,    kBlack},
  {"convgap",    "Conversion gap",       kShowerConvGapCut, kBlack},
  {"trklen",     "Track lenght",         kNueTrackLenCut,   kBlack},
  {"density",    "Shower density",       kShowerDensityCut, kBlack},
  {"energy",     "Shower energy",        kShowerEnergyCut,  kBlack},
  {"recocut",    "Reconstruction (all)", kRecoCut,          kBlack},
  {"everything", "Full selection",       kFullCut,          kBlack},
  // N-1 cuts
  {"N1cont",    "N1 Containment",          kN1Contained,  kBlack},
  {"N1flash",   "N1 Flash score",          kN1Flash,      kBlack},
  {"N1pandnu",  "N1 Neutrino score",       kN1Pandora,    kBlack},
  {"N1recoshw", "N1 Reconstructed shower", kN1RecoShower, kBlack},
  {"N1nshws",   "N1 Number of showers",    kN1NumShowers, kBlack},
  {"N1dedx",    "N1 Shower dE/dx",         kN1Dedx,       kBlack},
  {"N1convgap", "N1 Conversion gap",       kN1ConvGap,    kBlack},
  {"N1trklen",  "N1 Track length",         kN1TrkLen,     kBlack},
  {"N1density", "N1 Shower density",       kN1Density,    kBlack},
  {"N1recocut", "N1 Reconstruction (all)", kN1Reco,       kBlack},
  {"N1energy",  "N1 Shower energy",        kN1Energy,     kBlack}
  };

std::vector<SelDefSpill> sels_spill =
  {{"nocut_spill", "All Slices",         kNoSpillCut,   kBlack},
   {"crtveto_spill", "CRTVeto",          kCRTHitVetoFD, kBlue}
  };
