#pragma once 

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include "SBNAna/Vars/Binnings.h"
#include "SBNAna/Vars/NueVars.h"
#include "SBNAna/Vars/Vars.h"
#include "SBNAna/Cuts/Cuts.h"
#include "SBNAna/Cuts/NueCuts.h"
#include "SBNAna/Cuts/TruthCuts.h"

using namespace ana;

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

const Binning kGapBinning       = Binning::Simple(22,-1,10);
const Binning kDensityBinning 	= Binning::Simple(20,0.,10);
const Binning kOpenAngleBinning = Binning::Simple(20,0.,2);
const Binning kLengthBinning  	= Binning::Simple(40,0.,200);
const Binning kPEBinning        = Binning::Simple(70,0.,1400);
const Binning kTimeBinning      = Binning::Simple(70,0.,3500);


// // In this example, we are making the following Spectra
std::vector<PlotDef> plots = 
  {{"count", "", Binning::Simple(3,0,3), kCounting},
   {"bestenergy", "Best plane energy (MeV)", kLowEnergyBinning,	kRecoShower_BestEnergy},  
   {"conversion", "Conversion gap (cm)",     kGapBinning,       kRecoShower_ConversionGap},  
   {"density",    "Shower density (MeV/cm)", kDensityBinning,   kRecoShower_Density},
   {"energy",     "Energy (MeV)",            kNueEnergyBinning, kRecoShower_Energy},
   {"length",     "Length (cm)",             kLengthBinning,    kRecoShower_Length},
   {"openangle",  "Opening angle",           kOpenAngleBinning,	kRecoShower_OpenAngle},
   {"startx", "Shower start position X (cm)", kPositionXFDBinning, kRecoShower_StartX},
   {"starty", "Shower start position Y (cm)", kPositionYFDBinning, kRecoShower_StartY},
   {"startz", "Shower start position Z (cm)", kPositionZFDBinning, kRecoShower_StartZ},
   {"endx",   "Shower end position X (cm)",   kPositionXFDBinning, kRecoShower_EndX},
   {"endy",   "Shower end position Y (cm)",   kPositionYFDBinning, kRecoShower_EndY},
   {"endz",   "Shower end position Z (cm)",   kPositionZFDBinning, kRecoShower_EndZ},
   {"vtxx", "Slice vertex X (cm)", kPositionXFDBinning, kSlcVtxX},
   {"vtxy", "Slice vertes Y (cm)", kPositionYFDBinning, kSlcVtxY},
   {"vtxz", "Slice vertex Z (cm)", kPositionZFDBinning, kSlcVtxZ}
 };

std::vector<PlotDefSpill> plots_spill =
  {{"crtx",   "CRT Hit Position X (cm)", kPositionXFDBinning, kCRTHitX},
  {"crty",    "CRT Hit Position Y (cm)", kPositionYFDBinning, kCRTHitY},
  {"crtz",    "CRT Hit Position Z (cm)", kPositionZFDBinning, kCRTHitZ},
  {"crttime", "CRT Hit Time (xs)", kTimeBinning, kCRTHitTime},
  {"crtpe",   "CRT PE", kPEBinning, kCRTHitPE}
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
  SpillCut cut = kSpillNoCut;
  int color = kBlack;
};

// In this example, our signal is Nue cc
const Cut kSig = kIsNue && !kIsNC;
const Cut kContained = kNueContainedFD;

// Define a basic cut for testing
const Cut kShortShower(
    [](const caf::SRSliceProxy* slc){
      if(slc->reco.nshw>0){ // check if there is actually a shower
        return ( slc->reco.shw[0].len < 30.);
      }
      return false;
  }
  );


std::vector<SelDef> sels =
  {{"nocut", "All Slices",  kNoCut, kBlack},
   {"sig",   "True NumuCC", kSig,   kRed+1},
   {"bkg",   "Not NumuCC",  !kSig,  kAzure+2},
   {"cont",     "Contained",        kContained,          kBlack},
   {"sig_cont", "Contained Signal", kContained && kSig,  kRed+1},
   {"bkg_cont", "Contained Bkg",    kContained && !kSig, kAzure+2},
   {"uncont",     "Uncontained",        !kContained,          kBlack},
   {"sig_uncont", "Uncontained Signal", !kContained && kSig,  kRed+1},
   {"bkg_uncont", "Uncontained Bkg",    !kContained && !kSig, kAzure+2}
  };

std::vector<SelDefSpill> sels_spill =
  {{"nocut_spill", "All Slices",         kSpillNoCut,   kBlack},
   {"event_spill", "First Events",       kFirstEvents,  kRed+1},
   {"flash_spill", "Pass flash Trigger", kFlashTrigger, kAzure+2}
  };
