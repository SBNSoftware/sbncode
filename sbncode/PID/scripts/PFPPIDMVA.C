#include "TChain.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Tools.h"

void PFPPIDMVA()
{
  gStyle->SetOptStat(0);
    
  TString instance = "";
    
  TFile *outputFile = TFile::Open("Razzled" + instance + ".root", "RECREATE");

  TMVA::Factory *factory       = new TMVA::Factory("Razzled" + instance, outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=multiclass");
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("Razzled" + instance);

  TChain *pfpTree = new TChain("razzled/pfpTree");
  //  pfpTree->Add("/ADD/WHATEVER/SAMPLES/YOU'RE/USING/HERE");
  // Recommend using a combination of a rockbox sample & an intrinsic electron neutrino sample.
  pfpTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  pfpTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  pfpTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  TCut generic_cuts = "!unambiguousSlice && (trk_length > 3 || showerEnergy > 10)";

  TTree *electronTree = pfpTree->CopyTree("std::abs(truePDG)==11" + generic_cuts);
  TTree *muonTree     = pfpTree->CopyTree("std::abs(truePDG)==13" + generic_cuts);
  TTree *photonTree   = pfpTree->CopyTree("std::abs(truePDG)==22" + generic_cuts);
  TTree *pionTree     = pfpTree->CopyTree("std::abs(truePDG)==211" + generic_cuts);
  TTree *protonTree   = pfpTree->CopyTree("std::abs(truePDG)==2212" + generic_cuts);

  dataloader->AddTree(electronTree, "Electron");
  dataloader->AddTree(muonTree, "Muon");
  dataloader->AddTree(photonTree, "Photon");
  dataloader->AddTree(pionTree, "Pion");
  dataloader->AddTree(protonTree, "Proton");

<<<<<<< HEAD
  dataloader->AddVariable("pfp_numDaughters", "PFP N Daughters", "", 'F', 0, 5);
  dataloader->AddVariable("pfp_maxDaughterHits", "PFP Max Daughter Hits", "", 'F', 0, 500);
  dataloader->AddVariable("pfp_trackScore", "PFP Track Score", "", 'F', 0, 1);

  dataloader->AddVariable("trk_length", "Track Length", "cm", 'F', 0, 250);
  dataloader->AddVariable("trk_chi2PIDMuon", "Track Chi2 PID Muon", "", 'F', 0, 80);
  dataloader->AddVariable("trk_chi2PIDProton", "Track Chi2 PID Proton", "", 'F', 0, 300);
  dataloader->AddVariable("trk_chi2PIDMuonPionDiff", "Track Chi2 PID Muon-Pion", "", 'F', 0, 80);
  dataloader->AddVariable("trk_mcsScatterMean", "Track Mean MCS Scattering Angle", "mRad", 'F', 0, 600);
  dataloader->AddVariable("trk_mcsScatterMaxRatio", "Track Max/Mean MCS Scattering Angle Ratio", "", 'F', 0, 800);
  dataloader->AddVariable("trk_meanDCA", "Track Mean DCA", "cm", 'F', 0, 10);
  dataloader->AddVariable("trk_stoppingdEdxChi2Ratio", "Track Stopping Chi2Ratio", "", 'F', 0, 5);
  dataloader->AddVariable("trk_chi2Pol0dEdxFit", "Track Fitted Pol0 dE/dx", "MeV/cm", 'F', 0, 20);
  dataloader->AddVariable("trk_momDiff", "Track Momentum Agreement", "", 'F', 0, 10);

  dataloader->AddVariable("shw_bestdEdx", "Shower Best Plane dEdx", "MeV/cm", 'F', 0, 10);
  dataloader->AddVariable("shw_convGap", "Shower Conversion Gap", "cm", 'F', 0, 10);
  dataloader->AddVariable("shw_openAngle", "Shower Opening Angle", "rad", 'F', 0, 10);
  dataloader->AddVariable("shw_modHitDensity", "Shower Modified Hit Density", "", 'F', 0, 10);
  dataloader->AddVariable("shw_sqrtEnergyDensity > 2.5 ? 2.5 : shw_sqrtEnergyDensity", "Shower Sqrt Energy Density", "", 'F', 0, 10);
=======
  dataloader->AddVariable("pfp_numDaughters", "PFP N Daughters", "", 'F');
  dataloader->AddVariable("pfp_maxDaughterHits", "PFP Max Daughter Hits", "", 'F');
  dataloader->AddVariable("pfp_trackScore", "PFP Track Score", "", 'F');

  dataloader->AddVariable("trk_length", "Track Length", "cm", 'F');
  dataloader->AddVariable("trk_chi2PIDMuon", "Track Chi2 PID Muon", "", 'F');
  dataloader->AddVariable("trk_chi2PIDProton", "Track Chi2 PID Proton", "", 'F');
  dataloader->AddVariable("trk_chi2PIDMuonPionDiff", "Track Chi2 PID Muon-Pion", "", 'F');
  dataloader->AddVariable("trk_mcsScatterMean", "Track Mean MCS Scattering Angle", "mRad", 'F');
  dataloader->AddVariable("trk_mcsScatterMaxRatio", "Track Max/Mean MCS Scattering Angle Ratio", "", 'F');
  dataloader->AddVariable("trk_meanDCA", "Track Mean DCA", "cm", 'F');
  dataloader->AddVariable("trk_stoppingdEdxChi2Ratio", "Track Stopping Chi2Ratio", "", 'F');
  dataloader->AddVariable("trk_chi2Pol0dEdxFit", "Track Fitted Pol0 dE/dx", "MeV/cm", 'F');
  dataloader->AddVariable("trk_momDiff", "Track Momentum Agreement", "", 'F');

  dataloader->AddVariable("shw_bestdEdx", "Shower Best Plane dEdx", "MeV/cm", 'F');
  dataloader->AddVariable("shw_convGap", "Shower Conversion Gap", "cm", 'F');
  dataloader->AddVariable("shw_openAngle", "Shower Opening Angle", "rad", 'F');
  dataloader->AddVariable("shw_modHitDensity", "Shower Modified Hit Density", "", 'F');
  dataloader->AddVariable("shw_sqrtEnergyDensity > 2.5 ? 2.5 : shw_sqrtEnergyDensity", "Shower Sqrt Energy Density", "", 'F');
>>>>>>> f817b6d838f3b16bc31ac46c95a5c2145a49b155

  const TCut baseCut("(abs(trackStartX) < 180 && abs(trackStartY) < 180 && trackStartZ > 10"
                     " && trackStartZ < 450 && abs(showerStartX) < 180 && abs(showerStartY) < 180"
                     " && showerStartZ > 10 && showerStartZ < 450 && recoPrimary == 1"
                     " && energyPurity > 0.5 && energyComp > 0.5)");

  dataloader->PrepareTrainingAndTestTree(baseCut, "");

  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                      "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.50::BaggedSampleFraction=0.60"
                      ":nCuts=100:MaxDepth=3:DoBoostMonitor");

  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  outputFile->Close();

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch())
    TMVA::TMVAMultiClassGui("Razzled" + instance + ".root");
}
