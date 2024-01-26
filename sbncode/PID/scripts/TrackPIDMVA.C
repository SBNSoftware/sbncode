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

void TrackPIDMVA()
{
  gStyle->SetOptStat(0);

  TFile *outputFile = TFile::Open("TrackPIDMVA.root", "RECREATE");

  TMVA::Factory *factory       = new TMVA::Factory("TrackPIDMVA", outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=multiclass");
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

  TChain *trackTree = new TChain("dazzle/trackTree");
  //  trackTree->Add("/ADD/WHATEVER/SAMPLES/YOU'RE/USING/HERE");
  trackTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  trackTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  trackTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  TTree *muonTree   = trackTree->CopyTree("std::abs(truePdg)==13");
  TTree *pionTree   = trackTree->CopyTree("std::abs(truePdg)==211");
  TTree *protonTree = trackTree->CopyTree("std::abs(truePdg)==2212");
  TTree *otherTree  = trackTree->CopyTree("std::abs(truePdg)!=13 && std::abs(truePdg)!=211 && std::abs(truePdg)!=2212");

  dataloader->AddTree(muonTree, "Muon");
  dataloader->AddTree(pionTree, "Pion");
  dataloader->AddTree(protonTree, "Proton");
  dataloader->AddTree(otherTree, "Other");

  dataloader->AddVariable("recoLen", "Reco. Length", "cm", 'F', 0, 250);
  dataloader->AddVariable("chi2PIDMuon", "Chi2 PID Muon", "", 'F', 0, 80);
  dataloader->AddVariable("chi2PIDProton", "Chi2 PID Proton", "", 'F', 0, 300);
  dataloader->AddVariable("chi2PIDMuonPionDiff", "Chi2 PID Muon-Pion", "", 'F', 0, 80);
  dataloader->AddVariable("mcsScatterMean", "Mean MCS Scattering Angle", "mRad", 'F', 0, 600);
  dataloader->AddVariable("mcsScatterMaxRatio", "Max/Mean MCS Scattering Angle Ratio", "", 'F', 0, 800);
  dataloader->AddVariable("meanDCA", "Mean DCA", "cm", 'F', 0, 10);
  dataloader->AddVariable("stoppingChi2Ratio", "Stopping CHi2Ratio", "", 'F', 0, 5);
  dataloader->AddVariable("chi2Pol0Fit", "Fitted Pol0 dE/dx", "MeV/cm", 'F', 0, 20);
  dataloader->AddVariable("pDiff", "Momentum Agreement", "", 'F', 0, 10);
  dataloader->AddVariable("numDaughters", "Num Daughters", "", 'I', 0, 5);
  dataloader->AddVariable("maxDaughterHits", "Daughter Hits", "", 'I', 0, 300);

  const TCut baseCut("(abs(startX) < 175 && abs(startY) < 175 && startZ > 25 && startZ < 450"
                     "&& recoPrimary == 1 && recoLen>10 && trackScore > 0.5 && "
                     "energyPurity > 0.5 && energyComp > 0.5 && recoContained)");

  dataloader->PrepareTrainingAndTestTree(baseCut, "");


  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                      "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.50::BaggedSampleFraction=0.60:nCuts=100:MaxDepth=3:DoBoostMonitor");
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  outputFile->Close();

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch())
    TMVA::TMVAMultiClassGui("TrackPIDMVA.root");
}
