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

void ShowerPIDMVA()
{
    gStyle->SetOptStat(0);

    TFile *outputFile = TFile::Open("ShowerPIDMVA.root", "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("ShowerPIDMVA", outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=multiclass");
    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

    TChain *showerTree = new TChain("razzle/showerTree");
    showerTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root");
    showerTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intrnue.root");
    showerTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intime.root");

    TTree *electronTree = showerTree->CopyTree("std::abs(truePdg)==11");
    TTree *photonTree   = showerTree->CopyTree("std::abs(truePdg)==22");
    TTree *otherTree    = showerTree->CopyTree("std::abs(truePdg)!=11 && std::abs(truePdg)!=22");

    dataloader->AddTree(electronTree, "Electron");
    dataloader->AddTree(photonTree, "Photon");
    dataloader->AddTree(otherTree, "Other");

    dataloader->AddVariable("bestdEdx", "Best Plane dEdx", "MeV/cm", 'F', 0, 10);
    dataloader->AddVariable("convGap", "Conversion Gap", "cm", 'F', 0, 10);
    dataloader->AddVariable("openAngle", "Opening Angle", "rad", 'F', 0, 10);
    dataloader->AddVariable("modHitDensity", "Modified Hit Density", "", 'F', 0, 10);
    dataloader->AddVariable("sqrtEnergyDensity", "Sqrt Energy Density", "", 'F', 0, 10);
    // dataloader->AddVariable("logEnergyDensity", "Sqrt Energy Density", "cm^-1", 'F', 0, 10);

    // dataloader->AddSpectator("trueP", "True Momentum", "GeV");
    // dataloader->AddSpectator("bestEnergy", "Reco Energy", "MeV");

    const TCut baseCut("(abs(startX) < 175 && abs(startY) < 175 && startZ > 25 && startZ < 450 && recoPrimary == 1 && bestEnergy>100 && "
                       "energyPurity > 0.5 && energyComp > 0.5 && recoContained)");

    dataloader->PrepareTrainingAndTestTree(baseCut, "");

    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
        "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.50::BaggedSampleFraction=0.60:nCuts=100:MaxDepth=3:DoBoostMonitor");

    factory->TrainAllMethods();

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

    outputFile->Close();

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch())
        TMVA::TMVAMultiClassGui("ShowerPIDMVA.root");
}
