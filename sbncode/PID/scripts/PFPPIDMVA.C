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
    pfpTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root");
    pfpTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intrnue.root");
    pfpTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intime.root");

    TCut generic_cuts = "!unambiguousSlice && trk_length > 5 && showerEnergy > 10";

    TTree *electronTree = pfpTree->CopyTree("std::abs(truePDG)==11" + generic_cuts);
    TTree *muonTree     = pfpTree->CopyTree("std::abs(truePDG)==13" + generic_cuts);
    TTree *photonTree   = pfpTree->CopyTree("std::abs(truePDG)==22" + generic_cuts);
    TTree *pionTree     = pfpTree->CopyTree("std::abs(truePDG)==211" + generic_cuts);
    TTree *protonTree   = pfpTree->CopyTree("std::abs(truePDG)==2212" + generic_cuts);
    TTree *otherTree    = pfpTree->CopyTree("std::abs(truePDG)!=11 && std::abs(truePDG)!=13 && std::abs(truePDG)!=22 && std::abs(truePDG)!=211 && std::abs(truePDG)!=2212" + generic_cuts);

    dataloader->AddTree(electronTree, "Electron");
    dataloader->AddTree(muonTree, "Muon");
    dataloader->AddTree(photonTree, "Photon");
    dataloader->AddTree(pionTree, "Pion");
    dataloader->AddTree(protonTree, "Proton");
    dataloader->AddTree(otherTree, "Other");

    dataloader->AddVariable("pfp_numDaughters", "PFP N Daughters", "", 'F', 0, 5);
    dataloader->AddVariable("pfp_maxDaughterHits", "PFP Max Daughter Hits", "", 'F', 0, 500);
    dataloader->AddVariable("pfp_trackScore", "PFP Track Score", "", 'F', 0, 1);
    dataloader->AddVariable("pfp_chargeEndFrac", "PFP Charge End Fraction", "", 'F', 0, 1);
    dataloader->AddVariable("pfp_chargeFracSpread", "PFP Charge Fraction Spread", "", 'F', 0, 3);
    dataloader->AddVariable("pfp_linearFitDiff", "PFP Linear Fit Difference", "", 'F', 0, .5);
    dataloader->AddVariable("pfp_linearFitLength", "PFP Linear Fit Length", "cm", 'F', 0, 500);
    dataloader->AddVariable("pfp_linearFitGapLength", "PFP Linear Fit Gap Length", "cm", 'F', 0, 1);
    dataloader->AddVariable("pfp_linearFitRMS", "PFP Linear Fit RMS", "", 'F', 0, 4);
    dataloader->AddVariable("pfp_openAngleDiff", "PFP Opening Angle Difference", "#circ", 'F', 0, 50);
    dataloader->AddVariable("pfp_secondaryPCARatio", "PFP Secondary PCA Ratio", "", 'F', 0, 1);
    dataloader->AddVariable("pfp_tertiaryPCARatio", "PFP Tertiary PCA Ratio", "", 'F', 0, 1);
    dataloader->AddVariable("pfp_vertexDist", "PFP Vertex Distance", "cm", 'F', 0, 500);

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
    dataloader->AddVariable("shw_sqrtEnergyDensity", "Shower Sqrt Energy Density", "", 'F', 0, 10);


    const TCut baseCut("(abs(trackStartX) < 175 && abs(trackStartY) < 175 && trackStartZ > 25 && trackStartZ < 450"
		       " && abs(showerStartX) < 175 && abs(showerStartY) < 175 && showerStartZ > 25 && showerStartZ < 450"
		       " && recoPrimary == 1 && energyPurity > 0.5 && energyComp > 0.5 && trackContained && showerContained)");

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
        TMVA::TMVAMultiClassGui("Razzled" + instance + ".root");
}
