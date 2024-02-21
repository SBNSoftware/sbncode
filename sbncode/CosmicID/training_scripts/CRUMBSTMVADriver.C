#include <cstdlib>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void TrainCRUMBSInstance(const TString outDirName, TTree *inputTree, 
                         const TCut sigCut, const TCut backCut,
                         const bool useOpT0Finder = true, const bool useSimpleFlash = false,
                         const double signalWeight = 1.0, const double backgroundWeight = 1.0);

int CRUMBSTMVADriver(const bool useOpT0Finder = true, const bool useSimpleFlash = false)
{
  if(useOpT0Finder && useSimpleFlash)
    {
      std::cout << "Cannot use both flash matchers, please chose one..." << std::endl;
      return 0;
    }
  else if(!useOpT0Finder && !useSimpleFlash)
    {
      std::cout << "Would recommend using one of the flash matchers..." << std::endl;
      return 0;
    }

  TMVA::Tools::Instance();

  std::cout << std::endl;
  std::cout << "==> Start CRUMBS TMVA Training" << std::endl;

  TChain *inputTree = new TChain("crumbs/SliceTree");
  //  inputTree->Add("/ADD/WHATEVER/SAMPLES/YOU'RE/USING/HERE");
  // Recommend using a combination of a rockbox sample, an intrinsic electron neutrino sample
  // and an intime cosmics sample in roughly a 2:2:1 ratio
  inputTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  inputTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  inputTree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  if (!inputTree) {
    std::cout << "ERROR: could not access tree" << std::endl;
    exit(1);
  }
  std::cout << "--- CRUMBS TMVA Training     : Accessed tree: " << inputTree->GetName() << std::endl;

  TCut background      = "!strstr(matchedType,\"Nu\")";
  TCut inclusiveSignal = "strstr(matchedType,\"Nu\") && !strstr(matchedType,\"DirtNu\") && matchedPurity > 0.8 && matchedCompleteness > 0.8";
  TCut ccnumuSignal    = inclusiveSignal + "ccnc == 0 && abs(nutype) == 14";
  TCut ccnueSignal     = inclusiveSignal + "ccnc == 0 && abs(nutype) == 12";
  TCut ncSignal        = inclusiveSignal + "ccnc == 1";

  TrainCRUMBSInstance("CRUMBS_Inclusive", inputTree, inclusiveSignal, background, useOpT0Finder, useSimpleFlash);
  TrainCRUMBSInstance("CRUMBS_CCNuMu", inputTree, ccnumuSignal, background, useOpT0Finder, useSimpleFlash);
  TrainCRUMBSInstance("CRUMBS_CCNuE", inputTree, ccnueSignal, background, useOpT0Finder, useSimpleFlash);
  TrainCRUMBSInstance("CRUMBS_NC", inputTree, ncSignal, background, useOpT0Finder, useSimpleFlash);

  std::cout << std::endl;
  std::cout << "==> Finish CRUMBS TMVA Training" << std::endl;

  return 0;
}

void TrainCRUMBSInstance(const TString outDirName, TTree *inputTree, 
                         const TCut sigCut, const TCut backCut,
                         const bool useOpT0Finder, const bool useSimpleFlash,
                         const double signalWeight, const double backgroundWeight)
{
  std::cout << std::endl;
  std::cout << "--- CRUMBS TMVA Training     : Beginning instance: " << outDirName << std::endl;

  gSystem->Exec("mkdir -v " + outDirName);
  TMVA::DataLoader *dataloader=new TMVA::DataLoader(outDirName);

  TString outfileName( outDirName + "/CRUMBSTMVA.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "CrumbsTMVAClassification", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );

  dataloader->AddVariable("tpc_CRFracHitsInLongestTrack","Fraction of Hits in Longest Track (Cosmic Reco)","",'F');
  dataloader->AddVariable("tpc_CRLongestTrackDeflection","Longest Track Deflection (Cosmic Reco)","",'F');
  dataloader->AddVariable("tpc_CRLongestTrackDirY","Longest Track Y Direction (Cosmic Reco)","",'F');
  dataloader->AddVariable("tpc_CRNHitsMax","nHits (Cosmic Reco)","",'F');
  dataloader->AddVariable("tpc_NuEigenRatioInSphere","Eigen Ratio in sphere (Nu Reco)","",'F');
  dataloader->AddVariable("tpc_NuNFinalStatePfos","nPFOs (Nu Reco)","",'F');
  dataloader->AddVariable("tpc_NuNHitsTotal","nHits (Nu Reco)","",'F');
  dataloader->AddVariable("tpc_NuNSpacePointsInSphere","nSpacePoints in sphere (Nu Reco)","",'F');
  dataloader->AddVariable("tpc_NuVertexY","Vertex Y (Nu Reco)","cm",'F');
  dataloader->AddVariable("tpc_NuWeightedDirZ","Weighted Z Direction (Nu Reco)","",'F');
  dataloader->AddVariable("tpc_StoppingChi2CosmicRatio","Stopping Fits Chi2 Ratio","",'F');

  if(useOpT0Finder)
    {
      dataloader->AddVariable("pds_OpT0Score","OpT0 Score","",'F');
      dataloader->AddVariable("isinf(pds_OpT0MeasuredPE) ? -10000 : pds_OpT0MeasuredPE", "OpT0 Measured PE","",'F');
    }

  if(useSimpleFlash)
    {
      dataloader->AddVariable("pds_FMTotalScore","Total FM Score","",'F');
      dataloader->AddVariable("pds_FMPE","nPE in flash","",'F');
      dataloader->AddVariable("pds_FMTime","FM Time","#mu s",'F');
    }

  dataloader->AddVariable("crt_TrackScore","CRT Track Match Score","",'F');
  dataloader->AddVariable("crt_SPScore","CRT SpacePoint Match Score","",'F');
  dataloader->AddVariable("crt_TrackTime","CRT Track Match Time","#mu s",'F');
  dataloader->AddVariable("crt_SPTime","CRT SpacePoint Match Time","#mu s",'F');

  dataloader->AddSignalTree    (inputTree, signalWeight);
  dataloader->AddBackgroundTree(inputTree, backgroundWeight);

  dataloader->PrepareTrainingAndTestTree( sigCut, backCut, "SplitMode=random:!V" );
  
  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                       "!H:!V:NTrees=100:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=100");

  factory->TrainAllMethods();
  
  factory->TestAllMethods();
  
  factory->EvaluateAllMethods();
  
  outputFile->Close();
  
  std::cout << "--- CRUMBS TMVA Training     : Wrote root file: " << outputFile->GetName() 
            << std::endl;
  std::cout << "--- CRUMBS TMVA Training     : Finishing instance: " << outDirName << std::endl;

  delete factory;
  delete dataloader;

  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
}
