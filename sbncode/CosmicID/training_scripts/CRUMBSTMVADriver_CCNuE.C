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

int CRUMBSTMVADriver_CCNuE()
{
  TMVA::Tools::Instance();

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  TFile *input(0);
  TString fname = "PATH/TO/YOUR/TREE/FILE/HERE";
  if (!gSystem->AccessPathName( fname )) {
    input = TFile::Open( fname );
  }
  if (!input) {
    std::cout << "ERROR: could not open data file" << std::endl;
    exit(1);
  }
  std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;

  TTree *signalTree     = (TTree*)input->Get("crumbs/SliceTree");
  TTree *background     = (TTree*)input->Get("crumbs/SliceTree");

  TString outDirName("CRUMBSDataset_CCNuE");
  gSystem->Exec("mkdir " + outDirName);
  TMVA::DataLoader *dataloader=new TMVA::DataLoader(outDirName);

  TString outfileName( outDirName + "/CRUMBSTMVA_CCNuE.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "CrumbsTMVAClassification", outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

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

  dataloader->AddVariable("pds_FMTotalScore","Total FM Score","",'F');
  dataloader->AddVariable("pds_FMPE","nPE in flash","",'F');
  dataloader->AddVariable("pds_FMTime","FM Time","#mu s",'F');

  dataloader->AddVariable("crt_TrackScore","CRT Track Match Score","",'F');
  dataloader->AddVariable("crt_HitScore","CRT Hit Match Score","",'F');
  dataloader->AddVariable("crt_TrackTime","CRT Track Match Time","#mu s",'F');
  dataloader->AddVariable("crt_HitTime","CRT Hit Match Time","#mu s",'F');

  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  dataloader->AddSignalTree    (signalTree, signalWeight);
  dataloader->AddBackgroundTree(background, backgroundWeight);

  TCut mycuts = "strstr(matchedType,\"Nu\") && !strstr(matchedType,\"DirtNu\") && matchedPurity > 0.8 && matchedCompleteness > 0.8 && ccnc == 0 && abs(nutype) == 12";
  TCut mycutb = "!strstr(matchedType,\"Nu\")";

  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "SplitMode=random:!V" );
  
  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
		       "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

  factory->TrainAllMethods();
  
  factory->TestAllMethods();
  
  factory->EvaluateAllMethods();
  
  outputFile->Close();
  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;
  delete dataloader;

  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

  return 0;
}
