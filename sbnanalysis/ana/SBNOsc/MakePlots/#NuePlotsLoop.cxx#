//Really Bad Macro for making the selection plots. 

//C++ Includes 
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

//Root Includes 
#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h" 
#include "TPad.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFileCollection.h"
#include "TMath.h"

#include "/sbnd/app/users/dbarker/larsoft_sbncode/sbncode/srcs/sbncode/sbnanalysis/core/Event.hh"
#include "/sbnd/app/users/dbarker/larsoft_sbncode/sbncode/srcs/sbncode/sbnanalysis/core/SubRun.hh"

void NuePlotsLoop(){

  //g++ `root-config --cflags` -L$SBN_LIB_DIR /sbnd/app/users/dbarker/larsoft_sbncode/sbncode/srcs/sbncode/sbnanalysis/ana/SBNOsc/MakePlots/NuePlotsLoop.cxx -lsbnanalysis_Event `root-config --libs` -o NuePlotsLoop

  bool fOnlySamples = true;

  //  int nbinsvals = 6;
  //  int nbinsvals = 9;
  int nbinsvals = 11;

  //Stack to go in
  Double_t nue_bins[nbinsvals+1] = {0.200,0.350,0.500,0.650,0.800,0.950,1.100,1.300,1.500,1.750,2.000,3.000};
  //  Double_t nue_bins[nbinsvals+1] = //{0.200,0.450,0.650,0.750,1050,1.5,1.4,1.7,2.2,3};
  //     {0.1,0.450,0.650,0.900,1.250,2.000,5.000};

  //  Double_t nue_bins[nbinsvals+1] = {0.225,5};



  // nue_bins.push_back(0.100);
  // nue_bins.push_back(0.200);
  // nue_bins.push_back(0.350);
  // nue_bins.push_back(0.500);
  // nue_bins.push_back(0.650);
  // nue_bins.push_back(0.800);
  // nue_bins.push_back(0.950);
  // nue_bins.push_back(1.100);
  // nue_bins.push_back(1.300);
  // nue_bins.push_back(1.500);
  // nue_bins.push_back(1.750);
  // nue_bins.push_back(2.000);
  // nue_bins.push_back(3.000);

  TH1D* hbin =  new TH1D("hbin", "hbin", nbinsvals,nue_bins);
  THStack*  VisibleEnergy_Selection_Stack = new THStack("VisibleEnergy_Selection_Stack","VisibleEnergy_Selection_Stack");
  VisibleEnergy_Selection_Stack->SetHistogram(hbin);


  //get the text file with the sbn files in 
  std::string inputfile_name;
  std::cout << "Please enter the name of the list of the input files: ";
  std::cin >> inputfile_name;
  const char* inputfilelist = inputfile_name.c_str();
  
  //Flavours of events 
  std::vector<std::string> HistTypes = {"NuMu","InNuE","OscNuE","NCInNuE","NCOscNuE","NCNuMu","DirtNuMu","DirtInNuE","DirtOscNuE","DirtNCInNuE","DirtNCOscNuE","DirtNCNuMu","Cosmic"};

  //Merge NC and Dirt Events.
  std::vector<std::string> HistTypesFinal = {"InNuE F","NC F","NuMu F","Dirt F", "Cosmic F","OscNuE F"};

  //  std::vector<int>         Colours   = {419,418,792,603};
  //Sort out the colours 
  std::map<std::string,int> Colours;
  Colours["OscNuE F"] = 419;
  Colours["InNuE F"] = 418;
  Colours["NuMu F"] = 603;
  Colours["NC F"] = 792;
  Colours["Dirt F"] = 920;
  Colours["Cosmic F"] = 2;

  //Histrograms for the prdocutions
  std::map<std::string,TH1D*> VisibleEnergyHist;
  for(int i=0; i<HistTypes.size(); ++i){
    const char* hist_name = HistTypes[i].c_str();
    VisibleEnergyHist[HistTypes[i]] = new TH1D(hist_name, hist_name, nbinsvals,nue_bins);
  }

  std::map<std::string,TH1D*> VisibleEnergyHistFinal;
  for(int i=0; i<HistTypesFinal.size(); ++i){
    const char* hist_name = HistTypesFinal[i].c_str();
    VisibleEnergyHistFinal[HistTypesFinal[i]] = new TH1D(hist_name, hist_name, nbinsvals,nue_bins);
  }

  TH1D* VisibleEnergyHistTotal = new TH1D("totalhist", "totalhist", nbinsvals,nue_bins);

  float ProposalPOT = 6.6e20;

  float POTWeight;
  float POTTotal = 0;
  std::map<std::string,float> POTMap;
  std::string inputfile;
  std::string line;


  //Loop over the text files and make a root file.              
  std::ifstream txtfile_stream;
  txtfile_stream.open(inputfilelist);
  if (txtfile_stream.is_open()){

    //First Get the POT norm for each element.
    while (txtfile_stream.good()){

      std::getline(txtfile_stream, line);

      if(line.size() == 0){
        continue;
      }

      std::cout << "New file with name: " << line << std::endl;

      std::size_t pos = line.find(" ");

      inputfile = line.substr(0,pos);

      std::string HistType = " ";

      std::size_t isuboone = line.find("uboone");
      if(isuboone !=std::string::npos){
	ProposalPOT = 1.32e21;
      }
      isuboone = line.find("Uboone");
      if(isuboone !=std::string::npos){
        ProposalPOT = 1.32e21;
      }


      //Find which sample it came from. we did a dedicated intrinsic, osc and muon run. 
      std::size_t found_intrinsic = line.find("nue_intrinsic");
      if(found_intrinsic!=std::string::npos){
    	HistType = "InNuE";
      }
      found_intrinsic = line.find("Int");
      if(found_intrinsic!=std::string::npos){
	std::cout << "is nue" << std::endl;
        HistType = "InNuE";
      }
      found_intrinsic = line.find("IntBackUp");
      if(found_intrinsic!=std::string::npos){
        HistType = "InNuE";
      }

      //idiot 
      std::size_t found_osc = line.find("nue_osc");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
      } 
      //idiot 
      found_osc = line.find("_Osc.");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
      } 
      found_osc = line.find("_Osc__");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
      } 

      found_osc = line.find("OscBackup");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
      } 


      std::size_t found_numu = line.find("nu/");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
      }
      found_numu = line.find("_Numu");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
      }
      found_numu = line.find("nu_");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
      }
      found_numu = line.find("Numu__Nue_");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
      }

      std::size_t found_cosdirt = line.find("nu_cosdirt");
      if(found_cosdirt!=std::string::npos){
	HistType = "CosDirt";
      }

      found_cosdirt = line.find("CosDirt.");
      if(found_cosdirt!=std::string::npos){
        HistType = "CosDirt";
      }

      found_cosdirt = line.find("_Dirt.");
      if(found_cosdirt!=std::string::npos){
        HistType = "Dirt";
      }

	  
      
      if(HistType == " "){
	HistType = "NuMu";
	std::cout << "It went wrong" << std::endl;
      }

      
    
      //Get the inputfile                                                                                       
      const char* InputFileName = inputfile.c_str();
      TFile *InputFile = TFile::Open(InputFileName);

      if(!InputFile->IsOpen()){ std::cout << "File: " << InputFileName << " failed to open. exiting" << std::endl;return;}

      //Get the POT                        

      TTree *subruntree;

      gDirectory->GetObject("sbnsubrun",subruntree);
      gDirectory->ls();
      double totgoodpot;
      float totalpot = 0;

      SubRun *subruns  = 0;
      subruntree->SetBranchAddress("subruns",&subruns);

      std::string exp;
      std::size_t issbnd = line.find("sbnd");
      std::size_t isSBNDBIG = line.find("SBND");
      if(issbnd!=std::string::npos){
	exp = "sbnd";
      }
      else if(isSBNDBIG!=std::string::npos){
	exp = "sbnd";
      }

      std::size_t isuboonebig = line.find("Uboone");
      if(isuboone!=std::string::npos){
	exp = "uboone";
      }
      else if(isuboonebig!=std::string::npos){
	exp = "uboone";
      }


      std::size_t isicarus = line.find("icarus");
      std::size_t isicarusbig = line.find("Icarus");
      if(isicarus!=std::string::npos){
	exp = "icarus"; 
      }
      else if(isicarusbig!=std::string::npos){
	exp = "icarus"; 
      }


      std::cout << "HistType: " << HistType << std::endl;

      Long64_t nentries = subruntree->GetEntries();
      for (Long64_t i=0;i<nentries;i++) {

	//subruntree->LoadTree(i);
        subruntree->GetEntry(i);
	//	HistType = "OscNuE";
	//std::cout << "totgoodpot: " << subruns->totgoodpot << " i: " << i  << std::endl;
	//	std::cout << "HistType: " << HistType << std::endl;
	double POT =  subruns->totgoodpot;
	
	if(HistType == "NuMu" && (exp == "uboone" ||  exp== "icarus")){
	  POT /= 10;
	}

        POTMap[HistType] += POT;
      }
      InputFile->Close();
    }
  }

  std::map<std::string,float> POTFlavourMap; 
  //Now we need to add up the POT for each flavour
  for(int i=0; i<HistTypes.size(); ++i){

    std::cout << "Pot for Hist type: " << HistTypes[i] << " is " <<  POTMap[HistTypes[i]] << std::endl;
    
    if(HistTypes[i] == "OscNuE"){
      //Only comes from the osc nue sample son
      POTFlavourMap["OscNuE"] += POTMap["OscNuE"];
    }
    
    if(HistTypes[i] == "InNuE"){
      //Comes from the intrinsic sample and numu sample.
      POTFlavourMap["InNuE"] += POTMap["InNuE"];
      //      POTFlavourMap["InNuE"] = 4.788e+19 + 5.12255e+21;
      //      POTFlavourMap["InNuE"]  = 983150000000000065536.000000 + 1.23124e+23;
      //      POTFlavourMap["InNuE"] = 118884300000000000000.000000 + 6.48597e+21;
    }
    
    if(HistTypes[i] == "NuMu"){
      //Comes intrinic
      POTFlavourMap["NuMu"] += POTMap["NuMu"];
      //POTFlavourMap["NCNuMu"] =  4.788e+19;
      //POTFlavourMap["NuMu"] += 983150000000000065536.000000;
      //      POTFlavourMap["NuMu"] = 118884300000000000000.000000;

    }
    
    if(HistTypes[i] == "NCInNuE"){
      POTFlavourMap["NCInNuE"] += POTMap["InNuE"];
      //POTFlavourMap["NCInNuE"] = 4.788e+19 + 5.12255e+21;
      // POTFlavourMap["NCInNuE"] = 983150000000000065536.000000 + 1.23124e+23;
      //      POTFlavourMap["NCInNuE"] =  118884300000000000000.000000 + 6.48597e+21

    }

    if(HistTypes[i] == "NCOscNuE"){
      POTFlavourMap["NCOscNuE"] += POTMap["OscNuE"];
    }
    
    if(HistTypes[i] == "NCNuMu"){
      POTFlavourMap["NCNuMu"] += POTMap["NuMu"];
      //      POTFlavourMap["NCNuMu"] =  4.788e+19; 
      //      POTFlavourMap["NuMu"] =983150000000000065536.000000;
	    //POTFlavourMap["NuMu"] = 118884300000000000000.000000;

    }

    if(HistTypes[i] =="DirtNuMu"){
      POTFlavourMap["DirtNuMu"] += POTMap["Dirt"];
    }

    if(HistTypes[i] == "DirtInNuE"){
      //POTFlavourMap["DirtInNuE"] += POTMap["InNuE"];
      //      POTFlavourMap["DirtInNuE"] += POTMap["NuMu"];
      POTFlavourMap["DirtInNuE"] += POTMap["Dirt"];
    }

    if(HistTypes[i] == "DirtOscNuE"){
      POTFlavourMap["DirtOscNuE"] += POTMap["Dirt"];
    }
     
    if(HistTypes[i] == "DirtNCInNuE"){
      //POTFlavourMap["DirtNCInNuE"] += POTMap["InNuE"]; 
      //POTFlavourMap["DirtNCInNuE"] += POTMap["NuMu"];
      POTFlavourMap["DirtNCInNuE"] += POTMap["Dirt"];
    }

    if(HistTypes[i] ==  "DirtNCOscNuE"){
      POTFlavourMap["DirtNCOscNuE"] += POTMap["Dirt"];
    }

    if(HistTypes[i] ==  "DirtNCNuMu"){
      POTFlavourMap["DirtNCNuMu"] += POTMap["Dirt"];
    }

    if(HistTypes[i] == "Cosmic"){
      POTFlavourMap["Cosmic"] += POTMap["CosDirt"];
    }

  }

  for(int i=0; i<HistTypes.size(); ++i){
    std::cout << "POT for type: " << HistTypes[i] << " is: " << POTFlavourMap[HistTypes[i]] << std::endl;;
  }

  
  txtfile_stream.close();
  txtfile_stream.open(inputfilelist);

  if (txtfile_stream.is_open()){

    while (txtfile_stream.good()){
    
      std::getline(txtfile_stream, line);

      if(line.size() == 0){
	continue;
      }

      std::string HistType = " ";
      double POT = 0;

      //Find which sample it came from. we did a dedicated intrinsic, osc and muon run. 
      std::size_t found_intrinsic = line.find("nue_intrinsic");
      if(found_intrinsic!=std::string::npos){
    	HistType = "InNuE";
	POT = POTFlavourMap[HistType];
      }
      found_intrinsic = line.find("Int");
      if(found_intrinsic!=std::string::npos){
        HistType = "InNuE";
	POT = POTFlavourMap[HistType];
      }
      found_intrinsic = line.find("IntBackUp");
      if(found_intrinsic!=std::string::npos){
        HistType = "InNuE";
	POT = POTFlavourMap[HistType];
      }

      //idiot 
      std::size_t found_osc = line.find("nue_osc");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
	POT = POTFlavourMap[HistType];
      } 
      //idiot 
      found_osc = line.find("_Osc.");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
	POT = POTFlavourMap[HistType];
      } 
      found_osc = line.find("_Osc__");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
	POT = POTFlavourMap[HistType];
      } 

      found_osc = line.find("OscBackup");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
	POT = POTFlavourMap[HistType];
      } 

      std::size_t found_numu = line.find("nu/");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
	POT = POTFlavourMap[HistType];
      }
      found_numu = line.find("Numu.");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
	POT = POTFlavourMap[HistType];
      }
      found_numu = line.find("nu_");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
	POT = POTFlavourMap[HistType];
      }
      found_numu = line.find("Numu__Nue_");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
      }


      std::size_t found_cosdirt_files = line.find("nu_cosdirt");
      if(found_cosdirt_files!=std::string::npos){
	HistType = "Cosmic";
	POT = POTFlavourMap[HistType];
      }

      std::size_t found_cosdirt_file = line.find("CosDirt.");
      if(found_cosdirt_file!=std::string::npos){
        HistType = "Cosmic";
	POT = POTFlavourMap[HistType];
      }

      std::size_t found_dirt_file = line.find("_Dirt.");
      if(found_dirt_file!=std::string::npos){
        HistType = "Dirt";
	POT = POTFlavourMap[HistType];
      }

      std::cout << "HistType: " << HistType << std::endl;

      if(HistType == " " ){
	std::cout << " Something went wrong in sorting the POT: " << HistType << std::endl;
	//	continue;
      }

      //      double POTWeight = ProposalPOT/POT;
      POTWeight = 1;
      std::size_t pos = line.find(" ");    
      
      inputfile = line.substr(0,pos);

      //Get the inputfile
      const char* InputFileName = inputfile.c_str();
      TFile *InputFile = TFile::Open(InputFileName);

      
      std::cout << "file: " << inputfile << std::endl;

      TTree *sbnreco;
      gDirectory->GetObject("sbnreco",sbnreco);
      event::RecoEvent *reco_events = 0;
      sbnreco->SetBranchAddress("reco_events",&reco_events);

      double totgoodpot;
      float totalpot = 0;

      //      double L = 100.*5076142.;
      double L = 600.*5076142.;
      std::size_t issbnd = line.find("sbnd");
      std::size_t isSBNDBIG = line.find("SBND");
      if(issbnd!=std::string::npos){
	L = 100.*5076142.;
	//L = 90.*5076142.; 
      }
      else if(isSBNDBIG!=std::string::npos){
	L = 100.*5076142.;
		//	L = 90.*5076142.;
      }

      std::size_t isuboone = line.find("uboone");
      std::size_t isuboonebig = line.find("Uboone");
      if(isuboone!=std::string::npos){
	L = 470.*5076142.;
      }
      else if(isuboonebig!=std::string::npos){
	L = 470.*5076142.;
      }


      std::size_t isicarus = line.find("icarus");
      std::size_t isicarusbig = line.find("Icarus");
      if(isicarus!=std::string::npos){
	L = 600.*5076142.;
      }
      else if(isicarusbig!=std::string::npos){
	L = 600.*5076142.;
      }

      Long64_t nentries = sbnreco->GetEntries();
      std::cout << "nentries: " << nentries << std::endl;
      for (Long64_t i=0;i<nentries;i++) {


	sbnreco->LoadTree(i);
        sbnreco->GetEntry(i);

	event::RecoInteraction * reco = &(reco_events->reco);
      	std::vector<event::Interaction> * truth = &(reco_events->truth);

	if(reco->reco_energy < 0){continue;}
	//	if(truth->size() != 1){VisibleEnergyHist["Cosmic"]->Fill(reco->reco_energy, reco->weight*POTWeight); std::cout << "Why is this not 1? " << reco->weight << std::endl; continue;}
	//	if(truth->at(0).neutrino.pdg == 0){VisibleEnergyHist["Cosmic"]->Fill(reco->reco_energy, reco->weight*POTWeight); std::cout << "Why is this not 1? " << reco->weight << std::endl; continue;}

	if(truth->at(0).neutrino.genie_intcode == 5){continue;}
					  
	//Osc Weight 
	double Theta = 0.001;// 0.013; 
	double dm    = 1.32; //0.43; 
	double osc_w = 1;
	if(HistType == "OscNuE"){
	  osc_w = Theta*TMath::Sin(dm*L/(4*truth->at(0).neutrino.energy*1e9))*TMath::Sin(dm*L/(4*truth->at(0).neutrino.energy*1e9)); 
	}
	
	//std::cout << " energy: " << reco->reco_energy << " weight: " << reco->weight << " HistType: " << HistType  <<  " POTWeight: " << POTWeight  << " oscw: " <<  osc_w<< std::endl;
	
	//	bool  wasCosmic = reco->wasCosmic;
	//	bool  wasDirt   = reco->wasDirt;

	bool wasCosmic = false;
	bool wasDirt = false;

	// //Are we cosmic
	if(HistType == "Cosmic"){
	  VisibleEnergyHist["Cosmic"]->Fill(reco->reco_energy, reco->weight*POTWeight);
	  continue;
	}
	
	//Are we dirt backround
	if(HistType == "Dirt"){
	  
	  
	  if(truth->at(0).neutrino.isnc == false){

	    //Are we oscillated 
	    if(truth->at(0).neutrino.initpdg != truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){
	      VisibleEnergyHist["DirtOscNuE"]->Fill(reco->reco_energy,reco->weight*osc_w*POTWeight);
	      continue; 
	    }
	      
	    //Are we intrinsic?
	    if(truth->at(0).neutrino.initpdg == truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){
	      VisibleEnergyHist["DirtInNuE"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	      
	    //Are we charged current background 
	    if(TMath::Abs(truth->at(0).neutrino.pdg == 14) && truth->at(0).neutrino.iscc == true){
	      VisibleEnergyHist["DirtNCNuMu"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	  }
	    
	  //Are we NC background 
	  if(truth->at(0).neutrino.isnc == true){
	    //Are we oscillated 
	    if(truth->at(0).neutrino.initpdg != truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){
	      VisibleEnergyHist["DirtNCOscNuE"]->Fill(reco->reco_energy,reco->weight*osc_w*POTWeight);
	      continue; 
	    }
	    //Are we intrinsic?
	    if(truth->at(0).neutrino.initpdg == truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){
	      VisibleEnergyHist["DirtNCInNuE"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	    //Are we charged current background 
	    if(TMath::Abs(truth->at(0).neutrino.pdg) == 14 && truth->at(0).neutrino.iscc == true){
	      VisibleEnergyHist["DirtNCNuMu"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	      }
	  }
	  std::cout << "we have got to the end" << std::endl;
	  continue;
	}
	  
	  if(truth->at(0).neutrino.isnc == false){
	    //Are we oscillated?
	    if(truth->at(0).neutrino.initpdg != truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){	    
	      if(HistType == "OscNuE" || fOnlySamples)VisibleEnergyHist["OscNuE"]->Fill(reco->reco_energy,reco->weight*osc_w*POTWeight);
	      continue;
	    } 
	    
	    
	    //Are we intrinsic?
	    if(truth->at(0).neutrino.initpdg == truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){	
	      //	      std::cout << "filling" << std::endl;
	      if(HistType == "InNuE" || fOnlySamples) VisibleEnergyHist["InNuE"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	    
	    //Are we charged current background 
	    if(TMath::Abs(truth->at(0).neutrino.pdg) == 14 && truth->at(0).neutrino.iscc == true){
	      if(HistType == "NuMu" || fOnlySamples) VisibleEnergyHist["NuMu"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	  }
      
	  //Are we NC background 
	  if(truth->at(0).neutrino.isnc == true){
	    
	    //	    std::cout << " energy: " << reco->reco_energy << " weight: " << reco->weight << " mode: " << truth->at(0).neutrino.genie_intcode<< std::endl;
	    //	    if(reco->weight > 0.1){std::cout << "WARNING DODGY WEIGHT" << std::endl;}
	    //	    if(reco->weight > 0.1){continue;}

	    //Are we oscillated?
	    if(truth->at(0).neutrino.initpdg != truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){
	      //VisibleEnergyHist["NCOscNuE"]->Fill(reco->reco_energy,reco->weight*osc_w*POTWeight);
	      continue;
	    } 
	    
	    //Are we intrinsic?
	    if(truth->at(0).neutrino.initpdg == truth->at(0).neutrino.pdg && TMath::Abs(truth->at(0).neutrino.pdg) == 12){
	      if(HistType == "InNuE" || fOnlySamples) VisibleEnergyHist["NCInNuE"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	  
	    //Are we charged current background 
	    if(TMath::Abs(truth->at(0).neutrino.pdg) == 14){
	      if(HistType == "NuMu" || fOnlySamples) VisibleEnergyHist["NCNuMu"]->Fill(reco->reco_energy,reco->weight*POTWeight);
	      continue;
	    }
	  }
	  std::cout << "end" <<  truth->at(0).neutrino.pdg << " " << truth->at(0).neutrino.iscc <<  std::endl;
      }
    
      InputFile->Close();
    }
    
    std::cout << "total POT Wieght: " << POTTotal << std::endl;

  }
  else{
    std::cerr << "Could not open the input text file" << std::endl;
    return;
  }

  //TFile to hold the  selection plots.
  TFile *Figures = new TFile("Figures.root","RECREATE");


  for(int i=0; i<HistTypes.size(); ++i){

    float POTWeight = 1; 
    if(POTFlavourMap[HistTypes[i]] != 0){
      POTWeight = ProposalPOT/POTFlavourMap[HistTypes[i]];
    }

       std::cout << "POTWeight: " << POTWeight << " for sample:" << HistTypes[i]  << " proposal POT: " << ProposalPOT<<  std::endl;

       VisibleEnergyHist[HistTypes[i]]->Scale(POTWeight,"width");

    VisibleEnergyHistTotal->Add(VisibleEnergyHist[HistTypes[i]]);

    VisibleEnergyHist[HistTypes[i]]->Write();

    float Totalint = 0;
    for(int bin=0; bin< VisibleEnergyHist[HistTypes[i]]->GetNbinsX(); ++bin){                                    
      //      VisibleEnergyHistTotal->Fill(VisibleEnergyHistTotal->GetXaxis()->GetBinCenter(bin+1),POTWeight*VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin));

      //      Totalint += VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin);
      std::cout << " VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin): " << VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin) << " POTWeight: " << POTWeight << " width: " <<  VisibleEnergyHist[HistTypes[i]]->GetBinWidth(1+bin)  << std::endl;
      //      VisibleEnergyHist[HistTypes[i]]->SetBinContent(1+bin, POTWeight*VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin) / VisibleEnergyHist[HistTypes[i]]->GetBinWidth(1+bin));
      //      float Error =  VisibleEnergyHist[HistTypes[i]]->GetBinError(1+bin)*POTWeight/VisibleEnergyHist[HistTypes[i]]->GetBinWidth(1+bin);
      //      VisibleEnergyHist[HistTypes[i]]->SetBinError(1+bin,Error);

	//      if(HistTypes[i] != "OscNuE"){ VisibleEnergyHist[HistTypes[i]]->SetBinError(1+bin,0);
      //  }

      std::cout << " Bin Entry is: " << VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin) << " +- " << VisibleEnergyHist[HistTypes[i]]->GetBinError(1+bin) << std::endl;
      
    }

    for(int bin=0; bin< VisibleEnergyHist[HistTypes[i]]->GetNbinsX(); ++bin){
      std::cout << "Total Bin content: " << VisibleEnergyHistTotal->GetBinContent(1+bin) << " total bin error: " << VisibleEnergyHistTotal->GetBinError(1+bin) << std::endl;
    }

    double error = 0;
    //    Totalint = VisibleEnergyHist[HistTypes[i]]->IntegralAndError((int) 0,(int) VisibleEnergyHist[HistTypes[i]]->GetNbinsX(),error,"width");

    std::cout << "total flavour enteries is: " <<  Totalint << " +- " << error << std::endl;
    std::cout << " total: " << VisibleEnergyHist[HistTypes[i]]->Integral() << std::endl;
    // VisibleEnergyHist[HistTypes[i]]->SetFillColor(Colours[HistTypes[i]]);
    // VisibleEnergyHist[HistTypes[i]]->SetMarkerStyle(21);
    // VisibleEnergyHist[HistTypes[i]]->SetMarkerColor(colours[HistTypes[i]]);
    // VisibleEnergyHist[HistTypes[i]]->SetLineColor(Colours[HistTypes[i]]);

    //Merge NC and Dirt Events.                                                                                                                                                     
    //Add to the final historams
    if(HistTypes[i] == "OscNuE"){
      VisibleEnergyHistFinal["OscNuE F"]->Add(VisibleEnergyHist[HistTypes[i]]); 
    }

    if(HistTypes[i] == "NuMu"){
      VisibleEnergyHistFinal["NuMu F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "InNuE"){
      VisibleEnergyHistFinal["InNuE F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "NCInNuE" || HistTypes[i] == "NCNuMu" || HistTypes[i] == "NCOscNuE"){
      VisibleEnergyHistFinal["NC F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "DirtInNuE" || HistTypes[i] == "DirtNuMu" || HistTypes[i] == "DirtNCInNuE" || HistTypes[i] == "DirtNCNuMu"){
      VisibleEnergyHistFinal["Dirt F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "Cosmic"){
      VisibleEnergyHistFinal["Cosmic F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

  //    VisibleEnergy_Selection_Stack->Add(VisibleEnergyHist[HistTypes[i]]);
  }

  //Add to the stack 
  for(int i=0;i< HistTypesFinal.size(); ++i){
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetFillColor(Colours[HistTypesFinal[i]]);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetMarkerStyle(21);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetMarkerColor(kBlue);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetLineColor(Colours[HistTypesFinal[i]]);
    VisibleEnergy_Selection_Stack->Add(VisibleEnergyHistFinal[HistTypesFinal[i]]);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->Write();

    if(HistTypesFinal[i] == "OscNuE F"){
      VisibleEnergyHistFinal[HistTypesFinal[i]]->SetFillColor(0);
      VisibleEnergyHistFinal[HistTypesFinal[i]]->SetLineColor(kBlack);
      VisibleEnergyHistFinal[HistTypesFinal[i]]->SetLineWidth(5);
    }
  }

  VisibleEnergyHistTotal->Write();

  VisibleEnergy_Selection_Stack->Write();

  TCanvas *StackCanvas = new TCanvas("StackCanvas","StackCanvas",10,10,700,900);
  VisibleEnergy_Selection_Stack->Draw("eHIST");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  StackCanvas->Write();

    for(int bin=0; bin< VisibleEnergyHistTotal->GetNbinsX(); ++bin){
      std::cout << "Total Bin content: " << VisibleEnergyHistTotal->GetBinContent(1+bin) << " total bin error: " << VisibleEnergyHistTotal->GetBinError(1+bin)  << " percentage: " << (VisibleEnergyHistTotal->GetBinError(1+bin)/VisibleEnergyHistTotal->GetBinContent(1+bin))*100 << std::endl; ;
    }


}

# ifndef __CINT__
int main() {
  NuePlotsLoop();
  return 0;
}

# endif

