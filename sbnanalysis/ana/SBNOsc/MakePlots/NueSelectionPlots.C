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

#include "/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/core/Event.hh"
#include "/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/core/SubRun.hh"

void NueSelectionPlots(){

  //Stack to go in
  Double_t nue_bins[12] = {0.200,0.350,0.500,0.650,0.800,0.950,1.100,1.300,1.500,1.750,2.000,3.000};

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

  TH1D* hbin =  new TH1D("hbin", "hbin", 11,nue_bins);
  THStack*  VisibleEnergy_Selection_Stack = new THStack("VisibleEnergy_Selection_Stack","VisibleEnergy_Selection_Stack");
  VisibleEnergy_Selection_Stack->SetHistogram(hbin);


  //get the text file with the sbn files in 
  std::string inputfile_name;
  cout << "Please enter the name of the list of the input files: ";
  cin >> inputfile_name;
  const char* inputfilelist = inputfile_name.c_str();
  
  //Flavours of events 
  std::vector<std::string> HistTypes = {"NuMu","InNuE","OscNuE","NCInNuE","NCOscNuE","NCNuMu","DirtNuMu","DirtInNuE","DirtOscNuE","DirtNCInNuE","DirtNCOscNuE","DirtNCNuMu","Cosmic"};

  //Merge NC and Dirt Events.
  std::vector<std::string> HistTypesFinal = {"InNuE F","NC F","NuMu F"};

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
    VisibleEnergyHist[HistTypes[i]] = new TH1D(hist_name, hist_name, 11,nue_bins);
  }

  std::map<std::string,TH1D*> VisibleEnergyHistFinal;
  for(int i=0; i<HistTypesFinal.size(); ++i){
    const char* hist_name = HistTypesFinal[i].c_str();
    VisibleEnergyHistFinal[HistTypesFinal[i]] = new TH1D(hist_name, hist_name, 11,nue_bins);
  }

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

      //Find which sample it came from. we did a dedicated intrinsic, osc and muon run. 
      std::size_t found_intrinsic = line.find("nue_intrinsic/");
      if(found_intrinsic!=std::string::npos){
    	HistType = "InNuE";
      }

      //idiot 
      std::size_t found_osc = line.find("nue_osc/");
      if(found_osc!=std::string::npos){
        HistType = "OscNuE";
      } 

      std::size_t found_numu = line.find("nu/");
      if(found_numu!=std::string::npos){
        HistType = "NuMu";
      }

      std::size_t found_cosdirt = line.find("nu_cosmic_dirt/");
      if(found_cosdirt!=std::string::npos){
	HistType = "CosDirt";
      }
	  
      
      if(HistType == " "){
	std::cout << "It went wrong" << std::endl;
      }

      //Get the inputfile                                                                                       
      const char* InputFileName = inputfile.c_str();
      TFile *InputFile = new TFile(InputFileName);

      //Get the POT                                                                                                         
      TTree *subruntree = (TTree*)InputFile->Get("sbnsubrun");

      double totgoodpot;
      float totalpot = 0;

      SubRun *subruns = new SubRun;

      subruntree->SetBranchAddress("subruns",&subruns);



      Long64_t nentries = subruntree->GetEntries();
      for (Long64_t i=0;i<nentries;i++) {
        subruntree->GetEntry(i);
	//	HistType = "OscNuE";
	//	std::cout << "totgoodpot: " << subruns->totgoodpot << " i: " << i  << std::endl;
	//	std::cout << "HistType: " << HistType << std::endl;
        POTMap[HistType] += subruns->totgoodpot;
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
      POTFlavourMap["InNuE"] += POTMap["NuMu"];
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
      POTFlavourMap["NCInNuE"] += POTMap["NuMu"];
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
      POTFlavourMap["DirtNuMu"] += POTMap["NuMu"];
      POTFlavourMap["DirtNuMu"] += POTMap["CosDirt"];
    }

    if(HistTypes[i] == "DirtInNuE"){
      POTFlavourMap["DirtInNuE"] += POTMap["InNuE"];
      POTFlavourMap["DirtInNuE"] += POTMap["NuMu"];
      POTFlavourMap["DirtInNuE"] += POTMap["CosDirt"];
    }

    if(HistTypes[i] == "DirtOscNuE"){
      POTFlavourMap["DirtOscNuE"] += POTMap["OscNuE"];
    }
     
    if(HistTypes[i] == "DirtNCInNuE"){
      POTFlavourMap["DirtNCInNuE"] += POTMap["InNuE"]; 
      POTFlavourMap["DirtNCInNuE"] += POTMap["NuMu"];
      POTFlavourMap["DirtNCInNuE"] += POTMap["CosDirt"];
    }

    if(HistTypes[i] ==  "DirtNCOscNuE"){
      POTFlavourMap["DirtNCOscNuE"] += POTMap["OscNuE"];
    }

    if(HistTypes[i] ==  "DirtNCNuMu"){
      POTFlavourMap["DirtNCNuMu"] += POTMap["NuMu"];
      POTFlavourMap["DirtNCNuMu"] += POTMap["CosDirt"];
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


      std::size_t pos = line.find(" ");    
      
      inputfile = line.substr(0,pos);
      //      POTWeight = stof(line.substr(pos));
      //std::cout << "POTWeight: " << POTWeight <<std::endl;
            
      //Get the inputfile
      const char* InputFileName = inputfile.c_str();
      TFile *InputFile = new TFile(InputFileName);

      // //Get the POT
      // TTree *subruntree = (TTree*)InputFile->Get("sbnsubrun");

      // double totgoodpot;
      // float totalpot = 0;

      // SubRun *subruns = new SubRun;

      // subruntree->SetBranchAddress("subruns",&subruns);

      // Long64_t nentries = subruntree->GetEntries();
      // for (Long64_t i=0;i<nentries;i++) {
      // 	subruntree->GetEntry(i);
      // 	std::cout << "totgoodpot: " << subruns->totgoodpot << " i: " << i  << std::endl;
      // 	totalpot += subruns->totgoodpot;
      // }

      // std::cout << "totalpot: " << totalpot << std::endl;
      
      //move to the Hist folder
      InputFile->cd("Histograms");
 
      for(int i=0; i<HistTypes.size(); ++i){
	
	//	std::cout << "New hist type" << std::endl;
		
	//float POTWeight = (6.6e20)/POTFlavourMap[HistTypes[i]];
	//std::cout << "POTWeight: " << POTWeight <<std::endl;

	std::string  VisibleEnergy_Selection_String  = HistTypes[i] + " VisibleEnergy_Selection";
	const char*  VisibleEnergy_Selection_Name    = VisibleEnergy_Selection_String.c_str();
	
	//Get the histogram
	if(gDirectory->Get(VisibleEnergy_Selection_Name) == NULL){continue;}
	TH1D* VisibleEnergyHist1 = (TH1D*)(gDirectory->Get(VisibleEnergy_Selection_Name))->Clone();
	TH1D* VisibleEnergyHist2 = (TH1D*)VisibleEnergyHist1->Rebin(11,VisibleEnergy_Selection_Name,nue_bins);

	//Account for /Gev
	// for(int bin=0; bin< VisibleEnergyHist2->GetNbinsX(); ++bin){
	//   //	  std::cout << "VisibleEnergyHist2->GetBinWidth(1+bin): " << VisibleEnergyHist2->GetBinWidth(1+bin) << std::endl;
	//   VisibleEnergyHist2->SetBinContent(1+bin, POTWeight*VisibleEnergyHist2->GetBinContent(1+bin) / VisibleEnergyHist2->GetBinWidth(1+bin));
	// }

	if(VisibleEnergyHist2 == NULL){
	  std::cerr << "Invalid Histogram initalisation" << std::endl;
	  return;
	}

	// VisibleEnergyHist2->SetFillColor(Colours[HistTypes[i]]);
	// VisibleEnergyHist2->SetMarkerStyle(21);
	// VisibleEnergyHist2->SetMarkerColor(Colours[HistTypes[i]]);
	// VisibleEnergyHist2->SetLineColor(Colours[HistTypes[i]]);

	VisibleEnergyHist[HistTypes[i]]->Add(VisibleEnergyHist2);
      }
      InputFile->Close();
    }
    
    std::cout << "total POT Wieght: " << POTTotal << std::endl;

  }
  else{
    std::cerr << "Could not open the input text file" << std::endl;
    return;
  }

  for(int i=0; i<HistTypes.size(); ++i){

    float POTWeight = 1; 
    if(POTFlavourMap[HistTypes[i]] != 0){
      POTWeight = ProposalPOT/POTFlavourMap[HistTypes[i]];
    }

    std::cout << "POTWeight: " << POTWeight << " for sample:" << HistTypes[i]  << " proposal POT: " << ProposalPOT<<  std::endl;
    
    for(int bin=0; bin< VisibleEnergyHist[HistTypes[i]]->GetNbinsX(); ++bin){                                    
      VisibleEnergyHist[HistTypes[i]]->SetBinContent(1+bin, POTWeight*VisibleEnergyHist[HistTypes[i]]->GetBinContent(1+bin) / VisibleEnergyHist[HistTypes[i]]->GetBinWidth(1+bin));
    }

    // VisibleEnergyHist[HistTypes[i]]->SetFillColor(Colours[HistTypes[i]]);
    // VisibleEnergyHist[HistTypes[i]]->SetMarkerStyle(21);
    // VisibleEnergyHist[HistTypes[i]]->SetMarkerColor(colours[HistTypes[i]]);
    // VisibleEnergyHist[HistTypes[i]]->SetLineColor(Colours[HistTypes[i]]);

    //Merge NC and Dirt Events.                                                                                                                                                     
    //Add to the final historams
    //    if(HistTypes[i] == "OscNuE"){
    //  VisibleEnergyHistFinal["OscNuE"]->Add(VisibleEnergyHist[HistTypes[i]]); 
    // }

    if(HistTypes[i] == "NuMu"){
      VisibleEnergyHistFinal["NuMu F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "InNuE"){
      VisibleEnergyHistFinal["InNuE F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "NCInNuE" || HistTypes[i] == "NCNuMu"){
      VisibleEnergyHistFinal["NC F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "DirtInNuE" || HistTypes[i] == "DirtNuMu" || HistTypes[i] == "DirtNCInNuE" || HistTypes[i] == "DirtNCNuMu"){
      //      VisibleEnergyHistFinal["Dirt F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

    if(HistTypes[i] == "Cosmic"){
      //      VisibleEnergyHistFinal["Cosmic F"]->Add(VisibleEnergyHist[HistTypes[i]]);
    }

  //    VisibleEnergy_Selection_Stack->Add(VisibleEnergyHist[HistTypes[i]]);
  }

  //Add to the stack 
  for(int i=0;i< HistTypesFinal.size(); ++i){
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetFillColor(Colours[HistTypesFinal[i]]);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetMarkerStyle(21);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetMarkerColor(Colours[HistTypesFinal[i]]);
    VisibleEnergyHistFinal[HistTypesFinal[i]]->SetLineColor(Colours[HistTypesFinal[i]]);
    VisibleEnergy_Selection_Stack->Add(VisibleEnergyHistFinal[HistTypesFinal[i]]);
  }

  //TFile to hold the  selection plots. 
  TFile *Figures = new TFile("Figures.root","RECREATE");
  VisibleEnergy_Selection_Stack->Write();

  TCanvas *StackCanvas = new TCanvas("StackCanvas","StackCanvas",10,10,700,900);
  VisibleEnergy_Selection_Stack->Draw("HIST");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  StackCanvas->Write();

}
