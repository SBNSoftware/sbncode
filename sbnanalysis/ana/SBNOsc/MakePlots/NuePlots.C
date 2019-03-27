//C++ Includes 
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

//Root Includes 
#include "TH1.h"
#include "THStack.h"

void NuePlots(){

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
  std::vector<std::string> HistTypes = {"OscNuE","InNuE","NC","NuMu"};
  std::vector<int>         Colours   = {419,418,792,603};

  //Histrograms for the prdocutions
  std::map<std::string,TH1D*> VisibleEnergyHist;
  for(int i=0; i<HistTypes.size(); ++i){
    const char* hist_name = HistTypes[i].c_str();
    VisibleEnergyHist[HistTypes[i]] = new TH1D(hist_name, hist_name, 11,nue_bins);
  }

  //Loop over the text files and make a root file.                                                          
  std::ifstream txtfile_stream;
  txtfile_stream.open(inputfilelist);
  if (txtfile_stream.is_open()){
    float POTWeight;
    std::string inputfile;

    std::string line;

    while (txtfile_stream.good()){
    
      std::getline(txtfile_stream, line);

      if(line.size() == 0){
	continue;
      }


      std::size_t pos = line.find(" ");    
      
      inputfile = line.substr(0,pos);
      POTWeight = stof(line.substr(pos));
      std::cout << "POTWeight: " << POTWeight <<std::endl;
            
      //Get the inputfile
      const char* InputFileName = inputfile.c_str();
      TFile *InputFile = new TFile(InputFileName);
     
      //move to the Hist folder
      InputFile->cd("Histograms");
 
      for(int i=0; i<HistTypes.size(); ++i){
	
	std::cout << "New hist type" << std::endl;

	std::string  VisibleEnergy_Selection_String  = HistTypes[i] + " VisibleEnergy_Selection";
	const char*  VisibleEnergy_Selection_Name    = VisibleEnergy_Selection_String.c_str();
	
	//Get the histogram
	TH1D* VisibleEnergyHist1 = (TH1D*)(gDirectory->Get(VisibleEnergy_Selection_Name))->Clone();
	TH1D* VisibleEnergyHist2 = (TH1D*)VisibleEnergyHist1->Rebin(11,VisibleEnergy_Selection_Name,nue_bins);

	//Account for /Gev
	for(int bin=0; bin< VisibleEnergyHist2->GetNbinsX(); ++bin){
	  std::cout << "VisibleEnergyHist2->GetBinWidth(1+bin): " << VisibleEnergyHist2->GetBinWidth(1+bin) << std::endl;
	  VisibleEnergyHist2->SetBinContent(1+bin, POTWeight*VisibleEnergyHist2->GetBinContent(1+bin) / VisibleEnergyHist2->GetBinWidth(1+bin));
	}

	if(VisibleEnergyHist2 == NULL){
	  std::cerr << "Invalid Histogram initalisation" << std::endl;
	  return;
	}

	VisibleEnergyHist2->SetFillColor(Colours[i]);
	VisibleEnergyHist2->SetMarkerStyle(21);
	VisibleEnergyHist2->SetMarkerColor(Colours[i]);
	VisibleEnergyHist2->SetLineColor(Colours[i]);

	VisibleEnergyHist[HistTypes[i]]->Add(VisibleEnergyHist2);
      }
    }
  }
  else{
    std::cerr << "Could not open the input text file" << std::endl;
    return;
  }

  for(int i=0; i<HistTypes.size(); ++i){
    VisibleEnergyHist[HistTypes[i]]->SetFillColor(Colours[i]);
    VisibleEnergyHist[HistTypes[i]]->SetMarkerStyle(21);
    VisibleEnergyHist[HistTypes[i]]->SetMarkerColor(Colours[i]);
    VisibleEnergyHist[HistTypes[i]]->SetLineColor(Colours[i]);
    VisibleEnergy_Selection_Stack->Add(VisibleEnergyHist[HistTypes[i]]);
  }

  //TFile to hold the  selection plots. 
  TFile *Figures = new TFile("Figures.root","RECREATE");
  VisibleEnergy_Selection_Stack->Write();

  TCanvas *StackCanvas = new TCanvas("StackCanvas","StackCanvas",10,10,700,900);
  VisibleEnergy_Selection_Stack->Draw("HIST");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  StackCanvas->Write();

}
