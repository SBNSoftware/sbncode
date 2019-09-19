#include <cmath>
#include <vector>

#include "OscLib/OscCalculator.h"
#include "OscLib/OscCalculatorGeneral.h"
#include "OscLib/OscCalculatorPMNS.h"
#include "OscLib/OscCalculatorPMNSOpt.h"
#include "OscLib/OscCalculatorPMNS_CPT.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <fenv.h>

int GetRows(int prod);
double GetMax(double array[],int Nentries);
double GetMin(double array[],int Nentries);

int main()
{
#ifndef DARWINBUILD
  feenableexcept(FE_INVALID); // Spot any infs or nans early
#else
  std::cerr << "WARNING: OscLib/test/testOscRndm.cc was built on OS X where feenableexcept is unavailable" << std::endl;
#endif

  TRandom3 rnd;
  rnd.SetSeed(0);

  double dmsq21, dmsq32, th12, th23, th13, delta, E, rho, L;
  int anti, from, to;

  const int kNumCalcs = 5;
  const int kNumDiffs = kNumCalcs*(kNumCalcs-1)/2;
  const int kRows = GetRows(kNumDiffs);
  const int kCols = kNumDiffs/kRows;

//  std::cout << kRows << " rows and " << kCols << " columns." << std::endl;

  TCanvas* canvs = new TCanvas("canv_all","",350*kCols,350*kRows);
  canvs->Divide(kCols, kRows);

  osc::OscCalculator osc1;
  osc::OscCalculatorGeneral osc2;
  osc::OscCalculatorPMNS osc3;
  osc::OscCalculatorPMNSOpt osc4;
  osc::OscCalculatorPMNS_CPT osc5;

  osc::IOscCalculatorAdjustable* oscs[kNumCalcs] = {&osc1, &osc2, &osc3, &osc4, &osc5};
  const TString names[kNumCalcs] = {"Approx", "General", "PMNS", "PMNSOpt", "PMNS_CPT"};

  const int kIter = 10000;

  double dv[kNumDiffs][kIter];
  TString diffnames[kNumDiffs];

  for(int it=0; it<kIter; it++){ 

    if(it%(kIter/20)==0) std::cout << 100*it/kIter << "% done" << std::endl; 

    dmsq21 = rnd.Gaus(7.59e-5,0.2e-5);
    dmsq32 = rnd.Gaus(2.39e-3,0.1e-3);
    th12 = rnd.Gaus(0.6,0.02);
    th23 = TMath::ASin(sqrt(rnd.Gaus(0.957,0.035)))/2;
    th13 = TMath::ASin(sqrt(rnd.Gaus(0.098,0.013)))/2;
    if(rnd.Integer(2)==0) th23 = TMath::Pi()/2 - th23;
    if(rnd.Integer(2)==0) dmsq32 *= -1;
    delta = rnd.Rndm()*2*TMath::Pi();
    E = rnd.Rndm()*20 + 0.1;  // At least 100 MeV because of approximations
    rho = rnd.Rndm()*20;
//    rho = 0;
//    L = rnd.Rndm()*13000;
    L = 800;                  // Fix at 800 km because of approximations
    anti = 2*((int)rnd.Integer(2)) - 1;
    from = 2*((int)rnd.Integer(3)) + 12;
    to = 2*((int)rnd.Integer(3)) + 12;

    double Ps[kNumCalcs];
    for(int n = 0; n < kNumCalcs; ++n){
      oscs[n]->SetL(L);
      oscs[n]->SetRho(rho);
      oscs[n]->SetDmsq21(dmsq21);
      oscs[n]->SetDmsq32(dmsq32);
      oscs[n]->SetTh12(th12);
      oscs[n]->SetTh13(th13);
      oscs[n]->SetTh23(th23);
      oscs[n]->SetdCP(delta);
      if (names[n] == "PMNS_CPT"){
          osc5.SetDmsq21Bar(dmsq21);
          osc5.SetDmsq32Bar(dmsq32);
          osc5.SetTh12Bar(th12);
          osc5.SetTh13Bar(th13);
          osc5.SetTh23Bar(th23);
          osc5.SetdCPBar(delta);
      }
      const double P = oscs[n]->P(anti*from, anti*to, E);
      Ps[n] = P;
    }

    int k = 0;
    for(int i = 0; i < kNumCalcs; ++i){
    for(int j = i+1; j < kNumCalcs; ++j){
      dv[k][it] = Ps[i]-Ps[j];
      diffnames[k] = names[i] + " - " + names[j];
      k++;
    }} // end for i,j

  }// iterations

  TH1D *hd[kNumDiffs];

  for(int i=0; i<kNumDiffs; i++){
    double max = GetMax(dv[i],kIter);
    double min = GetMin(dv[i],kIter);
    double dx = max - min;
    min -= dx/2;
    max += dx/2;

    hd[i] = new TH1D("","",100,min,max);

    for(int j=0;j<kIter;j++){
      hd[i]->Fill(dv[i][j]);
    }

    hd[i]->SetLineWidth(2);
    hd[i]->SetTitle(";"+diffnames[i]+";Random Tests");
    canvs->cd(i+1);
    hd[i]->DrawCopy("hist");
    canvs->GetPad(i+1)->SetLogy();

  }

  TFile* f = new TFile("testOscRndm.root", "RECREATE");
  canvs->Write();
  f->Close();

  std::cout << "See output file testOscRndm.root for difference histograms" << std::endl;

  return 0;
}

int GetRows(int prod){

  int row = 1;

  for(int irow=2; irow<sqrt(prod); irow++){
    if( prod%irow == 0 ){
      row = irow;
    }
  }

  return row;

}

double GetMax(double array[],int Nentries){

  double out = array[0];

  for(int i=1;i<Nentries;i++){
    if(array[i]>out) out = array[i];
  }

  return out;

}

double GetMin(double array[],int Nentries){

  double out = array[0];

  for(int i=1;i<Nentries;i++){
    if(array[i]<out) out = array[i];
  }

  return out;

}
