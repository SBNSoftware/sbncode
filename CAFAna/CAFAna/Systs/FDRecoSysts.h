// FDRecoSysts.h
// Uses genie and neut file kinematic discrepancies to simulate reconstruction inefficiencies
// Use make_FD_reco_systs.C to make the histograms
#pragma once

#include "CAFAna/Core/ISyst.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <cassert>

namespace ana {

 //  class FDRecoNumuSyst: public ISyst
 //  {
 //  public:
 //  FDRecoNumuSyst() : ISyst("FDRecoNumuSyst", "Far Detector Numu Reconstruction Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr,
	//        double& weight) const override 
 //    {
 //      // Load histograms if they have not been loaded already
 //      if (!hist) {
	// TFile f((FindCAFAnaDir()+"/Systs/modelComp.root").c_str());
	// assert(!f.IsZombie());
	// hist = (TH2*)f.Get("hYratio_neutfhc_geniefhc");
	// hist->SetDirectory(0);
	// assert(hist);
 //      }
 //      // Passes FD selection cut
 //      if (sr->dune.isFD && sr->dune.cvnnumu >= 0.5) {
	// int EBin   = hist->GetXaxis()->FindBin(sr->dune.Ev);
	// int VarBin = hist->GetYaxis()->FindBin(sr->dune.Y);
	// double w   = hist->GetBinContent(EBin, VarBin);
	// weight    *= 1. + sigma * (1. - w) ;
 //      }
 //    }
    
 //  protected:
 //    mutable TH2* hist;
 //  }; 

 //  extern const FDRecoNumuSyst kFDRecoNumuSyst;

 //  // Nue reco syst
 //  class FDRecoNueSyst: public ISyst
 //  {
 //  public:
 //  FDRecoNueSyst() : ISyst("FDRecoNueSyst", "Far Detector Nue Reconstruction Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr,
	//        double& weight) const override 
 //    {
 //      // Load histograms if they have not been loaded already
 //      if (!hist) {
	// TFile f((FindCAFAnaDir()+"/Systs/modelComp.root").c_str());
	// assert(!f.IsZombie());
	// hist = (TH2*)f.Get("hYratio_neutfhc_geniefhc");
	// hist->SetDirectory(0);
	// assert(hist);
 //      }
 //      // Passes FD nue selection
 //      if (sr->dune.isFD && sr->dune.cvnnue >= 0.5) {
	// int EBin   = hist->GetXaxis()->FindBin(sr->dune.Ev);
	// int VarBin = hist->GetYaxis()->FindBin(sr->dune.Y);
	// double w   = hist->GetBinContent(EBin, VarBin);
	// weight    *= 1. + sigma * (1. - w) ;
 //      }
 //    }
    
 //  protected:
 //    mutable TH2* hist;
 //  };

 //  extern const FDRecoNueSyst kFDRecoNueSyst;

 //  struct FDRecoSystVector: public std::vector<const ISyst*>
 //  {
 //  };

 //  FDRecoSystVector GetFDRecoSysts();
  
}
