// NDRecoSysts.h
// Systematics to simulate reconstruction systematics in the ND
#pragma once

#include "CAFAna/Core/ISyst.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <cassert>

namespace ana {
  
 //  // Take ND events which pass the CC selection cuts but are NC and reweight by 20%
 //  class RecoNCSyst: public ISyst
 //  {
 //  public:
 //  RecoNCSyst() : ISyst("RecoNCSyst", "ND Neutral Current Reconstruction Syst") {}
 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      // Is ND
 //      if(!sr->dune.isFD) {
	// // CC event selection but is NC
	// if((sr->dune.reco_numu || sr->dune.reco_nue) && (sr->dune.muon_contained || sr->dune.muon_tracker) && (sr->dune.reco_q == -1 || sr->dune.reco_q == 1) && sr->dune.Ehad_veto<30 && !sr->dune.isCC) {
	//   weight *= 1 + .2*sigma;
	// }
 //      }
 //    }
 //  };
 //  extern const RecoNCSyst kRecoNCSyst;

 //  // Systematic designed to mimic effect of lepton acceptance in the ND
 //  class LeptonAccSyst: public ISyst
 //  {
 //  public:
 //  LeptonAccSyst() : ISyst("LeptonAccSyst", "ND Lepton Acceptance Syst") 
 //      {
	// hist = 0;
 //      }
          
 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      // Load hist if it hasn't been loaded already
 //      const double m_mu = 0.105658;
 //      if (!hist) {
	// TFile f("/dune/app/users/marshalc/ND_syst/ND_eff_syst.root", "read");
	// assert(!f.IsZombie());
	// hist = (TH2*)f.Get("unc");
 //        hist->SetDirectory(0);
 //      }
   
 //      // Is ND and is a true numu CC event
 //      if (!sr->dune.isFD && sr->dune.isCC && abs(sr->dune.nuPDG) == 14) {
	// double LepE = sr->dune.LepE;
	// int bin = hist->FindBin(sqrt(LepE*LepE - m_mu*m_mu) * cos(sr->dune.LepNuAngle), sqrt(LepE*LepE - m_mu*m_mu) * sin(sr->dune.LepNuAngle));
	// double w = hist->GetBinContent(bin);
	// weight *= 1. + w*sigma;
     
 //      }
 //    }
    
 //  protected:
 //    mutable TH2 *hist;
    
 //  };
 //  extern const LeptonAccSyst kLeptonAccSyst;

 //  // Systematic designed to mimic effect of hadron acceptance in the ND
 //  class HadronAccSyst: public ISyst
 //  {
 //  public:
 //  HadronAccSyst() : ISyst("HadronAccSyst", "ND Hadron Acceptance Syst") 
 //      {
	// hist = 0;
 //      }
      
 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      // Load hist if it hasn't been loaded already
 //      if (!hist) {
	// TFile f("/dune/app/users/marshalc/ND_syst/ND_eff_syst.root", "read");
	// assert(!f.IsZombie());
	// hist = (TH1*)f.Get("hunc");
 //        hist->SetDirectory(0);
 //      }
      
 //      // Is ND
 //      if (!sr->dune.isFD) {
	// double HadE = sr->dune.Ev - sr->dune.LepE;
	// if (HadE > 5.) {
	//   HadE = 5.;
	// }
	// int bin = hist->FindBin(HadE);
	// double w = hist->GetBinContent(bin);
	// weight *= 1. + w*sigma;  
 //      }
 //    }
   
 //  protected:
 //    mutable TH1 *hist;
 //    //    TFile *f;
 //  };
 //  extern const HadronAccSyst kHadronAccSyst;

 //  struct NDRecoSystVector: public std::vector<const ISyst*>
 //  {

 //  };
 //  NDRecoSystVector GetNDRecoSysts();

}
