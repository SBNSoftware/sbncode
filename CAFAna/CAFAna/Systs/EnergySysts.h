#pragma once

#include "CAFAna/Core/ISyst.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>

namespace ana
{ 
 //  // Slope energy scale systematics
 //  // Affect ND only
 //  class UncorrNDHadLinSyst: public ISyst
 //  {
 //  public:
 //  UncorrNDHadLinSyst() : ISyst("UncorrNDHadLinSyst", "Uncorrelated ND Linear Hadron Syst") {}
 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco);
 //      if (!sr->dune.isFD) {
	// const double scale = .01 * sigma;
	// double sumE = sr->dune.eP + sr->dune.ePip + sr->dune.ePim;
	// const double fracE = sumE / sr->dune.Ev;
	// sr->dune.Ev_reco += sr->dune.Ev_reco * sumE * scale * fracE;
 //      }
 //    }
 //  };

 //  extern const UncorrNDHadLinSyst kUncorrNDHadLinSyst;


 //  class UncorrNDPi0LinSyst: public ISyst
 //  {
 //  public:
 //  UncorrNDPi0LinSyst() : ISyst("UncorrNDPi0LinSyst", "Uncorrelated ND Linear Pi0 Syst") {}
 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco);
 //      if (!sr->dune.isFD) {
	// const double scale = .01 * sigma;
	// double sumE = sr->dune.ePi0;
	// const double fracE = sumE / sr->dune.Ev;
	// sr->dune.Ev_reco += sr->dune.Ev_reco * sumE * scale * fracE;
 //      }
 //    }
 //  };

 //  extern const UncorrNDPi0LinSyst kUncorrNDPi0LinSyst;

 //  // Neutron
 //  class UncorrNDNLinSyst: public ISyst
 //  {
 //  public:
 //  UncorrNDNLinSyst() : ISyst("UncorrNDNLinSyst", "Uncorrelated ND Linear N Syst") {}
 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco);
 //      if (!sr->dune.isFD) {
	// const double scale = .01 * sigma;
	// double visE = 0.25 * sr->dune.eN; // crude approximation
	// const double fracE = visE / sr->dune.Ev;
	// sr->dune.Ev_reco += sr->dune.Ev_reco * visE * scale * fracE;
 //      }
 //    }
 //  };

 //  extern const UncorrNDNLinSyst kUncorrNDNLinSyst;


 //  /// 1% systematic on muon energy for energy deposition in liquid argon
 //  /// 100% correlated between near and far detectors for those ND events that stop in the LAr
 //  class eScaleMuLArSyst: public ISyst
 //  {
 //  public:
 //  eScaleMuLArSyst() : ISyst("eScaleMuLAr", "Muon Energy Scale LAr") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco,
 //                  sr->dune.Elep_reco,
 //                  sr->dune.Ev_reco_numu,
 //                  sr->dune.RecoLepEnNumu);

 //      const double scale = 1 + .01*sigma;

 //      // Checks if ND
 //      if(!sr->dune.isFD){
	// // Select only CC muon neutrino events that stop in LAr
 //        if(abs(sr->dune.nuPDG) == 14 && sr->dune.isCC == 1 && sr->dune.muon_contained == 1){
 //          sr->dune.Ev_reco   = sr->dune.Ev_reco * ( (1 - sr->dune.Y) * scale + sr->dune.Y );
	//   sr->dune.Elep_reco = sr->dune.Elep_reco * scale;
 //        }
	// else { }
 //      }
 //      // Otherwise is FD
 //      else {
	// if(sr->dune.isCC && abs(sr->dune.nuPDG) == 14){
	//   sr->dune.Ev_reco_numu = sr->dune.Ev_reco_numu *( (1 - sr->dune.Y) * scale + sr->dune.Y);
	//   sr->dune.RecoLepEnNumu = sr->dune.RecoLepEnNumu * scale;
	// }
 //      }
 //    }
 //  };

 //  extern const eScaleMuLArSyst keScaleMuLArSyst;

 //  /// 1% systematics on muon energy for those tracks that are measured by the magnetic field
 //  // Uncorrelated between ND and FD
 //  // This is a temporary solution - need some momentum dependent function
 //  class EnergyScaleMuSystND: public ISyst
 //  {
 //  public:
 //  EnergyScaleMuSystND() : ISyst("eScaleMuND", "Muon Energy Scale Near Detector") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco,
 //                  sr->dune.Elep_reco);

 //      const double scale = 1 + .01*sigma;

 //      // Is a numu CC and enters the tracker
 //      if(!sr->dune.isFD && abs(sr->dune.nuPDG)==14 && sr->dune.isCC==1 && sr->dune.muon_tracker==1){
	// sr->dune.Ev_reco = sr->dune.Ev_reco * ( (1 - sr->dune.Y) * scale + sr->dune.Y );
	// sr->dune.Elep_reco = sr->dune.Elep_reco * scale;
 //      }
 //    }
 //  };

 //  extern const EnergyScaleMuSystND kEnergyScaleMuSystND;

 //  /// 2% energy scale systematic on electron energy
 //  /// 100% correlated between near and far detectors
 //  class EnergyScaleESyst: public ISyst
 //  {
 //  public:
 //  EnergyScaleESyst() : ISyst("eScaleE", "Electron Energy Scale") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco,
 //                  sr->dune.Ev_reco_nue,
 //                  sr->dune.Elep_reco,
 //                  sr->dune.RecoLepEnNue);

 //      const double scale = 1 + .02*sigma;

 //      // Checks if ND
 //      if(!sr->dune.isFD){
 //        if(abs(sr->dune.nuPDG) == 12 && sr->dune.isCC){
	//   sr->dune.Ev_reco = sr->dune.Ev_reco * ( (1 - sr->dune.Y) * scale + sr->dune.Y);
	//   sr->dune.Elep_reco = sr->dune.Elep_reco * scale;
 //        }
 //      }
 //      // Otherwise is FD
 //      else {
 //        if(sr->dune.isCC && abs(sr->dune.nuPDG) == 12){
 //          sr->dune.Ev_reco_nue = sr->dune.Ev_reco_nue * ( (1 - sr->dune.Y) * scale + sr->dune.Y );
	//   sr->dune.RecoLepEnNue = sr->dune.RecoLepEnNue * scale;
 //        }
 //      }
 //    }
 //  };

 //  extern const EnergyScaleESyst kEnergyScaleESyst;


 //  /// Energy scale systematics for hadronic final state particles
 //  // 5% near/far correlated part for charged hadrons
 //  class ChargedHadCorrSyst: public ISyst
 //  {
 //  public:
 //  ChargedHadCorrSyst() : ISyst("ChargedHadCorr", "Charged Hadron Correlated Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco,
 //                  sr->dune.Ev_reco_nue,
 //                  sr->dune.Ev_reco_numu,
 //                  sr->dune.RecoHadEnNumu,
 //                  sr->dune.RecoHadEnNue);

 //      const double scale = 1. + 0.05*sigma;
 //      double sumE = 0.;
 //      sumE = sr->dune.eP + sr->dune.ePim + sr->dune.ePip;
 //      const double fracE = sumE / sr->dune.Ev;
 //      const double fracEY = sumE / (sr->dune.Ev * sr->dune.Y);
 //      sr->dune.Ev_reco = sr->dune.Ev_reco * (fracE * scale + (1 - fracE));
 //      sr->dune.Ev_reco_numu = sr->dune.Ev_reco_numu * (fracE * scale + (1 - fracE));
 //      sr->dune.Ev_reco_nue = sr->dune.Ev_reco_nue * (fracE * scale + (1 - fracE));
 //      sr->dune.RecoHadEnNumu = sr->dune.RecoHadEnNumu * (fracEY * scale + (1 - fracEY));
 //      sr->dune.RecoHadEnNue = sr->dune.RecoHadEnNue * (fracEY * scale + (1 - fracEY));
 //      // Want to apply this syst to the reco lepton energy if we have a pion misID'd as a muon
 //      /*
 //      if (!sr->dune.isFD && sr->dune.reco_numu == 1 && !sr->dune.isCC) {
	// sr->dune.Elep_reco = sr->dune.Elep_reco * ( (1 - fracE) * scale + fracE );
 //      }
 //      else if (sr->dune.isFD && !sr->dune.isCC && sr->dune.cvnnumu > 0.5) {
	// sr->dune.Ev_reco_numu = sr->dune.Ev_reco_numu * ( (1 - fracE) * scale + fracE);
	// sr->dune.RecoLepEnNumu = sr->dune.RecoLepEnNumu * ( (1 - fracEY) * scale + fracEY);
 //      } 
 //      */
 //    }
 //  };

 //  extern const ChargedHadCorrSyst kChargedHadCorrSyst;

 //  // 1% uncorrelated FD syst for charged hadrons
 //  class ChargedHadUncorrFDSyst: public ISyst
 //  {
 //  public:
 //  ChargedHadUncorrFDSyst() : ISyst("ChargedHadUncorrFD", "Charged Hadron Uncorrelated FD Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco_nue,
 //                  sr->dune.Ev_reco_numu,
 //                  sr->dune.RecoHadEnNumu,
 //                  sr->dune.RecoHadEnNue);

 //      const double scale = 1. + 0.01*sigma;
      
 //      // TEMPORARY FIX: CHANGE BACK AFTER CAFs HAVE BEEN RERUN
 //      if(sr->dune.isFD) { 
	// const double sumE = sr->dune.eP + sr->dune.ePim + sr->dune.ePip;

	// const double fracE = sumE / sr->dune.Ev;
	// const double fracEY = sumE / (sr->dune.Ev * sr->dune.Y);

	// sr->dune.Ev_reco_numu = sr->dune.Ev_reco_numu * (fracE * scale + (1 - fracE));
	// sr->dune.Ev_reco_nue = sr->dune.Ev_reco_nue * (fracE * scale + (1 - fracE));
	// sr->dune.RecoHadEnNumu = sr->dune.RecoHadEnNumu * (fracEY * scale + (1 - fracEY));
	// sr->dune.RecoHadEnNue = sr->dune.RecoHadEnNue * (fracEY * scale + (1 - fracEY));
 //      }
 //    }
 //  };
  
 //  extern const ChargedHadUncorrFDSyst kChargedHadUncorrFDSyst;

 //  /// 1% uncorrelated ND syst for charged hadrons
 //  class ChargedHadUncorrNDSyst: public ISyst
 //  {
 //  public:
 //  ChargedHadUncorrNDSyst() : ISyst("ChargedHadUncorrND", "Charged Hadron Uncorrelated ND Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco);

 //      const double scale = 1. + 0.01*sigma;
      
 //      if(!sr->dune.isFD) { 
	// const double sumE = sr->dune.eP + sr->dune.ePim + sr->dune.ePip;
	// const double fracE = sumE / sr->dune.Ev;
	// sr->dune.Ev_reco = sr->dune.Ev_reco * (fracE * scale + (1 - fracE));
 //      }
 //    }
 //  };
  
 //  extern const ChargedHadUncorrNDSyst kChargedHadUncorrNDSyst;

 //  // Assume 25% of the neutron energy is visible - this is fairly crude and should be changed later
 //  // Neutron energy scale
 //  class NUncorrNDSyst: public ISyst
 //  {
 //  public:
 //  NUncorrNDSyst() : ISyst("eScaleN_ND", "Neutron Energy Scale ND") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco);

 //      const double scale = .20 * sigma;

 //      double visE = 0.; // neutron visible energy

 //      if(!sr->dune.isFD) {
	// visE = sr->dune.eN * .25; // crude assumption
	
	// sr->dune.Ev_reco       += (visE * scale);
	// sr->dune.Ev_reco_numu  += (visE * scale);
	// sr->dune.Ev_reco_nue   += (visE * scale);
 //      }   
 //    }
 //  };

 //  extern const NUncorrNDSyst kNUncorrNDSyst;  


 //  // Assume 25% of the neutron energy is visible - this is fairly crude and should be changed later
 //  // Neutron energy scale for FD
 //  class NUncorrFDSyst: public ISyst
 //  {
 //  public:
 //  NUncorrFDSyst() : ISyst("eScaleN_FD", "Neutron Energy Scale FD") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco_numu,
 //                  sr->dune.Ev_reco_nue,
 //                  sr->dune.RecoHadEnNumu,
 //                  sr->dune.RecoHadEnNue);

 //      const double scale = .20 * sigma;

 //      double visE = 0.; // neutron visible energy

 //      if(sr->dune.isFD) {
	// // CHANGE THIS ONCE CAFs ARE RERUN
	// visE = sr->dune.eN * .25; // crude assumption
	
	// double recoNueTmp = sr->dune.RecoHadEnNue;
	// double recoNumuTmp = sr->dune.RecoHadEnNumu;

	// if (sr->dune.RecoHadEnNumu < 0) { 
	//   sr->dune.RecoHadEnNumu = 0.;
	//   sr->dune.Ev_reco_numu -=recoNumuTmp;
	// }
	// else if (sr->dune.RecoHadEnNue < 0) {
	//   sr->dune.RecoHadEnNue = 0.;
	//   sr->dune.Ev_reco_nue -=recoNueTmp;
	// }
	// else {
	//   sr->dune.Ev_reco_numu  += (visE * scale);
	//   sr->dune.Ev_reco_nue   += (visE * scale);
	// }
 //      }
 //    }
 //  };

 //  extern const NUncorrFDSyst kNUncorrFDSyst;  

 //  // Pi0 energy scale correlated between near and far
 //  // 5% on reconstructed energy
 //  class Pi0CorrSyst: public ISyst
 //  {
 //  public:
 //  Pi0CorrSyst() : ISyst("eScalePi0Corr", "Pi0 Correlated Energy Scale") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco,
 //                  sr->dune.Ev_reco_nue,
 //                  sr->dune.Ev_reco_numu,
 //                  sr->dune.RecoHadEnNumu,
 //                  sr->dune.RecoHadEnNue);

 //      const double scale = 1 + .05 * sigma;
 //      double fracPi0 = 0;
 //      double fracPi0Y = 0;
           
 //      fracPi0 = (sr->dune.ePi0 / sr->dune.Ev);
 //      fracPi0Y = (sr->dune.ePi0 / (sr->dune.Ev*sr->dune.Y));
      
 //      sr->dune.Ev_reco      = sr->dune.Ev_reco * (fracPi0 * scale + (1 - fracPi0));
 //      sr->dune.Ev_reco_numu = sr->dune.Ev_reco_numu * (fracPi0 * scale + (1 - fracPi0));
 //      sr->dune.Ev_reco_nue  = sr->dune.Ev_reco_nue * (fracPi0 * scale + (1 - fracPi0));
 //      sr->dune.RecoHadEnNumu = sr->dune.RecoHadEnNumu * (fracPi0Y * scale + (1 - fracPi0Y));
 //      sr->dune.RecoHadEnNue  = sr->dune.RecoHadEnNue * (fracPi0Y * scale + (1 - fracPi0Y));
 //    }
 //  };

 //  extern const Pi0CorrSyst kEnergyScalePi0Syst;

 //  // 2% uncorrelated FD syst for pi0
 //  class Pi0UncorrFDSyst: public ISyst
 //  {
 //  public:
 //  Pi0UncorrFDSyst() : ISyst("Pi0UncorrFD", "Pi0 Uncorrelated FD Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco_nue,
 //                  sr->dune.Ev_reco_numu,
 //                  sr->dune.RecoHadEnNumu,
 //                  sr->dune.RecoHadEnNue);

 //      const double scale = 1. + 0.02*sigma;
      
 //      if(sr->dune.isFD) { 
	// const double fracPi0 = sr->dune.ePi0 / sr->dune.Ev;
	// const double fracPi0Y = sr->dune.ePi0 / (sr->dune.Ev * sr->dune.Y);

	// sr->dune.Ev_reco_numu = sr->dune.Ev_reco_numu * (fracPi0 * scale + (1 - fracPi0));
	// sr->dune.Ev_reco_nue = sr->dune.Ev_reco_nue * (fracPi0 * scale + (1 - fracPi0));
	// sr->dune.RecoHadEnNumu = sr->dune.RecoHadEnNumu * (fracPi0Y * scale + (1 - fracPi0Y));
	// sr->dune.RecoHadEnNue = sr->dune.RecoHadEnNue * (fracPi0Y * scale + (1 - fracPi0Y));
 //      }
 //    }
 //  };
  
 //  extern const Pi0UncorrFDSyst kPi0UncorrFDSyst;


 //  /// 2% uncorrelated ND syst for pi0
 //  class Pi0UncorrNDSyst: public ISyst
 //  {
 //  public:
 //  Pi0UncorrNDSyst() : ISyst("Pi0UncorrND", "Pi0Uncorrelated ND Syst") {}

 //    void Shift(double sigma,
	//        Restorer& restore,
	//        caf::StandardRecord* sr, double& weight) const override
 //    {
 //      restore.Add(sr->dune.Ev_reco);

 //      const double scale = 1. + 0.02*sigma;
      
 //      if(!sr->dune.isFD) { 
	// const double fracPi0 = sr->dune.ePi0 / sr->dune.Ev;
	// sr->dune.Ev_reco = sr->dune.Ev_reco * (fracPi0 * scale + (1 - fracPi0));
 //      }
 //    }
 //  };
  
 //  extern const Pi0UncorrNDSyst kPi0UncorrNDSyst;

 //  // Vector of energy scale systematics
 //  struct EnergySystVector: public std::vector<const ISyst*>
 //  {
 //    /*
 //    operator std::vector<const ISyst*>()
 //    {
 //      return std::vector<const ISyst*>(begin(), end());
 //    }
 //    */
 //  };

 //  EnergySystVector GetEnergySysts();

}
