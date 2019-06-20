#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>
#include "CrossSections.h"

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

#include <TFile.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TMath.h>
#include <core/Event.hh>

#include <fstream>
#include <iostream>

namespace ana {
namespace SBNOsc {

void CrossSections::Initialize(fhicl::ParameterSet* config) {
    
    // must have config
    assert(config != NULL);
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("CrossSections");
    fOutputFile = pconfig.get<std::string>("OutputFile", "");

    fMinX = pconfig.get<double>("MinX");
    fMaxX = pconfig.get<double>("MaxX");
    fMinY = pconfig.get<double>("MinY");
    fMaxY = pconfig.get<double>("MaxY");
    fMinZ = pconfig.get<double>("MinZ");
    fMaxZ = pconfig.get<double>("MaxZ");

    fMinContainedLength = pconfig.get<double>("MinContainedLength");
    fMinExitingLength   = pconfig.get<double>("MinExitingLength");

    fProtonThreshold = pconfig.get<double>("ProtonThreshold");
    fPionThreshold   = pconfig.get<double>("PionThreshold");
    fMuonThreshold   = pconfig.get<double>("MuonThreshold");
    fPi0Threshold    = pconfig.get<double>("Pi0Threshold");
    
    fProtonPidEff = pconfig.get<double>("ProtonPidEff");
    fPionPidEff   = pconfig.get<double>("PionPidEff");

    // Declare histograms
    std::vector<std::string> recoChannels {"Inc", "0pi", "1pi", "geq2pi", "pi0", "Oth"};
    std::vector<std::string> trueTypes {"QE", "RES", "DIS", "COH", "MEC", "NC", "NuE", "Oth"};
    std::vector<std::string> trueFsis {"0pi", "1pi", "geq2pi", "pi0", "Oth"};

    std::vector<std::string> keys;
    for(auto const& channel : recoChannels){
      for(auto const& type : trueTypes){
        keys.push_back("Reco" + channel + "TrueType" + type);
      }
      for(auto const& fsi : trueFsis){
        keys.push_back("Reco" + channel + "TrueFsi" + fsi);
      }
    }
    for(auto const& key : keys){
      hRecoNuE[key]     = new TH1D(Form("hRecoNuE%s", key.c_str()),     "", 20, 0,    3);
      hRecoMuP[key]     = new TH1D(Form("hRecoMuP%s", key.c_str()),     "", 20, 0,    3);
      hRecoMuTheta[key] = new TH1D(Form("hRecoMuTheta%s", key.c_str()), "", 20, -3.2, 3.2);
      hTrueNuE[key]     = new TH1D(Form("hTrueNuE%s", key.c_str()),     "", 20, 0,    3);
      hTrueMuP[key]     = new TH1D(Form("hTrueMuP%s", key.c_str()),     "", 20, 0,    3);
      hTrueMuTheta[key] = new TH1D(Form("hTrueMuTheta%s", key.c_str()), "", 20, -3.2, 3.2);
    }
    
    hNuETrue   = new TH1D("hNuETrue", "", 100, 0, 5);
    hNuEReco   = new TH1D("hNuEReco", "", 100, 0, 5);

    fRandom = new TRandom2();

}

void CrossSections::ProcessEvent(const Event *event) {
    // iterate over each interaction in the event
    for (int n = 0; n < event->reco.size(); n++) {
        unsigned truth_ind = event->reco[n].truth_index;

        // Get energy
        double true_nuE = event->reco[n].truth.neutrino.energy;
        hNuETrue->Fill(true_nuE);
        // Add up visible reco energy
        double reco_nuE = event->reco[n].reco_energy;
        hNuEReco->Fill(reco_nuE);

        // Get lepton momentum and angle
        double true_lepP = (event->reco[n].truth.lepton.momentum).Mag();
        // FIXME start and end are the same
        double true_lepTheta = (event->reco[n].truth.lepton.end - event->reco[n].truth.lepton.start).Theta();

        // Determine true interaction type (QE, RES, MEC, DIS, COH) (CC, NC) (NuMu, NuE)
        // 0 (NC), 1 (CC)
        int true_cc = event->reco[n].truth.neutrino.iscc;
        // 0 (QE), 1 (RES), 2 (DIS), 3/4 (COH), 10 (MEC)...
        int true_int = event->reco[n].truth.neutrino.genie_intcode;
        // (-)12 NuE, (-)14 NuMu
        int true_nuPdg = event->reco[n].truth.neutrino.pdg;

        // Label the true interaction type
        std::string type = "Oth";
        if(std::abs(true_nuPdg) == 14){
          if(true_cc){
            if(true_int == 0) type = "QE";
            if(true_int == 1) type = "RES";
            if(true_int == 2) type = "DIS";
            if(true_int == 3 || true_int == 4) type = "COH";
            if(true_int == 10) type = "MEC";
          }
          else{
            type = "NC";
          }
        }
        else if(std::abs(true_nuPdg) == 12) type = "NuE";

        // Determine true FSI type (based on stable final states) (CC, NC) (0pi, 1pi, 2+pi, pi0) (Np)
        int true_nP = 0, true_nPi = 0, true_nPi0 = 0;
        int nFs = event->reco[n].truth.finalstate.size();
        for(int i = 0; i < nFs; i++){
          if(event->reco[n].truth.finalstate[i].status_code != 1) continue;
          int fs_pdg = event->reco[n].truth.finalstate[i].pdg;
          if(fs_pdg == 2212) true_nP++;
          if(std::abs(fs_pdg) == 211) true_nPi++;
          if(fs_pdg == 111) true_nPi0++;
        }

        // Label the true FSI
        std::string fsi = "Oth";
        if(true_nPi0 > 0) fsi = "pi0";
        else if(true_nPi == 0) fsi = "0pi";
        else if(true_nPi == 1) fsi = "1pi";
        else if(true_nPi >= 2) fsi = "geq2pi";

        // Determine the reco FSI type
        // Run proposal CC inclusive selection
        bool cc_selected = true;
        // Get neutrino vertex
        TVector3 true_vtx = event->reco[n].truth.neutrino.position;
        // Apply fiducial volume cut
        if(true_vtx.X() < fMinX || true_vtx.X() > fMaxX || 
           true_vtx.Y() < fMinY || true_vtx.Y() > fMaxY || 
           true_vtx.Z() < fMinZ || true_vtx.Z() > fMaxZ) cc_selected = false;

        // Get length of longest track in event
        bool longest_contained = false;
        double longest_length = 0;
        int longest_i = -1;
        double reco_lepP = 0;
        double reco_lepTheta = 0;
        for(int i = 0; i < nFs; i++){
          if(event->reco[n].truth.finalstate[i].status_code != 1) continue;
          double fs_length = event->reco[n].truth.finalstate[i].contained_length;
          if(fs_length < longest_length) continue;
          longest_length = fs_length;
          longest_i = 0;

          // Get the reconstructed "muon" quantities
          reco_lepTheta = (event->reco[n].truth.finalstate[i].end - event->reco[n].truth.finalstate[i].start).Theta();

          // If length is same as contained length then particle is contained
          if(fs_length == event->reco[n].truth.finalstate[i].length){ 
            longest_contained = true;
            reco_lepP = SmearRangeMomentum(event->reco[n].truth.finalstate[i].momentum.Mag());
          }
          else{ 
            longest_contained = false;
            reco_lepP = SmearMcsMomentum(event->reco[n].truth.finalstate[i].momentum.Mag());
          }
        }

        // Is it contained or exiting
        // Apply corresponding length cut
        if(longest_contained && longest_length < fMinContainedLength) cc_selected = false;
        if(!longest_contained && longest_length < fMinExitingLength) cc_selected = false;

        // Count up inclusive rates
        int reco_nP = 0, reco_nPi = 0, reco_nPi0 = 0;
        bool tracks_contained = true;
        // Assume muon and smear variables for plotting
        // Loop over other final state particles
        for(int i = 0; i < nFs; i++){
          if(event->reco[n].truth.finalstate[i].status_code != 1) continue;
          // Don't consider the longest track
          if(i == longest_i) continue;
          int fs_pdg = event->reco[n].truth.finalstate[i].pdg;
          double fs_momentum = (event->reco[n].truth.finalstate[i].momentum).Mag();

          // Check track is contained, if not don't count (unless pi0 contained)
          bool contained = (event->reco[n].truth.finalstate[i].contained_length == event->reco[n].truth.finalstate[i].length);

          // Apply efficiency cut (energy threshold) P = 200 MeV, Pi/Mu = 50 MeV, Pi0 = 80 MeV
          if(fs_pdg == 2212 && fs_momentum < fProtonThreshold) continue;
          if(std::abs(fs_pdg) == 211 && fs_momentum < fPionThreshold) continue;
          if(std::abs(fs_pdg) == 13 && fs_momentum < fMuonThreshold) continue;
          if(fs_pdg == 111 && fs_momentum < fPi0Threshold) continue;
          else if(fs_pdg == 111){
            // TODO Apply flat efficiency after
            if(contained) reco_nPi0++;
          }

          // For tracks apply PID estimation
          if(fs_pdg == 2212){
            if(!contained) tracks_contained = false;
            // Generate random number
            double rand = fRandom->Rndm();
            // Apply PID effiecieny of 85% for protons
            if(rand < fProtonPidEff) reco_nP++;
            else reco_nPi++;
          }

          if(std::abs(fs_pdg) == 211 || std::abs(fs_pdg) ==13){
            if(!contained) tracks_contained = false;
            // Generate random number
            double rand = fRandom->Rndm();
            // Apply PID efficiency of 99% for mu/pi
            if(rand < fPionPidEff) reco_nPi++;
            else reco_nP++;
          }
        }

        // Label the reco channels
        std::string channel = "Oth";
        if(reco_nPi0 > 0) channel = "pi0";
        else if(reco_nPi == 0) channel = "0pi";
        else if(reco_nPi == 1) channel = "1pi";
        else if(reco_nPi >= 2) channel = "geq2pi";

        // Fill histograms
        if(!cc_selected) continue;
        std::vector<std::string> keys;
        keys.push_back( "RecoIncTrueType" + type);
        keys.push_back("RecoIncTrueFsi" + fsi);
        std::string typeKey = "RecoIncTrueType" + type;
        std::string fsiKey = "RecoIncTrueFsi" + fsi;

        typeKey = "Reco" + channel + "TrueType" + type;
        fsiKey = "Reco" + channel + "TrueFsi" + fsi;
        if(channel == "pi0"){
          keys.push_back("Reco" + channel + "TrueType" + type);
          keys.push_back("Reco" + channel + "TrueFsi" + fsi);
        }
        // If any other channel require all tracks contained for PID
        else if(tracks_contained){
          keys.push_back("Reco" + channel + "TrueType" + type);
          keys.push_back("Reco" + channel + "TrueFsi" + fsi);
        }

        for(auto const& key : keys){
          hTrueNuE[key]->Fill(true_nuE);
          hTrueMuP[key]->Fill(true_lepP);
          hTrueMuTheta[key]->Fill(true_lepTheta);
          hRecoNuE[key]->Fill(reco_nuE);
          hRecoMuP[key]->Fill(reco_lepP);
          hRecoMuTheta[key]->Fill(reco_lepTheta);
        }

    }
}


void CrossSections::FileCleanup(TTree *eventTree) {
}

void CrossSections::Write(){

  std::cout<<"Write called!\n";

  if(fOutputFile == ""){ std::cout<<"Error!\n"; return;}
  TFile *output = TFile::Open(fOutputFile.c_str(), "recreate");
  assert(output && output->IsOpen());

  for(auto const& nuE : hTrueNuE){
    nuE.second->Write();
  }
  for(auto const& muP : hTrueMuP){
    muP.second->Write();
  }
  for(auto const& muTheta : hTrueMuTheta){
    muTheta.second->Write();
  }
  for(auto const& nuE : hRecoNuE){
    nuE.second->Write();
  }
  for(auto const& muP : hRecoMuP){
    muP.second->Write();
  }
  for(auto const& muTheta : hRecoMuTheta){
    muTheta.second->Write();
  }

  hNuETrue->Write();
  hNuEReco->Write();

  output->Close();

}

double CrossSections::SmearMcsMomentum(double momentum){
  //For exiting muons use multiple coulomb scattering bias and resolution
  //Values from Fig 12 of https://arxiv.org/pdf/1703.06187.pdf
  double bias[] = {0.0273,0.0409,0.0352,0.0250,0.0227,0.0068,0.0364,0.0273,0.0227};
  double resolution[] = {0.127,0.145,0.143,0.141,0.164,0.177,0.250,0.266,0.341};
  int pos = 8; 
  for(int i=0; i<9; i++){
    if(momentum<(0.34+0.41*(i+1))) {pos = i; break;}
  }
  double momentum_smear = fRandom->Gaus(momentum, resolution[pos]*momentum) + bias[pos]*momentum;
  if(momentum_smear<0) {momentum_smear = 0;}
  return momentum_smear;

}

double CrossSections::SmearRangeMomentum(double momentum){
  //For contained muons use range based bias and resolution
  //Values from Fig 5 of https://arxiv.org/pdf/1703.06187.pdf
  double bias[] = {-0.0035,-0.0059,-0.0047,-0.0059,-0.0035,-0.0029,-0.0076,-0.0059,0.0006};
  double resolution[] = {0.017,0.021,0.023,0.026,0.025,0.030,0.030,0.040,0.032};
  int pos = 8; 
  for(int i=0; i<9; i++){
    if(momentum<(0.33+0.186*(i+1))) {pos = i; break;}
  }
  double momentum_smear = fRandom->Gaus(momentum, resolution[pos]*momentum)+bias[pos]*momentum;
  if(momentum_smear<0){momentum_smear = 0;}
  return momentum_smear;                                                                                                                                            
}

    
}   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::CrossSections);

