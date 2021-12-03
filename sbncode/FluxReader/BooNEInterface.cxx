#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <regex>

#include "TFile.h"
#include "TTree.h"

#include "Framework/ParticleData/PDGUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "BooNEInterface.h"


namespace fluxr {
  BooNEInterface::BooNEInterface()
  {
  }

  BooNEInterface::~BooNEInterface()
  {
  }

  void BooNEInterface::SetRootFile(TFile* fluxInputFile)
  {
    fBNBtree = dynamic_cast<TTree*>(fluxInputFile->Get("h201"));
    fBNBtree->SetBranchAddress("beamwgt",&fBooneNtp.beamwgt);
    fBNBtree->SetBranchAddress("ntp",&fBooneNtp.ntp);
    fBNBtree->SetBranchAddress("npart",&fBooneNtp.npart);
    fBNBtree->SetBranchAddress("id",fBooneNtp.id);
    fBNBtree->SetBranchAddress("ini_pos",&fBooneNtp.ini_pos[0][0]);
    fBNBtree->SetBranchAddress("ini_mom",&fBooneNtp.ini_mom[0][0]);
    fBNBtree->SetBranchAddress("ini_eng",fBooneNtp.ini_eng);
    fBNBtree->SetBranchAddress("ini_t",fBooneNtp.ini_t);
    fBNBtree->SetBranchAddress("fin_mom",&fBooneNtp.fin_mom[0][0]);
    fBNBtree->SetBranchAddress("fin_pol",&fBooneNtp.fin_pol[0][0]);

    fWindowtree = dynamic_cast<TTree*>(fluxInputFile->Get("h220"));
    fWindowtree->SetBranchAddress("tank_pos_beam",&fBeamNtp.tank_pos_beam[0]);
    fWindowtree->SetBranchAddress("targ_pos_beam",&fBeamNtp.targ_pos_beam[0]);
    fWindowtree->SetBranchAddress("pot",&fBeamNtp.pot);
    fWindowtree->SetBranchAddress("windowbase",&fBeamNtp.windowbase[0]);
    fWindowtree->SetBranchAddress("windowdir1",&fBeamNtp.windowdir1[0]);
    fWindowtree->SetBranchAddress("windowdir2",&fBeamNtp.windowdir2[0]);


    fNEntries = fBNBtree->GetEntries();
    mf::LogInfo("BooNEInterface") << "Reading " << fNEntries << " events" << std::endl;

    fMaxWeight = 0;
    fMinWeight = 100;
    fSumWeights = 0;

    std::string run_number = fluxInputFile->GetName();
    std::regex r("\\d\\d\\d\\d");
    std::smatch m;
    std::regex_search(run_number,m,r);

    std::istringstream buffer(m[0].str().c_str());
    buffer >> fRun;


    for ( int iev = 0 ; iev < fNEntries ; iev++ ) {
      fBNBtree->GetEntry(iev);
      fSumWeights += fBooneNtp.beamwgt;
      if( fMaxWeight < fBooneNtp.beamwgt) fMaxWeight = fBooneNtp.beamwgt;
      if( fMinWeight > fBooneNtp.beamwgt) fMinWeight = fBooneNtp.beamwgt;
    }

    mf::LogInfo("BooNEInterface") << "Min Weight: " << fMinWeight
                                  << ", Max Weight: " << fMaxWeight
                                  << ", Sum Weights: " << fSumWeights
                                  << std::endl;


    fWindowtree->GetEntry(0);
    fPOT = fBeamNtp.pot;

    std::cout << "Window:" << std::endl;
    std::cout << "\tBase: " << fBeamNtp.windowbase[0] << ", " << fBeamNtp.windowbase[1] << ", " << fBeamNtp.windowbase[2] << std::endl;
    std::cout << "\tDir 1: " << fBeamNtp.windowdir1[0] << ", " << fBeamNtp.windowdir1[1] << ", " << fBeamNtp.windowdir1[2] << std::endl;
    std::cout << "\tDir 2: " << fBeamNtp.windowdir2[0] << ", " << fBeamNtp.windowdir2[1] << ", " << fBeamNtp.windowdir2[2] << std::endl;
  }


  bool BooNEInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux)
  {
    if (!fBNBtree->GetEntry(ientry))
      return false;

    mf::LogDebug("BooNEInterface") << "At entry: " << ientry << std::endl;

    fWindowtree->GetEntry(0);


    if ( fBooneNtp.ntp == 1 ) {
      flux.fntype = 12; //nue
    }
    else if ( fBooneNtp.ntp == 2 ) {
      flux.fntype = -12; //nuebar
    }
    else if ( fBooneNtp.ntp == 3 ) {
      flux.fntype = 14; //numu
    }
    else if ( fBooneNtp.ntp == 4 ) {
      flux.fntype = -14; //numubar
    }
    else{
      mf::LogWarning("BooNEInterface") << "Neutrino type not recognized! ntp = " << fBooneNtp.ntp
                                       << std::endl;
    }

    flux.fnimpwt = fBooneNtp.beamwgt;

    // Convert to meters
    double nu_x = fBooneNtp.ini_pos[0][0] / 100;
    double nu_y = fBooneNtp.ini_pos[0][1] / 100;
    double nu_z = fBooneNtp.ini_pos[0][2] / 100;

    double nu_momx = fBooneNtp.ini_mom[0][0];
    double nu_momy = fBooneNtp.ini_mom[0][1];
    double nu_momz = fBooneNtp.ini_mom[0][2];

    double targ_x = fBeamNtp.targ_pos_beam[0] / 100;
    double targ_y = fBeamNtp.targ_pos_beam[1] / 100;
    double targ_z = fBeamNtp.targ_pos_beam[2] / 100;

    double tank_z = (fBeamNtp.tank_pos_beam[2] + fBeamNtp.windowbase[2]) / 100;

    // Calculate the neutrino x, y, and z position at the window location
    double gsimple_vtxx = nu_x + (nu_momx / nu_momz) * (tank_z - nu_z) + targ_x;
    double gsimple_vtxy = nu_y + (nu_momy / nu_momz) * (tank_z - nu_z) + targ_y;
    double gsimple_vtxz = tank_z+targ_z;

    // Calculate the distance from the neutrino production point to the window
    double dist = std::sqrt(std::pow(nu_x-gsimple_vtxx,2) +
                            std::pow(nu_y-gsimple_vtxy,2) +
                            std::pow(nu_z-tank_z,2));
    flux.fdk2gen = dist;

    // Calculate the neutrino time at the window
    float nu_time = fBooneNtp.ini_t[0]; // ns
    nu_time += dist / (TMath::C() / 1e9);

    mf::LogDebug("BooNEInterface") << "weight = " << flux.fnimpwt
                                   << ", nu time = " << nu_time
                                   << ", nu energy = " << fBooneNtp.ini_eng[0]
                                   << ", dist = " << dist
                                   << std::endl;

    fNuPos = TLorentzVector(gsimple_vtxx * 100,
                            gsimple_vtxy * 100,
                            gsimple_vtxz * 100,
                            nu_time);
    fNuMom = TLorentzVector(fBooneNtp.ini_mom[0][0],
                            fBooneNtp.ini_mom[0][1],
                            fBooneNtp.ini_mom[0][2],
                            fBooneNtp.ini_eng[0]);

    flux.fnenergyn = flux.fnenergyf = fBooneNtp.ini_eng[0];

    flux.frun      = fRun;
    flux.fevtno    = ientry;

    int npart      = fBooneNtp.npart;
    flux.ftpx      = fBooneNtp.ini_mom[npart-2][0]; // npart-2
    flux.ftpy      = fBooneNtp.ini_mom[npart-2][1]; // npart-2
    flux.ftpz      = fBooneNtp.ini_mom[npart-2][2]; // npart-2

    flux.fvx       = fBooneNtp.ini_pos[0][0]; //0
    flux.fvy       = fBooneNtp.ini_pos[0][1]; //0
    flux.fvz       = fBooneNtp.ini_pos[0][2]; //0

    flux.fpdpx     = fBooneNtp.fin_mom[1][0]; //1 final
    flux.fpdpy     = fBooneNtp.fin_mom[1][1]; //1 final
    flux.fpdpz     = fBooneNtp.fin_mom[1][2]; //1 final

    flux.fpppz     = fBooneNtp.ini_mom[1][2]; //1 init
    double pppx    = fBooneNtp.ini_mom[1][0]; //1 init
    double pppy    = fBooneNtp.ini_mom[1][1]; //1 init
    double apppz = flux.fpppz;
    if (TMath::Abs(flux.fpppz) < 1.0e-30) apppz = 1.0e-30;
    flux.fppdxdz   = pppx / apppz;
    flux.fppdydz   = pppy / apppz;


    flux.fppmedium = 0.; // empty


    int ptype_input = fBooneNtp.id[1];
    int tptype_input = fBooneNtp.id[npart-2];

    if (ptype_input != 0) ptype_input = genie::pdg::GeantToPdg(ptype_input);
    if (tptype_input != 0) tptype_input= genie::pdg::GeantToPdg(tptype_input);

    flux.fptype    = ptype_input;
    flux.ftptype   = tptype_input;


    /////
    //  Now need to calculate ndecay
    /////

    double Nenergy = fBooneNtp.ini_eng[0];
    double Ndxdz   = fBooneNtp.ini_mom[0][0] / fBooneNtp.ini_mom[0][2];
    double Ndydz   = fBooneNtp.ini_mom[0][1] / fBooneNtp.ini_mom[0][2];

    double ppenergy = fBooneNtp.ini_eng[1];
    double pdPx     = fBooneNtp.fin_mom[1][0];
    double pdPy     = fBooneNtp.fin_mom[1][1];
    double pdPz     = fBooneNtp.fin_mom[1][2];

    double ppdxdz   = fBooneNtp.ini_mom[1][0] / fBooneNtp.ini_mom[1][2];
    double ppdydz   = fBooneNtp.ini_mom[1][1] / fBooneNtp.ini_mom[1][2];
    double pppz     = fBooneNtp.ini_mom[1][2];

    // Get the neutrino energy in the parent decay cm
    double parent_mass = std::sqrt(ppenergy * ppenergy -
                                   pppz * pppz * (ppdxdz * ppdxdz +
                                   ppdydz * ppdydz +
                                   1.));

    double parent_energy = std::sqrt(pdPx * pdPx +
                                     pdPy * pdPy +
                                     pdPz * pdPz +
                                     parent_mass * parent_mass);

    double gamma = parent_energy / parent_mass;
    double beta[3];
    beta[0] = pdPx / parent_energy;
    beta[1] = pdPy / parent_energy;
    beta[2] = pdPz / parent_energy;

    double partial = fBooneNtp.ini_mom[0][2] * gamma *
                     (beta[0] * Ndxdz +
                     beta[1] * Ndydz +
                     beta[2]);

    double Necm = gamma * Nenergy - partial;

    if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 1) flux.fndecay = 1;
    else if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 2) flux.fndecay = 2;
    else if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 3) flux.fndecay = 3;
    else if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 4) flux.fndecay = 4;
    else if (fBooneNtp.id[1] == 11 && fBooneNtp.ntp == 3) {
      //check if it is a two or three body decay
      if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-Necm)/Necm <= 0.001)
        //two body decay (numu + mu+)
        flux.fndecay = 5;
      else {
        //three body decay (numu + pi0 + mu+)
        flux.fndecay = 7;
      }
    }
    else if (fBooneNtp.id[1] == 11 && fBooneNtp.ntp == 1) flux.fndecay = 6;
    else if (fBooneNtp.id[1] == 12 && fBooneNtp.ntp == 4) {
      if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-Necm)/Necm <= 0.001) {
        //two body decay (numu + mu+)
        flux.fndecay = 8;
      }
      else {
        //three body decay (numu + pi0 + mu+)
        flux.fndecay = 10;
      }
    }
    else if (fBooneNtp.id[1] == 12 && fBooneNtp.ntp == 2) flux.fndecay = 9;
    else if (fBooneNtp.id[1] == 5 ) flux.fndecay = 11;
    else if (fBooneNtp.id[1] == 6 ) flux.fndecay = 12;
    else if (fBooneNtp.id[1] == 8 ) flux.fndecay = 13;
    else if (fBooneNtp.id[1] == 9 ) flux.fndecay = 14;

    /////
    //  End calculation of ndecay
    /////

    double mupare;
    double muparpx;
    double muparpy;
    double muparpz;

    if ( fBooneNtp.id[1] == 5 || fBooneNtp.id[1] == 6) {
      mupare  = fBooneNtp.ini_eng[2];
      muparpx = fBooneNtp.fin_mom[2][0];
      muparpy = fBooneNtp.fin_mom[2][1];
      muparpz = fBooneNtp.fin_mom[2][2];
    } else {
      mupare  = -9999.;
      muparpx = -9999.;
      muparpy = -9999.;
      muparpz = -9999.;
    }

    flux.fmuparpx = muparpx;
    flux.fmuparpy = muparpy;
    flux.fmuparpz = muparpz;
    flux.fmupare = mupare;

    flux.fnwtnear = flux.fnwtfar = 1;

    return true;
  }
}
