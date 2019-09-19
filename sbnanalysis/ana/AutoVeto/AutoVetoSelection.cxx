#include <iostream>
#include <vector>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1F.h>
#include <set>
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "AutoVetoSelection.h"
#include "AutoVetoTools.h"
#include "core/Event.hh"
#include "util/Interaction.hh"
#include "core/ProviderManager.hh"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"

namespace ana {
  namespace AutoVetoAnalysis {

AutoVetoSelection::AutoVetoSelection()
    : SelectionBase(), fEventCounter(0),
      fNuCount(0), fCRTTag(0), fCCNC(-1), fNuReg(0) {}

void AutoVetoSelection::Initialize(fhicl::ParameterSet* config) {
  // Make a histogram
  fNuVertexXZHist = new TH2D("nu_vtx_XZ", "#nu vertex: XZ; Z [cm]; X [cm]",
                             100, -1000, 1000, 100, -450, 450);
  fNuVertexXYHist = new TH2D("nu_vtx_XY", "#nu vertex: XY; X [cm]; Y [cm]",
                             100, -450, 450, 100, -250, 250);
  fNuVertexYZHist = new TH2D("nu_vtx_YZ", "#nu vertex: YZ; Z [cm]; Y[cm]",
                             100, -1000, 1000, 100, -250, 250);
  fNuVertexXZAVHist = new TH2D("nu_vtx_AV_XZ", "#nu vertex: XZ; Z [cm]; X [cm]",
                             100, -1000, 1000, 100, -450, 450);
  fNuVertexXYAVHist = new TH2D("nu_vtx_AV_XY", "#nu vertex: XY; X [cm]; Y [cm]",
                             100, -450, 450, 100, -250, 250);
  fNuVertexYZAVHist = new TH2D("nu_vtx_AV_YZ", "#nu vertex: YZ; Z [cm]; Y[cm]",
                             100, -1000, 1000, 100, -250, 250);
  fNuVertexXZFVHist = new TH2D("nu_vtx_FV_XZ", "#nu vertex: XZ; Z [cm]; X [cm]",
                             100, -1000, 1000, 100, -450, 450);
  fNuVertexXYFVHist = new TH2D("nu_vtx_FV_XY", "#nu vertex: XY; X [cm]; Y [cm]",
                             100, -450, 450, 100, -250, 250);
  fNuVertexYZFVHist = new TH2D("nu_vtx_FV_YZ", "#nu vertex: YZ; Z [cm]; Y[cm]",
                             100, -1000, 1000, 100, -250, 250);
  fNuE            = new TH1D("nuE",       "#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEAV          = new TH1D("nuEAV",     "#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEAVNotFV     = new TH1D("nuEAVNotFV","#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEFV          = new TH1D("nuEFV",     "#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEVeto        = new TH1D("nuEVeto",       "#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEVetoAV      = new TH1D("nuEVetoAV",     "#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEVetoAVNotFV = new TH1D("nuEVetoAVNotFV","#nu energy; E_{#nu} [GeV]",100,0,10);
  fNuEVetoFV      = new TH1D("nuEVetoFV",     "#nu energy; E_{#nu} [GeV]",100,0,10);

  fCRTReg         = new TH1F("crtReg",        "CRT Region; region ID",23,29,52);
  fCRTRegAV       = new TH1F("crtRegAV",      "CRT Region; region ID",23,29,52);
  fCRTRegAVNotFV  = new TH1F("crtRegAVNotFV" ,"CRT Region; region ID",23,29,52);
  fCRTRegFV       = new TH1F("crtRegFV",      "CRT Region; region ID",23,29,52);
  fNCRTReg        = new TH1F("nCRTReg",       "No. CRT Regions Hit; no. regions",4,0,4);
  fNCRTRegAV      = new TH1F("nCRTRegAV",     "No. CRT Regions Hit; no. regions",4,0,4);
  fNCRTRegAVNotFV = new TH1F("nCRTRegAVNotFV","No. CRT Regions Hit; no. regions",4,0,4);
  fNCRTRegFV      = new TH1F("nCRTRegFV",     "No. CRT Regions Hit; no. regions",4,0,4);

  // Load configuration parameters
  fTruthTag = { "generator" };
  fDetTag = { "crtdaq" };


  if (config) {
    fhicl::ParameterSet pconfig = \
      config->get<fhicl::ParameterSet>("AutoVetoAnalysis");

    fFidXOut = pconfig.get<double>("FiducialXOuter", 25.0);
    fFidXIn  = pconfig.get<double>("FiducialXInner", 0.0);
    fFidYTop = pconfig.get<double>("FiducialYTop", 25.0);
    fFidYBot = pconfig.get<double>("FiducialYBottom",25.0);
    fFidZUp  = pconfig.get<double>("FiducialZUpstream",30.0);
    fFidZDown= pconfig.get<double>("FiducialZDownstream",50.0); 

    fTruthTag = \
      { pconfig.get<std::string>("MCTruthTag", "generator") };
  
    fDetTag = 
      { pconfig.get<std::string>("DetSimTag", "crtdaq") };
}

  // Add custom branches
  AddBranch("nuCount",           &fNuCount);
  AddBranch("crtTag",            &fCRTTag);
  AddBranch("ccnc",              &fCCNC);
  AddBranch("nuReg",             &fNuReg);

  // Use LArSoft service functions via the ProviderManager
  if (fProviderManager) {
    std::cout << "Detector: "
              << fProviderManager->GetGeometryProvider()->DetectorName()
              << std::endl;
  }
}//end Initialize


void AutoVetoSelection::Finalize() {
  // Output our histograms to the ROOT file
  fOutputFile->cd();
  fNuVertexXZHist->Write();
  fNuVertexXYHist->Write();
  fNuVertexYZHist->Write();
  fNuVertexXZAVHist->Write();
  fNuVertexXYAVHist->Write();
  fNuVertexYZAVHist->Write();
  fNuVertexXZFVHist->Write();
  fNuVertexXYFVHist->Write();
  fNuVertexYZFVHist->Write();

  fNuE->Write(); 
  fNuEAV->Write(); 
  fNuEAVNotFV->Write(); 
  fNuEFV->Write(); 
  fNuEVeto->Write(); 
  fNuEVetoAV->Write(); 
  fNuEVetoAVNotFV->Write(); 
  fNuEVetoFV->Write();

  /*fCRTReg->Write();
  fCRTRegAV->Write();
  fCRTRegAVNotFV->Write();
  fCRTRegFV->Write();
  fNCRTReg->Write();
  fNCRTRegAV->Write();
  fNCRTRegAVNotFV->Write();
  fNCRTRegFV->Write();*/
}


bool AutoVetoSelection::ProcessEvent(
    const gallery::Event& ev,
    const std::vector<event::Interaction>& truth,
    std::vector<event::RecoInteraction>& reco) {

  if (fEventCounter % 10 == 0) {
    std::cout << "AutoVeto: Processing event "
              << fEventCounter << std::endl;
  }
  fEventCounter++;

  const TVector3 fid1 = {fFidXOut, fFidYBot, fFidZUp};
  const TVector3 fid2 = {fFidXIn, fFidYTop, fFidZDown};

  fNuCount = 0;
  fCRTTag = 0;
  fNuReg = -1;
  fCCNC = -1;
  //std::set<int> crtregs;

  // Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);

  auto const &crtsim = 
    *ev.getValidHandle<std::vector<icarus::crt::CRTData>>(fDetTag);

  // Fill in the custom branches
  if (mctruths.size()>1)
      std::cout << "more than 1 neutrino found!" << std::endl;

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {

    auto const& mctruth = mctruths.at(i);

    // Fill neutrino vertex position histogram
    if (mctruth.NeutrinoSet()&&abs(mctruth.GetNeutrino().Nu().PdgCode())==14) {
      fCCNC = mctruth.GetNeutrino().CCNC();
      const TVector3 point = mctruth.GetNeutrino().Nu().Position(0).Vect();

      if(IsAV(*fProviderManager,point)) {
        if(IsFV(*fProviderManager,point,fid1,fid2)) {
            fNuReg=2;
        }
        else {
            fNuReg=1;
        }
      }

      double nue = mctruth.GetNeutrino().Nu().E();
      fNuCount++;
      if(!crtsim.empty()) {
          fCRTTag++;
          //for(const icarus::crt::CRTData& dat : crtsim) {
          //    crtregs.insert(MacToADReg(dat.Mac5()));
          //}
      }
      else {
          fNuEVeto->Fill(nue);
      }

      fNuVertexXZHist->Fill(point.Z(),point.X());
      fNuVertexXYHist->Fill(point.X(),point.Y());
      fNuVertexYZHist->Fill(point.Z(),point.Y());
      fNuE->Fill(nue);
      //if(crtregs.empty()) fCRTReg->Fill(0);
      //else fCRTReg->Fill(*(crtregs.begin()));
      //fNCRTReg->Fill(crtregs.size());

      if(fNuReg==1||fNuReg==2){ //AV, including FV
          fNuVertexXZAVHist->Fill(point.Z(),point.X());
          fNuVertexXYAVHist->Fill(point.X(),point.Y());
          fNuVertexYZAVHist->Fill(point.Z(),point.Y());
          fNuEAV->Fill(nue);
          if(fCRTTag==0) {
              fNuEVetoAV->Fill(nue);
              //fCRTRegAV->Fill(0);
          }
          //else
          //    fCRTRegAV->Fill(*(crtregs.begin()));
          //fNCRTRegAV->Fill(crtregs.size());
      }
      if(fNuReg==1) { //AV and NOT FV
          fNuEAVNotFV->Fill(nue);
          if(fCRTTag==0) {
              fNuEVetoAVNotFV->Fill(nue);
              //fCRTRegAVNotFV->Fill(0);
          }
          //else
          //    fCRTRegAVNotFV->Fill(*(crtregs.begin()));
          //fNCRTRegAVNotFV->Fill(crtregs.size());
      }
      if(fNuReg==2) { //FV
          fNuVertexXZFVHist->Fill(point.Z(),point.X());
          fNuVertexXYFVHist->Fill(point.X(),point.Y());
          fNuVertexYZFVHist->Fill(point.Z(),point.Y());
          fNuEFV->Fill(nue);
          if(fCRTTag==0) {
              fNuEVetoFV->Fill(nue);
              //fCRTRegFV->Fill(0);
          }
          //else
          //    fCRTRegFV->Fill(*(crtregs.begin()));
          //fNCRTRegFV->Fill(crtregs.size());
      }
    }//if we have a nuetrino and it's a nu_mu

    // Add in the "reconstructed" interaction
    //
    // Construct truth information from the provided vector
    //event::RecoInteraction interaction(i);

    // Get truth information from the event
    //TVector3 lepton_momentum = mctruth.GetNeutrino().Lepton().Momentum().Vect();

    // The "truth" vector also collects truth information for use
    // The i-th truth Interaction is the same as the i-th entry in the mctruths vector
    //double lepton_energy = truth[i].lepton.energy;

    // get "reconstructed" energy
    //interaction.reco_energy = util::ECCQE(lepton_momentum, lepton_energy);
    // Save the reconstructed interaction
    //reco.push_back(interaction);

  }//end loop over MCTruths

  return true;
}

  }  // namespace AutoVetoAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::AutoVetoAnalysis::AutoVetoSelection)

