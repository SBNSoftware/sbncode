#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

#include <TFile.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <core/Event.hh>

#include "Selection.h"

namespace ana {
namespace SBNOsc {
  void Selection::Initialize(fhicl::ParameterSet *config) {
    fOutputFile = new TFile("output.root", "CREATE");
    fOutputFile->cd();
    fCRTHitDistance = config->get<double>("CRTHitDistance", -1);
    fCutConfig = config->get<fhicl::ParameterSet>("Cuts");
    fRecoEvent = NULL;
  }

  void Selection::FileSetup(TTree *eventTree) {
    eventTree->SetBranchAddress("reco_event", &fRecoEvent);
    fCuts.Initialize(fCutConfig, fProviderManager->GetGeometryProvider());
  }

  void Selection::ProcessEvent(const event::Event *core_event) {
    fHistograms.Fill(*fRecoEvent, *core_event, fCuts);
  }

  void Selection::Finalize() {
    fOutputFile->cd();
    fHistograms.Write();
  }
  }   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Selection);

