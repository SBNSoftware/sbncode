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
    fCuts.Initialize(config->get<double>("trackMatchCompletionCut", -1));
    fCRTHitDistance = config->get<double>("CRTHitDistance", -1);
    fRecoEvent = NULL;
  }

  void Selection::FileSetup(TTree *eventTree) {
    eventTree->SetBranchAddress("reco_event", &fRecoEvent);
  }

  void Selection::ProcessEvent(const event::Event *core_event) {
    for (unsigned i = 0; i < fRecoEvent->reco.size(); i++) {
      // std::array<bool, Cuts::nCuts> cuts = fCuts.ProcessRecoCuts(*fRecoEvent, i);
       std::array<bool, Cuts::nCuts> cuts {false};
      fHistograms.Fill(*fRecoEvent, *core_event, cuts);
    }
  }

  void Selection::Finalize() {
    fOutputFile->cd();
    fHistograms.Write();
  }
  }   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Selection);

