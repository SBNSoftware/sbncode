#include <iostream>
#include <vector>
#include <TH1F.h>
#include <TH2D.h>
#include <TRandom.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ExampleSelection.h"

namespace ana {
  namespace ExampleAnalysis {

ExampleSelection::ExampleSelection() : SelectionBase(), fMyVar(0) {
  // Here you name the thing that "produced" the data product that you want to
  // look at. In our event dump we see two things:
  //
  // GenieGen.. | corsika..... | .... | std::vector<simb::MCTruth>.... | ....1
  // GenieGen.. | generator... | .... | std::vector<simb::MCTruth>.... | ....1
  //
  // This means that if you want to look at the simb::MCTruth data product
  // there are two "producers." We want to look at neutrinos so we choose
  // "generator" if you wanted to look at cosmics you could pick "coriska."
  mctruths_tag = { "generator" };
 
  // Define a few histograms
  n_nu_hist = new TH1F("Nnu", ";Neutrino count;Events", 50, 0, 50);
  nu_pdg_hist = new TH1F("pdg", ";PDG code;Events", 32, -16, 16);
  nu_vtx_YZ_hist = new TH2D("nu_vtx_YZ", "",
                            100, -1000, 1000, 100, -1000, 1000);
  nu_vtx_XZ_hist = new TH2D("nu_vtx_XZ", "",
                            100, -1000, 1000, 100, -1000, 1000);
}


ExampleSelection::~ExampleSelection() {
  // Write the histograms to the output file
  fOutputFile->cd();

  nu_vtx_XZ_hist->Write();
  nu_vtx_YZ_hist->Write();
  nu_pdg_hist->Write(); 
  n_nu_hist->Write(); 
}


void ExampleSelection::ProcessEvent(gallery::Event& ev) {
  // Grab the data product that you want from the event
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(mctruths_tag);

  // Dump how many MC truth entries there are. This is the number of neutrino
  // interactions there are in the event. Knowing what volume you are
  // generating your events in is important! 
  n_nu_hist->Fill(mctruths.size());
  
  // Now we'll iterate through these 
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);

    // Now for each simb::MCTruth we look at two things      
    nu_pdg_hist->Fill(mctruth.GetNeutrino().Nu().PdgCode());

    nu_vtx_YZ_hist->Fill(mctruth.GetNeutrino().Nu().Vy(),
                         mctruth.GetNeutrino().Nu().Vz());

    nu_vtx_XZ_hist->Fill(mctruth.GetNeutrino().Nu().Vx(),
                         mctruth.GetNeutrino().Nu().Vz());
  }

  // Fill in custom variables
  fMyVar++;
  fMyVector = { gRandom->Gaus(), gRandom->Gaus() };
}


  }  // namespace ExampleAnalysis
}  // namespace ana

