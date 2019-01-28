//C++ Includes 
#include <iostream>
#include <vector>

//Root Includes 
#include <TH2D.h>

//Framework Includes 
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"

//Larsoft Includes 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

//SBN Includes 
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

NueSelection::NueSelection() : SelectionBase(), EventCounter(0), NuCount(0) {}


void NueSelection::Initialize(fhicl::ParameterSet* config) {
  // Load configuration parameters
  fTruthTag = { "generator" };

  fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NueSelectione");

  // setup active volume bounding boxes
  std::vector<fhicl::ParameterSet> AVs =				\
    pconfig.get<std::vector<fhicl::ParameterSet> >("active_volumes");
  for (auto const& AV : AVs) {
    double xmin = AV.get<double>("xmin");
    double ymin = AV.get<double>("ymin");
    double zmin = AV.get<double>("zmin");
    double xmax = AV.get<double>("xmax");
    double ymax = AV.get<double>("ymax");
    double zmax = AV.get<double>("zmax");
    fConfig.active_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  
  std::vector<fhicl::ParameterSet> FVs =				\
    pconfig.get<std::vector<fhicl::ParameterSet> >("fiducial_volumes");
  for (auto const& FV : FVs) {
    double xmin = FV.get<double>("xmin");
    double ymin = FV.get<double>("ymin");
    double zmin = FV.get<double>("zmin");
    double xmax = FV.get<double>("xmax");
    double ymax = FV.get<double>("ymax");
    double zmax = FV.get<double>("zmax");
    fConfig.fiducial_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  
  if (config) {
    fTruthTag = { pconfig.get<std::string>("MCTruthTag", "generator") };
  }

  hello();
}


void NueSelection::Finalize() {}


bool NueSelection::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco){
  if (EventCounter % 10 == 0) {
    std::cout << "NueSelection: Processing event " << EventCounter << " "
              << "(" << NuCount << " neutrinos selected)"
              << std::endl;
  }
  EventCounter++;

  //Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);

  //Get tracks and showers 
  auto const& mctracks = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mcshowers = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  //Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    // build the interaction
    Event::Interaction interaction = truth[i];

    //Calculate the Energy 
    std::vector<double> visible_energy = FlavourEnergyDeposition(mctruth, mctracks, mcshowers);
    Event::RecoInteraction reco_interaction(interaction, i);
    reco_interaction.reco_energy = visible_energy[0];

    //Run the Selection
    bool selection = Select(ev, mctruth, i, fConfig, reco_interaction);

    //Calculate Energy - Remove Events below 200 MeV, remove 30% between 200 and 350 and 5% above 350 MeV.
    //Check for Tracks - Remove any track that come from the neutrino if 1m 
    //Find electron vertex candiates - Remove 80% due to fake dEdx cut 
    //If pi0 is made 2 photons come out as candiates in both are above 100 MeV remove.
    //If hadrons produce have more than 50 MeV this is visible if all photons pair produce more than 3cm away from the vertex remove.
    //Remove 94% pion of events as dEdx cut
    //Event could be CC0Pi Muons Interaction with photon background which mimics the CC1Pi Electron interatactions. these need to go throught the cuts above as well.  
    
    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::RecoInteraction interaction(truth[i], i);
      reco.push_back(interaction);
    }
  }

  bool selected = !reco.empty();

  if (selected) {
    NuCount++;
  }

  return selected;
}
  
bool NueSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth,unsigned i, Config fConfig,Event::RecoInteraction reco_interaction){

  return false;

}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)

