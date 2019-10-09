#include "Normalize.h"
#include "core/SubRun.hh"


void ana::SBNOsc::Normalize::Initialize(const fhicl::ParameterSet &fcl) {
  fPotPerReadout = fcl.get<double>("PotPerReadout");

  fLastCosmicFileID = -1;
  fLastCosmicEventID = -1;

  fLastNeutrinoFileID = -1;
  fLastNeutrinoEventID = -1;

  fNeutrinoPOT = 0.;
  fNNeutrinoEvents = 0;
  fNCosmicEvents = 0;

  // TODO: get number of cosmic events pre-filtering. For now, 
  // rely on knowing number of cosmic events per file. This
  // number is hard-coded, since we will want to change
  // how we do this long term
  fCosmicEventsPerFile = 666;

}

void ana::SBNOsc::Normalize::AddCosmicEvent(const event::Event &event) {
  // new file -- increment the total
  if (fLastCosmicFileID == -1 || event.metadata.fileEntry != fLastCosmicFileID) {
    fNCosmicEvents += fCosmicEventsPerFile;
  }
  fLastCosmicFileID = event.metadata.fileEntry;
  fLastCosmicEventID = event.metadata.eventID;
}

void ana::SBNOsc::Normalize::AddNeutrinoEvent(const event::Event &event) {
  // new file -- just count the event
  if (fLastNeutrinoFileID == -1 || event.metadata.fileEntry != fLastNeutrinoFileID) {
    fNNeutrinoEvents += 1;
  }
  else {
    fNNeutrinoEvents += event.metadata.eventID - fLastNeutrinoEventID;
  }
  fLastNeutrinoFileID = event.metadata.fileEntry;
  fLastNeutrinoEventID = event.metadata.eventID;
}

void ana::SBNOsc::Normalize::AddNeutrinoSubRun(const SubRun &subrun) {
  fNeutrinoPOT += subrun.totgoodpot;
}

double ana::SBNOsc::Normalize::ScaleCosmic(double goal_pot) const {
  // fraction of all POT accounted for by the neutrino sample
  double pot_per_neutrino = fNeutrinoPOT / fNNeutrinoEvents;

  // POT per cosmic event
  double pot_per_cosmic = fPotPerReadout - pot_per_neutrino; 

  // total POT equaivalent of cosmic
  double cosmic_pot = fNCosmicEvents * pot_per_cosmic;
  
  // scale
  return goal_pot / cosmic_pot;
}

double ana::SBNOsc::Normalize::ScaleNeutrino(double goal_pot) const {
  // simpler -- just scale neutrino POT to goal
  return goal_pot / fNeutrinoPOT;
}



