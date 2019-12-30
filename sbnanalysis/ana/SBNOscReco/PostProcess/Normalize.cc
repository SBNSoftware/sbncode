#include "Normalize.h"
#include "core/SubRun.hh"


void ana::SBNOsc::Normalize::Initialize(const fhicl::ParameterSet &fcl) {
  fPotPerReadout = fcl.get<double>("PotPerReadout", 5.e12);

  fLastCosmicFileID = -1;
  fLastCosmicEventID = -1;

  fLastNeutrinoFileID = -1;
  fLastNeutrinoEventID = -1;

  fNeutrinoPOT = 0.;
  fNNeutrinoEvents = 0;
  fNCosmicEvents = 0;
}

void ana::SBNOsc::Normalize::AddNeutrinoFile(const FileMeta &meta) {
  fNNeutrinoEvents += meta.n_gen_events;
}

void ana::SBNOsc::Normalize::AddCosmicFile(const FileMeta &meta) {
  fNCosmicEvents += meta.n_gen_events;  
}

void ana::SBNOsc::Normalize::AddNeutrinoSubRun(const SubRun &subrun) {
  std::cout << "New subrun: " << subrun.totgoodpot << std::endl;
  fNeutrinoPOT += subrun.totgoodpot;
}

double ana::SBNOsc::Normalize::ScaleCosmic(double goal_pot) const {

  std::cout << "Total N cosmics: " << fNCosmicEvents << std::endl;
  std::cout << "Total N neutrinos: " << fNNeutrinoEvents << std::endl;
  std::cout << "Total neutrino POT: " << fNeutrinoPOT << std::endl;

  // fraction of all POT accounted for by the neutrino sample
  double neutrino_per_pot = 0.;
  if (fNeutrinoPOT > 1e-3 && fNNeutrinoEvents > 0) {
    neutrino_per_pot = fNNeutrinoEvents / (fNeutrinoPOT);
  }

  std::cout << "Neutrino per POT: " << neutrino_per_pot << std::endl;

  // POT per cosmic event
  double cosmic_per_pot = (1./fPotPerReadout) - neutrino_per_pot;

  std::cout << "Cosmic per POT: " << cosmic_per_pot << std::endl;

  // total POT equaivalent of cosmic
  double cosmic_pot = fNCosmicEvents / cosmic_per_pot;

  std::cout << "Cosmic POT: " << cosmic_pot << std::endl;
  
  // scale
  return goal_pot / cosmic_pot;
}

double ana::SBNOsc::Normalize::ScaleNeutrino(double goal_pot) const {
  if (fNeutrinoPOT < 1e-3 || fNNeutrinoEvents == 0) return 0;
  // simpler -- just scale neutrino POT to goal
  return goal_pot / (fNeutrinoPOT);
}



