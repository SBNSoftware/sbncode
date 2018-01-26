#ifndef __sbnanalysis_io_BaseTree__
#define __sbnanalysis_io_BaseTree__

/**
 * The standard minimum output tree.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include <map>
#include <string>
#include <vector>
#include <TTree.h>
#include <TVector3.h>

class Event {
public:
  // Event information
  class Metadata {
  public:
    int run;
    int subrun;
    int eventID;
  };

  // Neutrino information
  class Neutrino {
  public:
    bool ccnc;
    int pdg;
    int targetPDG;
    int intcode;
    double bjorkenX;
    double inelasticityY;
    double q2;
    double w;
    double t;
    double energy;
    TVector3 momentum;
  };

  // Final state hadronic system
  class FinalStateParticle {
  public:
    int pdg;
    double energy;
    TVector3 momentum;
  };

  // Products associated with one neutrino
  class Interaction {
  public:
    Neutrino neutrino;
    FinalStateParticle lepton;
    std::vector<FinalStateParticle> hadrons;
  };

  Metadata metadata;
  std::vector<Interaction> interactions;
  std::map<std::string, std::vector<double> > weights;
};

#endif  // __sbnanalysis_io_BaseTree__

