#ifndef _sbnanalysis_Normalize_hh_
#define _sbnanalysis_Normalize_hh_

#include "fhiclcpp/ParameterSet.h"
#include "core/Event.hh"
#include "core/SubRun.hh"

namespace ana {
namespace SBNOsc {

class Normalize {
public:
  void Initialize(const fhicl::ParameterSet &fcl);

  void AddCosmicEvent(const event::Event &event);
  void AddNeutrinoEvent(const event::Event &event);

  void AddNeutrinoSubRun(const SubRun &subrun);

  double ScaleNeutrino(double goal_pot) const;
  double ScaleCosmic(double goal_pot) const;


private:
  int fLastCosmicFileID;
  int fLastCosmicEventID;

  int fLastNeutrinoFileID;
  int fLastNeutrinoEventID;

  double fPotPerReadout;
  double fNeutrinoPOT;
  unsigned fCosmicEventsPerFile;
  unsigned fNNeutrinoEvents;
  unsigned fNCosmicEvents;
};



} // end namespace ana
} // end namespace SBNOsc
#endif
