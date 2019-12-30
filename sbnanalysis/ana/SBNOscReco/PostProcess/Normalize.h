#ifndef _sbnanalysis_Normalize_hh_
#define _sbnanalysis_Normalize_hh_

#include "fhiclcpp/ParameterSet.h"
#include "core/Event.hh"
#include "core/SubRun.hh"
#include "core/FileMeta.hh"

namespace ana {
namespace SBNOsc {

/**
 * Class to handle normalization of overlay and intime-cosmic samples.
 */
class Normalize {
public:
  /**
 * Initialize the class
 * \param fcl The fhicl configuration
 */
  void Initialize(const fhicl::ParameterSet &fcl);

  /**
 * Add a neutrino subrunto the normalization
 * \param subrun The sbncode core Subrun information
 */
  void AddNeutrinoSubRun(const SubRun &subrun);

  void AddCosmicFile(const FileMeta &meta);
  void AddNeutrinoFile(const FileMeta &meta);

  /**
 * Scale the neutrino events to a certain POT. To be called after all events are processed
 * \param goal_pot The POT to scale the sample to.
 *
 * \return The scale factor for each neutrino event to the set POT
 */
  double ScaleNeutrino(double goal_pot) const;
  /**
 * Scale the cosmic events to a certain POT. To be called after all events are processed
 * \param goal_pot The POT to scale the sample to.
 *
 * \return The scale factor for each cosmic event to the set POT
 */
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
