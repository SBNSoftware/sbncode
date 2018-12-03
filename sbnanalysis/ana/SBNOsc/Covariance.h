#ifndef __sbnanalysis_ana_SBNOsc_Covariance__
#define __sbnanalysis_ana_SBNOsc_Covariance__

/**
 * \file Covariance.h
 */

#include <string>
#include <vector>

class TTree;

namespace ana {
  namespace SBNOsc {

class EventSample {
public:
  /** Constructor. */
  EventSample(TTree* _tree, float scaleFactor)
    : tree(_tree), fScaleFactor(scaleFactor) {}

  EventSample(std::vector<std::string> filenames, float fScaleFactor);

  TTree* tree;  //!< Event tree
  float fScaleFactor;  //!< Factor for POT (etc.) scaling
};


class Covariance {
public:
  Covariance(std::vector<EventSample> samples);
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Covariance__

