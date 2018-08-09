#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>
#include "Covariance.h"

namespace ana {
  namespace SBNOsc {

EventSample::EventSample(std::vector<std::string> filenames, float scaleFactor)
    : fScaleFactor(scaleFactor) {
  TChain* t = new TChain("sbnana");
  for (size_t i=0; i<filenames.size(); i++) {
    t->Add(filenames[i].c_str());
  }

  tree = dynamic_cast<TTree*>(t);
};


Covariance::Covariance(std::vector<EventSample> samples) {
  // Calculate covariance matrix here!
}

  }  // namespace SBNOsc
}  // namespace ana

