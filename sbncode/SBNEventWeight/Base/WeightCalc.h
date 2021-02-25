#ifndef _SBN_WEIGHTCALC_H_
#define _SBN_WEIGHTCALC_H_

// Adapted from LArSoft's larsim EventWeight

#include <cassert>
#include <string>
#include <map>
#include <vector>
#include "art/Framework/Principal/fwd.h"
#include "TMatrixD.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h"

namespace fhicl { class ParameterSet; }

namespace CLHEP { class HepRandomEngine; }

namespace sbn {
  namespace evwgh {

class WeightCalc {
public:
  virtual void Configure(fhicl::ParameterSet const& pset,
                         CLHEP::HepRandomEngine&) = 0;

  virtual std::vector<float> GetWeight(art::Event& e, size_t inu) = 0;

  void SetName(std::string name) { fName = name; }
  void SetType(std::string type) { fType = type; }

  std::string GetName() { return fName; }
  std::string GetType() { return fType; }

  EventWeightParameterSet fParameterSet;

private:
  std::string fName;
  std::string fType;
};

  }  // namespace evwgh
}  // namespace sbn

#endif  // _SBN_WEIGHTCALC_H_

