/**
 * Run sensitivity calculation.
 *
 * Document...
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "Covariance.h"
#include "Chi2Sensitivity.h"

int main(int argc, char* argv[]) {
  std::vector<ana::SBNOsc::EventSample> samples;
  for (int i=2; i<argc; i++) {
    // Build sample list
  }

  assert(!samples.empty());

  ana::SBNOsc::Covariance cov(samples);

  ana::SBNOsc::Chi2Sensitivity chi2();

  // Write contours...

  return 0;
}

