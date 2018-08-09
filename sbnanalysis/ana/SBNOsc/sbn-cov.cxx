/**
 * Generate covariance matrices.
 *
 * Document...
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "Covariance.h"

int main(int argc, char* argv[]) {
  std::vector<ana::SBNOsc::EventSample> samples;
  for (int i=2; i<argc; i++) {
    // Build sample list
  }

  assert(!samples.empty());

  ana::SBNOsc::Covariance cov(samples);

  // Write matrix out to a ROOT file...

  return 0;
}

