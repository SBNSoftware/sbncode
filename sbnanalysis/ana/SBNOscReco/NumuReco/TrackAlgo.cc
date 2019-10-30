#include "TrackAlgo.h"
#include <numeric>

float numu::MeanTruncateddQdx(const anab::Calorimetry &calo) {
  // copy the dQdx list for partial sorting to find median
  std::vector<float> dQdx = calo.dQdx();

  // if no or only 1 calo points, then no dQdx
  if (dQdx.size() <= 1) return -9999.;

  // calculate the median
  float median = -1;
  if (dQdx.size() % 2 == 0 /* even */) {
    const auto median_lo = dQdx.begin() + dQdx.size()/2 - 1;
    const auto median_hi = dQdx.begin() + dQdx.size()/2;
    std::nth_element(dQdx.begin(), median_lo, dQdx.end());
    float median_lo_val = *median_lo;
    std::nth_element(dQdx.begin(), median_hi, dQdx.end());
    float median_hi_val = *median_hi;
    median = (median_lo_val + median_hi_val) / 2.;    
  }
  else /* odd */ {
    const auto median_ptr = dQdx.begin() + dQdx.size()/2;
    std::nth_element(dQdx.begin(), median_ptr, dQdx.end());
    median = *median_ptr;
  }

  // get the mean and variance
  float mean = std::accumulate(dQdx.begin(), dQdx.end(), 0.) / dQdx.size();  
  float var = 0.;
  for (float d: dQdx) {
    var += (d - mean) * (d - mean);
  }
  var = var / dQdx.size();

  // truncated mean
  float trunc_sum = 0.;
  unsigned n_point = 0;
  for (float d: dQdx) {
    if (var < (d - mean) * (d - mean)) {
      trunc_sum += d;
      n_point ++;
    }
  }
  // if no points available, return garbage
  if (n_point == 0) return -9999.;

  return trunc_sum / n_point;
}
