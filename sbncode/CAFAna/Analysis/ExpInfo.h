#pragma once

namespace ana

{
  const double kPOTnominal = 6.6e20;
  const double kPOTuBoone = 1.3e21;
  const double kBaselineSBND = 0.11;
  const double kBaselineMicroBoone = 0.47;
  const double kBaselineIcarus = 0.6; 
  const int kSBND = 0;
  const int kMicroBoone = 1;
  const int kICARUS = 2;
  const std::vector<double> kBLs = {kBaselineSBND, kBaselineMicroBoone, kBaselineIcarus};
}
