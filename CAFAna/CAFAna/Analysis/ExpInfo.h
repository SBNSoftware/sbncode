#pragma once

namespace ana

{
  const double POTnominal = 6.6e20;
  const double BaselineSBND = 0.11;
  const double BaselineMicroBoone = 0.47;
  const double BaselineIcarus = 0.6; 
  const int kSBND = 0;
  const int kMicroBoone = 1;
  const int kICARUS = 2;
  const std::vector<double> kBLs = {BaselineSBND, BaselineMicroBoone, BaselineIcarus};
}
