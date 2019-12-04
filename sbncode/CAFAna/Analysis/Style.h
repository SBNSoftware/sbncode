#pragma once

#include "Rtypes.h"
#include "TAttLine.h"

namespace ana
{
  // Standardizing colors here means that if we change our minds in future it
  // should hopefully be easy to propagate the changes to all plots.

  // Colors for histograms

  const Color_t kTotalMCColor = kRed;
  const Color_t kTotalMCErrorBandColor = kRed-10;
  const Color_t kNueSignalColor = kViolet-5;
  const Color_t kTotalBackgroundColor = kAzure+2;
  const Color_t kNCBackgroundColor = kAzure;
  const Color_t kBeamNueBackgroundColor = kPink+9;
  const Color_t kCosmicBackgroundColor = kAzure+1;
  /// This is the color for the numu CC background to the nue analysis. I doubt
  /// numu plots will ever use this, preferring to just show the total MC.
  const Color_t kNumuBackgroundColor = kGreen+2;


  // Colors and styles for contours

  const Color_t kNormalHierarchyColor = kAzure+3;
  const Color_t kInvertedHierarchyColor = kOrange+9;

  const Style_t k90PercentConfidenceStyle = 9;
  const Style_t k68PercentConfidenceStyle = 7; ///< Dashed

  // Colors for slices and 123 sigma contours
  const Color_t kPrimColorNH = kAzure+2;
  const Color_t kSecoColorNH = kAzure+1;

  const Color_t kPrimColorIH = kOrange+9;
  const Color_t kSecoColorIH = kOrange+10;

  const Style_t kFillStyleNH = 3013;
  const Style_t kFillStyleIH = 3003;

  const Style_t k1SigmaConfidenceStyle = 7;
  const Style_t k2SigmaConfidenceStyle = 9;
  const Style_t k3SigmaConfidenceStyle = 10;

  const Color_t kCentralValueColorNH = kAzure+3;
  const Color_t k1SigmaConfidenceColorNH = kAzure+3;
  const Color_t k90PercConfidenceColorNH = kAzure+2;
  const Color_t k2SigmaConfidenceColorNH = kAzure+2;
  const Color_t k3SigmaConfidenceColorNH = kAzure+1;

  const Color_t kCentralValueColorIH = kOrange+9;
  const Color_t k1SigmaConfidenceColorIH = kOrange+9;
  const Color_t k2SigmaConfidenceColorIH = kOrange+10;
  const Color_t k3SigmaConfidenceColorIH = kOrange+8;
  const Color_t k90PercConfidenceColorIH = kOrange+10;

  const Color_t kCentralValueColor= kGray+2;
  const Style_t kCentralValueStyle = kSolid;

  const Float_t kBlessedTitleFontSize = 0.08; //not in pixels
  const Float_t kBlessedLabelFontSize = 0.06; //not in pixels
}
