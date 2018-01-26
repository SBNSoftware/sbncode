#ifndef __sbnanalysis_ana_ExampleAnalysis_ExampleProcessor__
#define __sbnanalysis_ana_ExampleAnalysis_ExampleProcessor__

/**
 * An example event selection processor.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"

class TH2D;

namespace ana {
  namespace ExampleAnalysis {

class ExampleSelection : public core::SelectionBase {
public:
  ExampleSelection();
  void Initialize(Json::Value* config=NULL);
  void Finalize();
  void ProcessEvent(gallery::Event& ev);

protected:
  art::InputTag mctruths_tag;

  // Configuration parameters
  int fMyParam;

  // Custom data branches
  int fNuCount;
  int fMyVar;

  // Histograms
  TH2D* fNuVertexXZHist;
};

  }  // namespace ExampleAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_ExampleAnalysis_ExampleProcessor__

