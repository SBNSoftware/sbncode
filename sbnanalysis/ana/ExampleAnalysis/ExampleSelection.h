#ifndef __sbnanalysis_ana_ExampleAnalysis_ExampleProcessor__
#define __sbnanalysis_ana_ExampleAnalysis_ExampleProcessor__

/**
 * An example event selection processor.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "io/SelectionBase.hh"

class TH1F;
class TH2D;

namespace ana {
  namespace ExampleAnalysis {

class ExampleSelection : public io::SelectionBase {
public:
  ExampleSelection();
  virtual ~ExampleSelection();
  void ProcessEvent(gallery::Event& ev);

  int fMyVar;
  std::vector<double> fMyVector;

protected:
  art::InputTag mctruths_tag;
  TH1F* n_nu_hist;
  TH1F* nu_pdg_hist;
  TH2D* nu_vtx_YZ_hist;
  TH2D* nu_vtx_XZ_hist;
};

  }  // namespace ExampleAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_ExampleAnalysis_ExampleProcessor__

