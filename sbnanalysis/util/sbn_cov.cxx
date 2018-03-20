#include <string>
#include <vector>
#include <Covariance.h>

int main(int argc, char* argv[]) {
  util::Covariance cov;
  cov.AddInputFile("sbnd/output_ccnum.root");
  cov.AddInputFile("sbnd/output_ccnue.root");
  cov.AddInputFile("uboone/output_ccnum.root");
  cov.AddInputFile("uboone/output_ccnue.root");
  cov.AddInputFile("icarus/output_ccnum.root");
  cov.AddInputFile("icarus/output_ccnue.root");
  //cov.AddWeight("genie_qema_Genie");

  cov.AddWeight("kminus_PrimaryHadronNormalization");
  cov.AddWeight("kplus_PrimaryHadronFeynmanScaling");
  cov.AddWeight("kzero_PrimaryHadronSanfordWang");
  cov.AddWeight("nucleoninexsec_FluxUnisim");
  cov.AddWeight("nucleonqexsec_FluxUnisim");
  cov.AddWeight("nucleontotxsec_FluxUnisim");
  cov.AddWeight("piminus_PrimaryHadronSWCentralSplineVariation");
  cov.AddWeight("pioninexsec_FluxUnisim");
  cov.AddWeight("pionqexsec_FluxUnisim");
  cov.AddWeight("piontotxsec_FluxUnisim");
  cov.AddWeight("piplus_PrimaryHadronSWCentralSplineVariation");
  cov.AddWeight("expskin_FluxUnisim");
  cov.AddWeight("horncurrent_FluxUnisim");

  cov.SetOutputFile("./out.root");
  cov.init();
  cov.analyze();

  return 0;
}

