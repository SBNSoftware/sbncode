#include "CAFAna/Systs/UniverseOracle.h"
#include "CAFAna/Core/Utilities.h"
using namespace ana;

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"

#include <iostream>

void test_oracle()
{
  UniverseOracle oracle;

  const unsigned int N = oracle.Systs().size();
  unsigned int Nx = sqrt(N);
  unsigned int Ny = sqrt(N);
  if(Nx*Ny < N) ++Nx;
  if(Nx*Ny < N) ++Ny;

  TCanvas* c = new TCanvas;
  c->Divide(Nx, Ny, 0, 0);

  int i = 0;
  for(const std::string& name: oracle.Systs()){
    const std::vector<double>& shifts = oracle.ShiftsForSyst(name);
    std::cout << name << ": " << shifts.size() << " universes" << std::endl;
    for(int sigma = -3; sigma <= +3; ++sigma){
      if(sigma == 0) continue;
      double truth;
      const unsigned int idx = oracle.ClosestIndex(name, sigma, &truth);
      std::cout << "  " << sigma << "sigma -> univ " << idx << " with " << truth << "sigma" << std::endl;
    }

    c->cd(++i);
    TH1* h = new TH1F(UniqueName().c_str(), (name+";Shift (#sigma);Universes").c_str(), 30, -3, +3);
    for(double v: shifts) h->Fill(v);
    h->Draw("hist");
  }

  c->cd(0);
  gPad->Print("univs.pdf");
}
