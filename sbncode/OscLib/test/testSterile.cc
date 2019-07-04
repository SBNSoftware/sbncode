#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "OscLib/OscCalculatorSterile.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"

#include <iostream>
#include <fenv.h>

int main(int argc, char* argv[])
{
  feenableexcept(FE_INVALID); // Spot any infs or nans early

  osc::OscCalculatorSterile oscStd;
  osc::OscCalculatorSterile oscSterile;

  double baselineNear(1.04);
  double baselineFar(810.0);
  double eDensity(2*1.45);
  
  double th12(0.6), th13(0.159), th23(0.681);
  
  double dm221(0.0000759), dm232(0.00239);

  double dm243 = std::atof( argv[1] );
  double th14  = std::atof( argv[2] );
  double th24  = std::atof( argv[3] );
  double th34  = std::atof( argv[4] );

  //double th14(0.1), th24(0.1), th34(0.1);
  //double dm243(1.0);

  double delta1(0.0), delta2(0.0), delta3(0.0);

  TGraph* gBeamNueDisND = new TGraph();
  TGraph* gNueAppND     = new TGraph();
  TGraph* gNueAppFD     = new TGraph();
  TGraph* gNumuDisND    = new TGraph();
  TGraph* gNumuDisFD    = new TGraph();
  TGraph* gNCDisND      = new TGraph();
  TGraph* gNCDisFD      = new TGraph();

  TGraph* g4FlvBeamNueDisND = new TGraph();
  TGraph* g4FlvNueAppND     = new TGraph();
  TGraph* g4FlvNueAppFD     = new TGraph();
  TGraph* g4FlvNumuDisND    = new TGraph();
  TGraph* g4FlvNumuDisFD    = new TGraph();
  TGraph* g4FlvNCDisND      = new TGraph();
  TGraph* g4FlvNCDisFD      = new TGraph();

  TGraph* g3FlvBeamNueDisND = new TGraph();
  TGraph* g3FlvNueAppND     = new TGraph();
  TGraph* g3FlvNueAppFD     = new TGraph();
  TGraph* g3FlvNumuDisND    = new TGraph();
  TGraph* g3FlvNumuDisFD    = new TGraph();
  TGraph* g3FlvNCDisND      = new TGraph();
  TGraph* g3FlvNCDisFD      = new TGraph();

  TGraph* gNue2NusND = new TGraph();
  
  oscStd.SetNFlavors(3);
  oscStd.SetRho( eDensity );
  oscStd.SetDm(2, dm221);
  oscStd.SetDm(3, dm232 + dm221);
  oscStd.SetAngle(1, 2, th12);
  oscStd.SetAngle(1, 3, th13);
  oscStd.SetAngle(2, 3, th23);
  oscStd.SetDelta(1, 3, delta1);

  oscSterile.SetNFlavors(4);
  oscSterile.SetL(baselineFar);
  oscSterile.SetRho( eDensity );
  oscSterile.SetDm(2, dm221);
  oscSterile.SetDm(3, dm232 + dm221);
  oscSterile.SetDm(4, dm232 + dm221 + dm243);
  oscSterile.SetAngle(1, 2, th12);
  oscSterile.SetAngle(1, 3, th13);
  oscSterile.SetAngle(2, 3, th23);
  oscSterile.SetAngle(1, 4, th14);
  oscSterile.SetAngle(2, 4, th24);
  oscSterile.SetAngle(3, 4, th34);
  oscSterile.SetDelta(1, 3, delta1);
  oscSterile.SetDelta(1, 4, delta2);
  oscSterile.SetDelta(2, 4, delta3);

  osc::OscCalculatorSterile oscTest;
  std::vector<double> state = oscSterile.GetState();
  oscTest.SetState(state);

  //Testing the get and set state functionality
  std::cout << "NFlavors: sterile= " << oscSterile.GetNFlavors() << " test= " << oscTest.GetNFlavors() << std::endl;
  std::cout << "L:        sterile= " << oscSterile.GetL()        << " test= " << oscTest.GetL()        << std::endl;
  std::cout << "Rho:      sterile= " << oscSterile.GetRho()      << " test= " << oscTest.GetRho()      << std::endl;
  std::cout << "Dm13:     sterile= " << oscSterile.GetDm(3)      << " test= " << oscTest.GetDm(3)      << std::endl;
  std::cout << "Dm14:     sterile= " << oscSterile.GetDm(4)      << " test= " << oscTest.GetDm(4)      << std::endl;
  std::cout << "Th23:     sterile= " << oscSterile.GetAngle(2,3) << " test= " << oscTest.GetAngle(2,3) << std::endl;
  std::cout << "Th24:     sterile= " << oscSterile.GetAngle(2,4) << " test= " << oscTest.GetAngle(2,4) << std::endl;
  std::cout << "Delta13:  sterile= " << oscSterile.GetDelta(1,3) << " test= " << oscTest.GetDelta(1,3) << std::endl;

  int iPoint(0);
  for (double energy = 0.5; energy < 10.0; energy += 0.0001)
    {
      //Beam NuE
      oscStd.SetL(baselineNear);
      oscSterile.SetL(baselineNear);

      g3FlvBeamNueDisND->SetPoint(iPoint, energy, oscStd.P(12, 12, energy));
      g4FlvBeamNueDisND->SetPoint(iPoint, energy, oscSterile.P(12, 12, energy));
      gBeamNueDisND->SetPoint(iPoint, energy, oscSterile.P(12, 12, energy)/oscStd.P(12, 12, energy));
      gNue2NusND->SetPoint(iPoint, energy, oscSterile.P(12, 0, energy));

      //NuMu: ND
      g4FlvNueAppND->SetPoint(iPoint, energy, oscSterile.P(14, 12, energy));
      g3FlvNueAppND->SetPoint(iPoint, energy, oscStd.P(14, 12, energy));
      gNueAppND->SetPoint(iPoint, energy, oscSterile.P(14, 12, energy)/oscStd.P(14, 12, energy));

      g4FlvNumuDisND->SetPoint(iPoint, energy, oscSterile.P(14, 14, energy));
      g3FlvNumuDisND->SetPoint(iPoint, energy, oscStd.P(14, 14, energy));
      gNumuDisND->SetPoint(iPoint, energy, oscSterile.P(14, 14, energy)/oscStd.P(14, 14, energy));

      g4FlvNCDisND->SetPoint(iPoint, energy, 1.0 - oscSterile.P(14, 0, energy));
      g3FlvNCDisND->SetPoint(iPoint, energy, 1.0 - oscStd.P(14, 0, energy));
      gNCDisND->SetPoint(iPoint, energy, 1.0 - oscSterile.P(14, 0, energy));

      //NuMu: FD
      oscStd.SetL(baselineFar);
      oscSterile.SetL(baselineFar);

      g4FlvNueAppFD->SetPoint(iPoint, energy, oscSterile.P(14, 12, energy));
      g3FlvNueAppFD->SetPoint(iPoint, energy, oscStd.P(14, 12, energy));
      gNueAppFD->SetPoint(iPoint, energy, oscSterile.P(14, 12, energy)/oscStd.P(14, 12, energy));

      g4FlvNumuDisFD->SetPoint(iPoint, energy, oscSterile.P(14, 14, energy));
      g3FlvNumuDisFD->SetPoint(iPoint, energy, oscStd.P(14, 14, energy));
      gNumuDisFD->SetPoint(iPoint, energy, oscSterile.P(14, 14, energy)/oscStd.P(14, 14, energy));

      g4FlvNCDisFD->SetPoint(iPoint, energy, 1.0 - oscSterile.P(14, 0, energy));
      g3FlvNCDisFD->SetPoint(iPoint, energy, 1.0 - oscStd.P(14, 0, energy));
      gNCDisFD->SetPoint(iPoint, energy, 1.0 - oscSterile.P(14, 0, energy));
      
      iPoint += 1;
    }

  TFile* file = new TFile("oscSterile.root", "recreate");
  file->cd();

  g4FlvBeamNueDisND->Write("4FlvBeamNueDisND");
  g4FlvNueAppND->Write("4FlvNueAppND");
  g4FlvNumuDisND->Write("4FlvNumuDisND");
  g4FlvNCDisND->Write("4FlvNCDisND");
  g4FlvNueAppFD->Write("4FlvNueAppFD");
  g4FlvNumuDisFD->Write("4FlvNumuDisFD");
  g4FlvNCDisFD->Write("4FlvNCDisFD");

  g3FlvBeamNueDisND->Write("3FlvBeamNueDisND");
  g3FlvNueAppND->Write("3FlvNueAppND");
  g3FlvNumuDisND->Write("3FlvNumuDisND");
  g3FlvNCDisND->Write("3FlvNCDisND");
  g3FlvNueAppFD->Write("3FlvNueAppFD");
  g3FlvNumuDisFD->Write("3FlvNumuDisFD");
  g3FlvNCDisFD->Write("3FlvNCDisFD");

  gBeamNueDisND->Write("BeamNueDisND");
  gNueAppND->Write("NueAppND");
  gNumuDisND->Write("NumuDisND");
  gNCDisND->Write("NCDisND");
  gNueAppFD->Write("NueAppFD");
  gNumuDisFD->Write("NumuDisFD");
  gNCDisFD->Write("NCDisFD");

  gNue2NusND->Write("Nue2NusND");

  return 0;
}
