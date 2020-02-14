#include <iostream>
#include "OscLib/EarthModel.h"
#include "OscLib/PMNS.h"
#include "TH1F.h"
#include "TFile.h"
using namespace osc;

// Fractional tolerance to model Earth density
double kTOL = 0.02;

void dumpDensityTable() 
{
  EarthModel prem  ("PREM",   kTOL);
  EarthModel stacey("STACEY", kTOL);
  TH1F* h1 = new TH1F("rho_prem",
		      "PREM; R [km]; #rho [g/cm^{3}]",
		      7000,0.0,7000);
  TH1F* h2 = new TH1F("ne_prem",
		      "PREM; R [km]; N_{e} [mole/cm^{3}]",
		      7000,0.0,7000);
  TH1F* h3 = new TH1F("rho_stacey",
		      "Stacey; R [km]; #rho [g/cm^{3}]",
		      7000,0.0,7000);
  TH1F* h4 = new TH1F("ne_stacey",
		      "Stacey; R [km]; N_{e} [mole/cm^{3}]",
		      7000,0.0,7000);
  double r;
  for (int i=0; i<7000; ++i) {
    r = (float)i+0.5;
    h1->Fill(r, prem.  Density(r));
    h2->Fill(r, prem.  Ne(r));
    h3->Fill(r, stacey.Density(r));
    h4->Fill(r, stacey.Ne(r));
  }
}

//......................................................................

void dumpLayers() 
{
  EarthModel prem  ("PREM",  kTOL);
  EarthModel stacey("STACEY",kTOL);

  std::vector<double> rlo;
  std::vector<double> rhi;
  std::vector<double> ne;
  
  TH1F* h1 = new TH1F("ne_prem_layers",
		      "Shells; R [km]; N_{e} [g/cm^{3}]",
		      7000,0.0,7000);
  TH1F* h2 = new TH1F("ne_stacey_layers",
		      "Shells; R [km]; N_{e} [g/cm^{3}]",
		      7000,0.0,7000);
  unsigned int n;
  double r;

  n = 0;
  rlo.clear();
  rhi.clear();
  ne. clear();
  prem.GetLayers(rlo, rhi, ne);
  for (int i=0; i<7000; ++i) {
    r = (float)i+0.5;
    if (r>rhi[n]) ++n;
    if (n<ne.size()) h1->Fill(r, ne[n]);
    else             h1->Fill(r, 0.0);
  }
  
  n = 0;
  rlo.clear();
  rhi.clear();
  ne. clear();
  stacey.GetLayers(rlo, rhi, ne);
  for (int i=0; i<7000; ++i) {
    r = (float)i+0.5;
    if (r>rhi[n]) ++n;
    if (n<ne.size()) h2->Fill(r, ne[n]);
    else             h2->Fill(r, 0.0);
  }
}

//......................................................................

void dumpLineProfiles() 
{
  EarthModel earth("PREM",kTOL);

  double prodh = 20.0;
  double rdet  = -0.3;
  double cosQ;

  std::list<double> L;
  std::list<double> Ne;

  for (int i=0; i<=20; ++i) {
    cosQ = -1.0+(float)i/10.0;

    std::cout << "prodh=" << prodh << "\tcosQ=" << cosQ << std::endl;
    earth.LineProfile(prodh, cosQ, rdet, L, Ne);

    std::list<double>::iterator itrL   (L.begin());
    std::list<double>::iterator itrLend(L.end());
    std::list<double>::iterator itrNe  (Ne.begin());
    for (; itrL!=itrLend; ++itrL, ++itrNe) {
      std::cout << (*itrL) << "\t" << (*itrNe) << std::endl;
    }   
    std::cout << std::endl;
  }
}

//......................................................................

void dumpOscP() 
{
  EarthModel earth("PREM",kTOL);
  PMNS pmns(35.0*M_PI/180.0,
	    45.0*M_PI/180.0,
	    10.0*M_PI/180.0,
	    0.0,
	    8.0E-5,
	    2.43E-3);
  
  TH1F *h1[21], *h2[21], *h3[21], *h4[21];

  double            prodh = 20.0; // production height in km
  double            rdet  = -0.3; // detector height above sea level (km)
  std::list<double> L;
  double            E;
  std::list<double> Ne;
  double            cosQ;
  
  for (int i=2; i<=20; i+=4) {
    E = (float)i;
    char b1[256], b2[256], b3[256], b4[256];
    sprintf(b1,"pmumu_nu_e%d", i);
    sprintf(b2,"pmue_nu_e%d",  i);
    sprintf(b3,"pmumu_anu_e%d",i);
    sprintf(b4,"pmue_anu_e%d", i);
    h1[i] = new TH1F(b1,";cos #theta;P(#mu#mu)",100,-1.0,1.0);
    h2[i] = new TH1F(b2,";cos #theta;P(#mu#e)", 100,-1.0,1.0);
    h3[i] = new TH1F(b3,";cos #theta;P(#mu#mu)",100,-1.0,1.0);
    h4[i] = new TH1F(b4,";cos #theta;P(#mu#e)", 100,-1.0,1.0);
    for (int j=0; j<100; ++j) {
      cosQ = -1.0+2.0*((float)j+0.5)/100.0;
      earth.LineProfile(prodh, cosQ, rdet, L, Ne);
      
      pmns.Reset();
      pmns.PropMatter(L, E, Ne, +1);
      h1[i]->Fill(cosQ, pmns.P(0,1));
      h2[i]->Fill(cosQ, pmns.P(1,1));
      
      pmns.Reset();
      pmns.PropMatter(L, E, Ne, -1);
      h3[i]->Fill(cosQ, pmns.P(0,1));
      h4[i]->Fill(cosQ, pmns.P(1,1));
    }
  }
}

//......................................................................

int main(void) 
{
  TFile f("testEarthModel.root","RECREATE");
  f.cd();

  dumpDensityTable();
  dumpLayers();
  dumpLineProfiles();
  dumpOscP();

  f.Write();
  f.Close();

  std::cout << "Wrote output to testEarthModel.root" << std::endl;
  
  return 0;
}

////////////////////////////////////////////////////////////////////////
