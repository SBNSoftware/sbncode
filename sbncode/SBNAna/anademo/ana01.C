// Make a plot with cuts 
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"

using namespace ana;

#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

void ana01()
{
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  const std::string fname = "larout.caf.root";
  // Source of events
  SpectrumLoader loader(fname);

  // ---- VARS -----
  // A Var returns a number per slice, a.k.a. variables to plot
  const Var kNTrk       ([](const caf::SRProxy* sr)
                         {
                           return sr->reco.ntrk;
                         });

  const Var kTrkLen   ([](const caf::SRProxy* sr)
		       { //length of 1st track
			 double len = -5.0;
			 if (sr->reco.ntrk > 0)
			   len = sr->reco.trk[0].len;
			 return len;
		       });

  // ---- CUTS -----
  // A Cut returns a boolean per slice 
  const Cut kLong      ([](const caf::SRProxy* sr)
			{
			  bool pass = false;
			  if (sr->reco.ntrk > 0) //check that info is filled!
			    pass = sr->reco.trk[0].len > 50;
			  return pass;
			});

  // ---- SPECTRA -----
  // A spectrum is a histogram with associated POT information
  const Binning binsN   = Binning::Simple(10, 0, 10);
  const Binning binsLen = Binning::Simple(50, 0, 500);

  //               axis( <Title>,Binning,Var) 
  const HistAxis axNTrk("Number of Tracks", binsN, kNTrk);
  const HistAxis axLen ("1st Track Length", binsLen, kTrkLen);

  //         spectrum(Spectrumloader,HistAxis,Cut) 
  Spectrum sNTracks  (loader, axNTrk, kNoCut);
  Spectrum sTrkLenAll(loader, axLen, kNoCut);
  Spectrum sTrkLen50 (loader, axLen, kLong);

  // This is the call that actually fills in the spectrum
  loader.Go();


  // ---- DRAW -----
  // Suppose you want to scale to 6e20
  const double pot = 6e20;

  // For plotting purposes we can convert spectra to a TH1
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1* hNTracks = sNTracks.ToTH1(pot);
  hNTracks->Draw("hist");

  TCanvas *c2 = new TCanvas("c2","c2");
  TH1* hLenAll = sTrkLenAll.ToTH1(pot);
  TH1* hLen50  =  sTrkLen50.ToTH1(pot);
  hLenAll->Draw("hist");
  hLen50->SetLineColor(kPink+2);
  hLen50->Draw("hist same");
  //Drawing a legend because we are not hobos
  TLegend *leg = new TLegend(0.6,0.7,0.85,0.8);
  leg->AddEntry(hLenAll,"All","l");
  leg->AddEntry(hLen50,"Pass length cut","l");
  leg->Draw();

}
