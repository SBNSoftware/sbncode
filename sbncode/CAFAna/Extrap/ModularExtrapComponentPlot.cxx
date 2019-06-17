#if 0
// Currently dead code

#include "CAFAna/Extrap/ModularExtrapComponent.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Utilities.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

namespace ana
{

  //---------------------------------------------------------------------------

  void ModularExtrapComponent::ComparisonPlot(
      Spectrum mc,
      Spectrum notMC,
      double pot,
      std::string notMCLabel,
      int notMCColor,
      std::string latex,
      std::string title,
      std::string saveAs,
      bool restrictRange
  ){
    TH1* histMC( mc.ToTH1(pot) );
    TH1* histNotMC( notMC.ToTH1(pot) );
    histMC->SetLineColor( kRed );
    histNotMC->SetLineColor( notMCColor );
    TCanvas c;
    histMC->SetMaximum( 1.5 * histMC->GetMaximum() );
    if (restrictRange) histMC->GetXaxis()->SetRangeUser( 0.5, 3.5);
    histMC->SetTitle( (latex+" "+title).c_str() );
    histMC->Draw("hist");
    histNotMC->Draw("same");
    TLegend legend( .15, .7, .4, .85 );
    legend.SetTextSize(0.05);
    legend.AddEntry( histNotMC, notMCLabel.c_str() );
    legend.AddEntry( histMC, "MC" );
    legend.SetFillStyle(0);
    legend.Draw();
    c.Write( saveAs.c_str() );
    delete histMC;
    delete histNotMC;
  }

  //---------------------------------------------------------------------------

  void NoReweight::SavePlots(TDirectory* dir, double potFD) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    TH1* hist( fRecoFD.Unoscillated().ToTH1(potFD) );
    hist->SetLineColor( kRed );
    TCanvas c;
    hist->SetMaximum( 1.5 * hist->GetMaximum() );
    hist->SetTitle( "FD Reco Spectrum" );
    hist->Draw("hist");
    c.Write("MC");
    delete hist;

    tmp->cd();
  }

  //---------------------------------------------------------------------------

  void TruthReweight::SavePlots(TDirectory* dir, double potFD) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    Spectrum decompresult( GetDecompResult(fDecomp,fDecompRes) );
    double potND( decompresult.POT() );

    ComparisonPlot(
      fRecoToTrueND.Unoscillated(), decompresult, potND,
      "Data", kBlack, fLatex, "ND Reco Spectrum", "NDReco" );

    TH2* histNDrtt( fRecoToTrueND.ToTH2(potND) );
    TCanvas cNDrtt;
    histNDrtt->GetYaxis()->SetRangeUser( 0.5, 3.5 );
    histNDrtt->SetTitle( (fLatex+" ND Reco to True Spectrum (MC)").c_str() );
    histNDrtt->Draw("colz");
    cNDrtt.Write("NDRecoToTrue");
    delete histNDrtt;

    OscillatableSpectrum recoToTrueND( fRecoToTrueND );
    recoToTrueND.ReweightToRecoSpectrum( decompresult );
    ComparisonPlot(
      fRecoToTrueND.TrueEnergy(), recoToTrueND.TrueEnergy(), potND,
      "Reweighted", kBlue, fLatex, "ND True Energy Spectrum", "NDTrue",
      true );

    ComparisonPlot(
      fTrueToRecoFD.TrueEnergy(), Return().TrueEnergy(), potFD,
      "Extrapolated", kBlue, fLatex, "FD True Energy Spectrum", "FDTrue",
      true );

    TH2* histFDttr( fTrueToRecoFD.ToTH2(potFD) );
    TCanvas cFDttr;
    histFDttr->GetYaxis()->SetRangeUser( 0.5, 3.5 );
    histFDttr->SetTitle( (fLatex+" FD True To Reco Spectrum (MC)").c_str() );
    histFDttr->Draw("colz");
    cFDttr.Write("FDTrueToReco");
    delete histFDttr;

    ComparisonPlot(
      fTrueToRecoFD.Unoscillated(), Return().Unoscillated(), potFD,
      "Extrapolated", kBlue, fLatex, "FD Reco Spectrum", "FDReco" );

    TH1* histNDDatReco = decompresult.ToTH1(potND);
    TH1* histNDRwtTrue = recoToTrueND.TrueEnergy().ToTH1(potND);
    TH1* histFDExtTrue = Return().TrueEnergy().ToTH1(potFD);
    TH1* histFDExtReco = Return().Unoscillated().ToTH1(potFD);

    TH1* histNDMCReco = fRecoToTrueND.Unoscillated().ToTH1(potND);
    TH1* histNDMCTrue = fRecoToTrueND.TrueEnergy().ToTH1(potND);
    TH1* histFDMCTrue = fTrueToRecoFD.TrueEnergy().ToTH1(potFD);
    TH1* histFDMCReco = fTrueToRecoFD.Unoscillated().ToTH1(potFD);

    histNDDatReco->Write("NDDatReco");
    histNDRwtTrue->Write("NDRwtTrue");
    histFDExtTrue->Write("FDExtTrue");
    histFDExtReco->Write("FDExtReco");

    histNDMCReco->Write("NDMCReco");
    histNDMCTrue->Write("NDMCTrue");
    histFDMCTrue->Write("FDMCTrue");
    histFDMCReco->Write("FDMCReco");

    delete histNDDatReco;
    delete histNDRwtTrue;
    delete histFDExtTrue;
    delete histFDExtReco;

    delete histNDMCReco;
    delete histNDMCTrue;
    delete histFDMCTrue;
    delete histFDMCReco;

    tmp->cd();
  }

  //---------------------------------------------------------------------------

  void RecoReweight::SavePlots(TDirectory* dir, double potFD) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    Spectrum fdMCReco( fTrueToRecoFD.Unoscillated() );
    Spectrum fdExtReco( Return().Unoscillated() );

    Spectrum decompresult( GetDecompResult(*fDecomp,fDecompRes) );
    double potND( decompresult.POT() );

    ComparisonPlot(
      fRecoND, decompresult, potND,
      "Data", kBlack, fLatex, "ND Reco Spectrum", "ND" );

    ComparisonPlot(
      fdMCReco, fdExtReco, potFD,
      "Extrapolated", kBlue, fLatex, "FD Reco Spectrum", "FD" );

    TH1* histDvMCND( (decompresult/fRecoND).ToTH1() );
    TH1* histEvMCFD( (fdExtReco / fdMCReco).ToTH1() );
    histDvMCND->SetLineColor( kRed );
    histEvMCFD->SetLineColor( kBlue );
    histEvMCFD->SetLineWidth( 4 );
    histEvMCFD->SetLineStyle( 2 );
    TCanvas cDvMC;
    histDvMCND->SetMaximum( 1.5 * histDvMCND->GetMaximum() );
    histDvMCND->SetTitle( (fLatex+" Data/MC Ratio").c_str() );
    histDvMCND->GetYaxis()->SetTitle("Ratio");
    histDvMCND->Draw("hist");
    histEvMCFD->Draw("same""hist");
    TLegend legendDvMC( .15, .7, .4, .85 );
    legendDvMC.SetTextSize(0.05);
    legendDvMC.AddEntry( histEvMCFD, "FD (Extrapolated)" );
    legendDvMC.AddEntry( histDvMCND, "ND" );
    legendDvMC.SetFillStyle(0);
    legendDvMC.Draw();
    cDvMC.Write("DataOverMC");
    delete histDvMCND;
    delete histEvMCFD;

    TH1* histFoNMC( (fdMCReco/fRecoND).ToTH1() );
    TH1* histFoNDE( (fdExtReco/decompresult).ToTH1() );
    histFoNMC->SetLineColor( kRed );
    histFoNDE->SetLineColor( kBlue );
    histFoNDE->SetLineWidth( 4 );
    histFoNDE->SetLineStyle( 2 );
    TCanvas cFoN;
    histFoNMC->SetMaximum( 1.5 * histFoNMC->GetMaximum() );
    histFoNMC->SetTitle( (fLatex+" Far/Near Ratio").c_str() );
    histFoNMC->GetYaxis()->SetTitle("Ratio");
    histFoNMC->Draw("hist");
    histFoNDE->Draw("same""hist");
    TLegend legendFoN( .15, .7, .4, .85 );
    legendFoN.SetTextSize(0.05);
    legendFoN.AddEntry( histFoNDE, "Extrapolation/Data" );
    legendFoN.AddEntry( histFoNMC, "MC" );
    legendFoN.SetFillStyle(0);
    legendFoN.Draw();
    cFoN.Write("FarOverNear");
    delete histFoNMC;
    delete histFoNDE;

    TH1* histNDDat = decompresult.ToTH1(potND);
    TH1* histFDExt = fdExtReco.ToTH1(potFD);
    TH1* histNDMC = fRecoND.ToTH1(potND);
    TH1* histFDMC = fdMCReco.ToTH1(potFD);

    histNDDat->Write("NDDat");
    histFDExt->Write("FDExt");
    histNDMC->Write("NDMC");
    histFDMC->Write("FDMC");

    delete histNDDat;
    delete histFDExt;
    delete histNDMC;
    delete histFDMC;

    tmp->cd();
  }

  //---------------------------------------------------------------------------

  TH1* RecoReweight::OptimalBinning(
    double potMCFD,
    double potMCND
  ) const {

    TH1* fd = fTrueToRecoFD.Unoscillated().ToTH1(potMCFD);
    TH1* nd = fRecoND.ToTH1(potMCND);

    TH1* ob = OptimalBinningHelper( fd, nd );

    delete fd;
    delete nd;

    return ob;

  }

  //---------------------------------------------------------------------------

  TH1* RecoReweight::OptimalBinningFit(
    double potMCFD,
    double potMCND
  ) const {

    TH1* fd = fTrueToRecoFD.Unoscillated().ToTH1(potMCFD);
    TH1* nd = fRecoND.ToTH1(potMCND);

    const double lower( fd->GetXaxis()->GetXmin() );
    const double upper( fd->GetXaxis()->GetXmax() );

    TF1* fdGauss = new TF1( UniqueName().c_str() ,"gaus(0)", lower, upper );
    TF1* ndGauss = new TF1( UniqueName().c_str() ,"gaus(0)", lower, upper );

    fdGauss->SetParameters( 100., 2., .5 );
    ndGauss->SetParameters( 100., 2., .5 );

    fd->Fit( fdGauss, "WW" );
    nd->Fit( ndGauss, "WW" );

    TH1* fdFine = new TH1D(UniqueName().c_str(),"fdFine",100,lower,upper);
    TH1* ndFine = new TH1D(UniqueName().c_str(),"ndFine",100,lower,upper);

    fdFine->Eval(fdGauss);
    ndFine->Eval(ndGauss);

    const double binScale = fdFine->GetBinWidth(1) / fd->GetBinWidth(1);

    fdFine->Scale(binScale);
    ndFine->Scale(binScale);

    TH1* ob(OptimalBinningHelper(fdFine,ndFine));

    delete fd;
    delete nd;

    return ob;

  }

  //---------------------------------------------------------------------------

  TH1* RecoReweight::OptimalBinningHelper( TH1* fd, TH1* nd)
  {

    TH1* ob( (TH1*)fd->Clone() );

    int nbins = fd->GetNbinsX();
    assert( (nbins == nd->GetNbinsX()) && "Bin Mismatch" );

    for (int bin(1); bin<=nbins; ++bin)
    {
      double invdensFD( fd->GetBinWidth(bin) / fd->GetBinContent(bin) );
      double invdensND( nd->GetBinWidth(bin) / nd->GetBinContent(bin) );
      double ratio( fd->GetBinContent(bin) / nd->GetBinContent(bin) );
      double ratioPrev( fd->GetBinContent(bin-1) / nd->GetBinContent(bin-1) );
      double ratioNext( fd->GetBinContent(bin+1) / nd->GetBinContent(bin+1) );
      double deltaRatio( ratioNext - ratioPrev );
      double deltaVar( fd->GetBinCenter(bin+1) - fd->GetBinCenter(bin-1) );
      double invslope( deltaVar / deltaRatio );
      double widthCubed( 6. * ( ratio*ratio * invslope*invslope )
                            * ( invdensFD + invdensND ) );
      double optbin( pow(widthCubed, 1./3.) );
      if ( std::isfinite(optbin) ) ob->SetBinContent( bin, optbin );
      else ob->SetBinContent( bin, 0. );
    }
    ob->SetBinContent( 0, 0. );
    ob->SetBinContent( nbins+1, 0. );

    return ob;

  }

}

#endif
