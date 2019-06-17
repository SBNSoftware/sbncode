#include "CAFAna/Analysis/Plots.h"
//#include "CAFAna/Systs/NueSystsSecondAna.h"

#include "CAFAna/Analysis/Style.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SystShifts.h"

#include "Utilities/func/MathUtil.h"

#include "TCanvas.h"
#include "TFeldmanCousins.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <ctime>
#include <functional>
#include <memory>

namespace ana
{
  //----------------------------------------------------------------------
  TH1* DataMCComparison(const Spectrum& data, const Spectrum& mc)
  {
    TH1* ret = 0;

    TH1* hMC = mc.ToTH1(data.POT());
    hMC->SetLineColor(kTotalMCColor);

    TH1* hData = data.ToTH1(data.POT());
    hData->Sumw2();
    hData->SetMarkerStyle(kFullCircle);

    // Need to draw the biggest thing first to get correct axis limits
    if(hMC->GetMaximum() > hData->GetMaximum()+sqrt(hData->GetMaximum())){
      ret = hMC;
      hMC->Draw("hist");
      hData->Draw("ep same");
    }
    else{
      ret = hData;
      hData->Draw("ep");
      hMC->Draw("hist same");
      hData->Draw("ep same"); // Data always goes on top
    }

    gPad->Update();

    return ret;
  }

  //----------------------------------------------------------------------
  TH1* DataMCComparisonAreaNormalized(const Spectrum& data, const Spectrum& mc)
  {
    TH1* ret = 0;

    TH1* hMC = mc.ToTH1(mc.POT());
    hMC->SetLineColor(kTotalMCColor);

    TH1* hData = data.ToTH1(data.POT());
    hData->Sumw2();
    hData->SetMarkerStyle(kFullCircle);

    hMC->Scale(hData->Integral()/hMC->Integral());

    // Need to draw the biggest thing first to get correct axis limits
    if(hMC->GetMaximum() > hData->GetMaximum()+sqrt(hData->GetMaximum())){
      ret = hMC;
      hMC->Draw("hist");
      hData->Draw("ep same");
    }
    else{
      ret = hData;
      hData->Draw("ep");
      hMC->Draw("hist same");
      hData->Draw("ep same"); // Data always goes on top
    }

    gPad->Update();

    return ret;
  }

  //----------------------------------------------------------------------
  TH1* DataMCComparisonComponents(const Spectrum& data,
                                  const IPrediction* mc,
                                  osc::IOscCalculator* calc)
  {
    TH1* ret = 0;

    TH1* hMC = mc->Predict(calc).ToTH1(data.POT());
    hMC->SetLineColor(kTotalMCColor);

    TH1* hNC = mc->PredictComponent(calc,
                                    Flavors::kAll,
                                    Current::kNC,
                                    Sign::kBoth).ToTH1(data.POT());
    hNC->SetLineColor(kNCBackgroundColor);

    TH1* hCC = mc->PredictComponent(calc,
                                    Flavors::kAllNuMu,
                                    Current::kCC,
                                    Sign::kBoth).ToTH1(data.POT());
    hCC->SetLineColor(kNumuBackgroundColor);

    TH1* hBeam = mc->PredictComponent(calc,
                                      Flavors::kNuEToNuE,
                                      Current::kCC,
                                      Sign::kBoth).ToTH1(data.POT());
    hBeam->SetLineColor(kBeamNueBackgroundColor);

    TH1* hData = data.ToTH1(data.POT());
    hData->SetMarkerStyle(kFullCircle);

    // Need to draw the biggest thing first to get correct axis limits
    if(hMC->GetMaximum() > hData->GetMaximum()+sqrt(hData->GetMaximum())){
      ret = hMC;
      hMC->Draw("hist");
      hData->Draw("ep same");
    }
    else{
      ret = hData;
      hData->Draw("ep");
      hMC->Draw("hist same");
      hData->Draw("ep same"); // Data always goes on top
    }

    hNC->Draw("hist same");
    hCC->Draw("hist same");
    hBeam->Draw("hist same");

    gPad->Update();

    return ret;
  }

  TH1* GetMCSystTotal(const IPrediction* mc,
		      osc::IOscCalculator* calc,
		      const SystShifts& shift,
		      std::string hist_name,
		      double pot,
		      bool force1D)
  {
    TH1* hTotal = mc->PredictSyst(calc, shift).ToTHX(pot, force1D);
    hTotal->SetNameTitle((hist_name+"_total").c_str(), (hist_name+"_total").c_str());
    return hTotal;
  }
  
  TH1* GetMCTotal(const IPrediction* mc,
		  osc::IOscCalculator* calc,
		  std::string hist_name,
		  double pot,
		  bool force1D)
  {
    return GetMCSystTotal(mc, calc, kNoShift, hist_name, pot, force1D);
  }

  std::vector<TH1*> GetMCComponents(const IPrediction* mc,
				    osc::IOscCalculator* calc,
				    std::string hist_name,
				    double pot,
				    bool force1D)
  {
    return GetMCSystComponents(mc, calc, kNoShift, hist_name, pot, force1D);
  }

  std::vector<TH1*> GetMCSystComponents(const IPrediction* mc,
					osc::IOscCalculator* calc,
					const SystShifts& shift,
					std::string hist_name,
					double pot,
					bool force1D)
  {
    std::vector<TH1*> ret;
    
    TH1* hTotal = mc->PredictSyst(calc, shift).ToTHX(pot, force1D);
    hTotal->SetNameTitle((hist_name+"_total").c_str(), (hist_name+"_total").c_str());
    ret .push_back(hTotal);
    
    TH1* hNC = mc->PredictComponentSyst(calc,shift,
					Flavors::kAll,
					Current::kNC,
					Sign::kBoth).ToTHX(pot, force1D);
    hNC->SetNameTitle((hist_name+"_NC").c_str(), (hist_name+"_NC").c_str());
    ret .push_back(hNC);
    
    TH1* hAllNumu = mc->PredictComponentSyst(calc,shift,
					     Flavors::kAllNuMu,
					     Current::kCC,
					     Sign::kBoth).ToTHX(pot, force1D);
    hAllNumu->SetNameTitle((hist_name+"_AllNumu").c_str(), (hist_name+"_AllNumu").c_str());
    ret .push_back(hAllNumu);

    TH1* hNumu = mc->PredictComponentSyst(calc,shift,
					  Flavors::kAllNuMu,
					  Current::kCC,
					  Sign::kNu).ToTHX(pot, force1D);
    hNumu->SetNameTitle((hist_name+"_Numu").c_str(), (hist_name+"_Numu").c_str());
    ret .push_back(hNumu);
    
    TH1* hNumubar = mc->PredictComponentSyst(calc,shift,
					     Flavors::kAllNuMu,
					     Current::kCC,
					     Sign::kAntiNu).ToTHX(pot, force1D);
    hNumubar->SetNameTitle((hist_name+"_Numubar").c_str(), (hist_name+"_Numubar").c_str());
    ret .push_back(hNumubar);

    TH1* hAllNutau = mc->PredictComponentSyst(calc,shift,
					      Flavors::kAllNuTau,
					      Current::kCC,
					      Sign::kBoth).ToTHX(pot, force1D);
    hAllNutau->SetNameTitle((hist_name+"_AllNutau").c_str(), (hist_name+"_AllNutau").c_str());
    ret .push_back(hAllNutau);    

    // Want AllSignalNue, SignalNue, SignalNuebar
    TH1* hAllNue = mc->PredictComponentSyst(calc,shift,
					    Flavors::kAllNuE,
					    Current::kCC,
					    Sign::kBoth).ToTHX(pot, force1D);
    hAllNue->SetNameTitle((hist_name+"_AllNue").c_str(), (hist_name+"_AllNue").c_str());
    ret .push_back(hAllNue);    
    
    TH1* hAllBeamNue = mc->PredictComponentSyst(calc,shift,
						Flavors::kNuEToNuE,
						Current::kCC,
						Sign::kBoth).ToTHX(pot, force1D);
    hAllBeamNue->SetNameTitle((hist_name+"_AllBeamNue").c_str(), (hist_name+"_AllBeamNue").c_str());
    ret .push_back(hAllBeamNue);

    TH1* hBeamNue = mc->PredictComponentSyst(calc,shift,
					     Flavors::kNuEToNuE,
					     Current::kCC,
					     Sign::kNu).ToTHX(pot, force1D);
    hBeamNue->SetNameTitle((hist_name+"_BeamNue").c_str(), (hist_name+"_BeamNue").c_str());
    ret .push_back(hBeamNue);

    TH1* hBeamNuebar = mc->PredictComponentSyst(calc,shift,
						Flavors::kNuEToNuE,
						Current::kCC,
						Sign::kAntiNu).ToTHX(pot, force1D);
    hBeamNuebar->SetNameTitle((hist_name+"_BeamNuebar").c_str(), (hist_name+"_BeamNuebar").c_str());
    ret .push_back(hBeamNuebar);
    
    TH1* hAllSignalNue = mc->PredictComponentSyst(calc,shift,
						  Flavors::kNuMuToNuE,
						  Current::kCC,
						  Sign::kBoth).ToTHX(pot, force1D);
    hAllSignalNue->SetNameTitle((hist_name+"_AllSignalNue").c_str(), (hist_name+"_AllSignalNue").c_str());
    ret .push_back(hAllSignalNue);
    
    TH1* hSignalNue = mc->PredictComponentSyst(calc,shift,
					       Flavors::kNuMuToNuE,
					       Current::kCC,
					       Sign::kNu).ToTHX(pot, force1D);
    hSignalNue->SetNameTitle((hist_name+"_SignalNue").c_str(), (hist_name+"_SignalNue").c_str());
    ret .push_back(hSignalNue);
    
    TH1* hSignalNuebar = mc->PredictComponentSyst(calc,shift,
						  Flavors::kNuMuToNuE,
						  Current::kCC,
						  Sign::kAntiNu).ToTHX(pot, force1D);
    hSignalNuebar->SetNameTitle((hist_name+"_SignalNuebar").c_str(), (hist_name+"_SignalNuebar").c_str());
    ret .push_back(hSignalNuebar);
    
    return ret;
  }

  //----------------------------------------------------------------------
  std::vector<TH1*> GetMCTotalForSystShifts(const IPrediction* mc,
					    osc::IOscCalculator* calc,
					    const ISyst* syst, 
					    std::string hist_base_name,
					    double pot,
					    bool force1D)
  {
    std::vector<TH1*> ret;
    std::string syst_name = syst->ShortName();
    for (int i = -3; i < 4; ++i){
      SystShifts s(syst, double(i));
      TH1* hTotal = mc->PredictSyst(calc, s).ToTHX(pot, force1D);
      hTotal->SetNameTitle((hist_base_name+"_total_"+syst_name+"_"+std::to_string(i)).c_str(), 
			   (hist_base_name+"_total_"+syst_name+"_"+std::to_string(i)).c_str());
      ret .push_back(hTotal);
    }
    return ret;
  }

  
  //----------------------------------------------------------------------
  void DataMCRatio(const Spectrum& data,
                   const IPrediction* mc,
                   osc::IOscCalculator* calc,
                   double miny, double maxy)
  {
    DataMCRatio(data, mc->Predict(calc), miny, maxy);
  }

  //----------------------------------------------------------------------
  void DataMCRatio(const Spectrum& data,
                   const Spectrum& mc,
                   double miny, double maxy)
  {
    Ratio ratio(data, mc);
    TH1* h = ratio.ToTH1();
    h->GetYaxis()->SetTitle("Data / MC");
    h->GetYaxis()->SetRangeUser(miny, maxy);
    h->Draw();

    TLine* one = new TLine(h->GetXaxis()->GetXmin(), 1,
                           h->GetXaxis()->GetXmax(), 1);
    one->SetLineWidth(2);
    one->SetLineStyle(7);
    one->Draw();

    gPad->Update();
  }

  //----------------------------------------------------------------------
  void DataMCAreaNormalizedRatio(const Spectrum& data,
                                 const IPrediction* mc,
                                 osc::IOscCalculator* calc,
                                 double miny, double maxy)
  {
    DataMCAreaNormalizedRatio(data, mc->Predict(calc), miny, maxy);
  }

  //----------------------------------------------------------------------
  void DataMCAreaNormalizedRatio(const Spectrum& data,
				 const Spectrum& mc,
                                 double miny, double maxy)
  {
    Spectrum mcScaled = mc;
    mcScaled.OverridePOT(mcScaled.POT() * mc.Integral(1) / data.Integral(1));

    DataMCRatio(data, mcScaled, miny, maxy);
  }

  //----------------------------------------------------------------------
  void RatioPlot(const Spectrum& data,
                 const Spectrum& expected,
                 const Spectrum& fit,
                 double miny, double maxy)
  {
    Ratio fitRatio(fit, expected);
    Ratio dataRatio(data, expected);

    TH1* h = fitRatio.ToTH1();
    h->GetYaxis()->SetTitle("Ratio to expectation");
    h->GetYaxis()->SetRangeUser(miny, maxy);
    h->SetLineColor(kTotalMCColor);
    h->Draw("][");

    h = dataRatio.ToTH1();
    h->SetMarkerStyle(kFullCircle);
    h->Draw("ep same");

    TLine* one = new TLine(h->GetXaxis()->GetXmin(), 1,
                           h->GetXaxis()->GetXmax(), 1);
    one->SetLineWidth(2);
    one->SetLineStyle(7);
    one->Draw();

    gPad->Update();
  }


  //----------------------------------------------------------------------
  void PlotWithSystErrorBand(IPrediction* pred,
                             const std::vector<const ISyst*>& systs,
                             osc::IOscCalculator* calc,
                             double pot,
                             int col, int errCol, float headroom,
			     bool newaxis)
  {

    Spectrum nom = pred->Predict(calc);

    std::vector<Spectrum> ups, dns;

    for(const ISyst* syst: systs){
      SystShifts shifts;
      shifts.SetShift(syst, +1);
      ups.push_back(pred->PredictSyst(calc, shifts));
      shifts.SetShift(syst, -1);
      dns.push_back(pred->PredictSyst(calc, shifts));
    }

    PlotWithSystErrorBand(nom, ups, dns, pot, col, errCol, headroom, newaxis);


  }


  //----------------------------------------------------------------------
  void PlotWithSystErrorBand(const Spectrum& nominal,
                             const std::vector<Spectrum>& upShifts,
                             const std::vector<Spectrum>& downShifts,
                             double pot, int col, int errCol,
                             float headroom, bool newaxis)
  {
    if(col == -1){
      col = kTotalMCColor;
      errCol = kTotalMCErrorBandColor;
    }
    else if(errCol == -1) errCol = col-7; // hopefully a lighter version

    TH1* nom =  nominal.ToTH1(pot);

    std::vector<TH1*> ups, dns;

    for(const auto& upShift:upShifts)     ups.push_back(upShift.ToTH1(pot));
    for(const auto& downShift:downShifts) dns.push_back(downShift.ToTH1(pot));

    nom->SetLineColor(col);
    nom->GetXaxis()->CenterTitle();
    nom->GetYaxis()->CenterTitle();
    if(newaxis) nom->Draw("hist ]["); // Set the axes up

    double yMax = nom->GetBinContent(nom->GetMaximumBin());

    TGraphAsymmErrors* g = new TGraphAsymmErrors;

    for(int binIdx = 0; binIdx < nom->GetNbinsX()+2; ++binIdx){
      const double y = nom->GetBinContent(binIdx);
      g->SetPoint(binIdx, nom->GetXaxis()->GetBinCenter(binIdx), y);

      const double w = nom->GetXaxis()->GetBinWidth(binIdx);

      double errUp = 0;
      double errDn = 0;

      for(unsigned int systIdx = 0; systIdx < ups.size(); ++systIdx){
        double hi = ups[systIdx]->GetBinContent(binIdx)-y;
        double lo = y-dns[systIdx]->GetBinContent(binIdx);

        if(lo <= 0 && hi <= 0) std::swap(lo, hi);

        errUp += hi*hi;
        errDn += lo*lo;

        // TODO: what happens if they're both high or both low?
      } // end for systIdx


      g->SetPointError(binIdx, w/2, w/2, sqrt(errDn), sqrt(errUp));
    } // end for i

    g->SetFillColor(errCol);
    g->Draw("e2 same");
    g->GetYaxis()->SetRangeUser(0, headroom*yMax);
    nom->GetYaxis()->SetRangeUser(0, headroom*yMax);

    nom->Draw("hist ][ same");

    for(TH1* up: ups) delete up;
    for(TH1* dn: dns) delete dn;
  }

  //----------------------------------------------------------------------
  THStack* ToTHStack(const std::vector<std::pair<Spectrum, Color_t>>& s,
                     double pot)
  {
    THStack* ret = new THStack;
    for(auto it: s){
      TH1* h = it.first.ToTH1(pot);
      h->SetFillColor(it.second);
      ret->Add(h, "hist");
    }
    return ret;
  }

  //----------------------------------------------------------------------
  /// Helper for \ref AutoPlaceLegend
  double PointDistanceToBox(double x, double y,
                            double x0, double y0, double x1, double y1)
  {
    // Inside
    if(x > x0 && x < x1 && y > y0 && y < y1) return 0;

    // Corners
    double d = util::sqr(x-x0)+util::sqr(y-y0);
    d = std::min(d, util::sqr(x-x1)+util::sqr(y-y0));
    d = std::min(d, util::sqr(x-x1)+util::sqr(y-y1));
    d = std::min(d, util::sqr(x-x0)+util::sqr(y-y1));

    // Top and bottom edges
    if(x > x0 && x < x1){
      d = std::min(d, util::sqr(y-y0));
      d = std::min(d, util::sqr(y-y1));
    }
    // Left and right
    if(y > y0 && y < y1){
      d = std::min(d, util::sqr(x-x0));
      d = std::min(d, util::sqr(x-x1));
    }

    return d;
  }

  //----------------------------------------------------------------------
  TLegend* AutoPlaceLegend(double dx, double dy, double yPin)
  {
    gPad->Update();

    // Convert requested width and height into physics coordinates
    dx *= (gPad->GetX2()-gPad->GetX1());
    dy *= (gPad->GetY2()-gPad->GetY1());

    // Range of axes in physics units
    const double x0 = gPad->GetUxmin();
    const double x1 = gPad->GetUxmax();
    const double y0 = gPad->GetUymin();
    const double y1 = gPad->GetUymax();

    const double X = x1-x0;
    const double Y = y1-y0;

    double bestd = 0;
    double bestx = 0;
    double besty = 0;
    // If we want to pin Y pos, set it to that now.
    if(yPin >= 0) besty = yPin * (gPad->GetY2()-gPad->GetY1());

    for(bool fallback: {false, true}){
      for(double x = x0+dx/2; x < x1-dx/2; x += X/50){
        for(double y = y0+dy/2; y < y1-dy/2; y += Y/50){
          double d = 999999;

          // Repel from edges
          d = std::min(d, util::sqr((x-dx/2-x0)/X));
          d = std::min(d, util::sqr((x+dx/2-x1)/X));
          d = std::min(d, util::sqr((y-dy/2-y0)/Y));
          d = std::min(d, util::sqr((y+dy/2-y1)/Y));
          if(d < bestd) continue;

          TIter next(gPad->GetListOfPrimitives());
          while(TObject* obj = next()){

            if(!obj->InheritsFrom(TH1::Class())) continue;
            if( obj->InheritsFrom(TH2::Class())) continue;

            TH1* h = (TH1*)obj;

            for(int n = 1; n <= h->GetNbinsX(); ++n){
              const double px = h->GetBinCenter(n);
              const double py = h->GetBinContent(n);

              if(fallback){
                d = std::min(d, util::sqr((px-x)/X)+util::sqr((py-y)/Y));
              }
              else{
                d = std::min(d, PointDistanceToBox(px/X, py/Y,
                                                   (x-dx/2)/X, (y-dy/2)/Y,
                                                   (x+dx/2)/X, (y+dy/2)/Y));
              }
              if(d < bestd) break;
            }
            if(d < bestd) break;
          } // end while

          if(d > bestd){
            bestd = d;
            bestx = x;
            // Update Y if we're not pinning it.
            if (yPin < 0) besty = y;
          }
        } // end for y
      } // end for x

      if(bestd != 0) break; // If we always collide, have to do fallback
    } // end for fallback

    // Convert to pad coordinates
    const double nx = (bestx-gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
    const double ny = (besty-gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());

    const double ndx = dx/(gPad->GetX2()-gPad->GetX1());
    const double ndy = dy/(gPad->GetY2()-gPad->GetY1());

    return new TLegend(nx-ndx/2, ny-ndy/2, nx+ndx/2, ny+ndy/2);
  }

  //----------------------------------------------------------------------
  /*
  void CountingExperimentErrorBarChart(const std::map<std::string, double>& systbars, double statErr, bool bkgdOrSig, bool shortchart)
  {
    std::map<std::string, double> systs = (shortchart) ? SumNueSecondAnaSysts(systbars) : systbars;

    int nbins = systs.size()+1;
    // Systematics are on top, total systematic and (optional) statistical
    // error bar are at the bottom.
    int firstSystBin = 1;
    if(statErr){
      ++nbins;
      ++firstSystBin;
    }

    // Sum all the systematics in qaudrature
    double systTot = 0;
    for(auto it: systs) systTot += it.second*it.second;
    systTot = sqrt(systTot);

    // x-axis range is a multiple of 5% with at least 2% padding
    double xrange = 0;
    while(xrange < std::max(systTot, statErr)+2) xrange += 5;

    TH2* axes = new TH2F("", ";Background uncertainty (%)",
                         10, -xrange, +xrange, nbins, 0, nbins);
    axes->GetXaxis()->CenterTitle();
    if(bkgdOrSig) axes->GetXaxis()->SetTitle("Signal uncertainty (%)");

    // Put the systematics in a vector so we can sort them by size
    std::vector<std::pair<double, std::string>> systList;
    for(auto it: systs) systList.push_back(std::make_pair(it.second, it.first));
    std::sort(systList.begin(), systList.end());

    // Label the y axis
    TAxis* yax = axes->GetYaxis();
    for(unsigned int i = 0; i < systList.size(); ++i)
      yax->SetBinLabel(i+firstSystBin+1, systList[i].second.c_str());

    yax->SetBinLabel(firstSystBin, "Total syst. error");
    if(statErr) yax->SetBinLabel(1, "Statistical error");

    // Make the labels approximately match the x axis title
    yax->SetLabelSize(.06);

    axes->Draw();

    // Blue boxes for each systematic
    for(unsigned int i = 0; i < systList.size(); ++i){
      TBox* box = new TBox(-systList[i].first, i+firstSystBin+.2,
                           +systList[i].first, i+firstSystBin+.8);
      std::cout<<"ERROR = "<<systList[i].first<<std::endl;
      box->SetFillColor(kAzure-8);
      box->Draw();
    }

    // Total systematic error
    TBox* syst = new TBox(-systTot, firstSystBin-.8,
                          +systTot, firstSystBin-.2);
    std::cout<<"TOTAL SYSTEMATIC ERROR = "<<systTot<<std::endl;
    syst->SetFillColor(kAzure-8);
    syst->Draw();

    // Red box for statistical error
    if(statErr){
      TBox* stat = new TBox(-statErr, .2, +statErr, .8);
      stat->SetFillColor(kRed-7);
      stat->Draw();
    }

    // Separate individual systematics from the total and statistical errors
    TLine* div = new TLine(-xrange, firstSystBin, +xrange, firstSystBin);
    div->SetLineWidth(2);
    div->Draw();

    // Vertical line marking the symmetry point
    TLine* zero = new TLine(0, 0, 0, nbins);
    zero->SetLineStyle(7);
    zero->SetLineWidth(2);
    zero->Draw();

    // Leave enough space for the systematic labels
    gPad->SetLeftMargin(.25);
  }
  */
  //----------------------------------------------------------------------
  TGraphAsymmErrors* GraphWithPoissonErrors(const TH1* h, bool noErrorsXaxis, bool drawEmptyBins)
  {
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(h);

    TFeldmanCousins fc(0.6827);//1 sigma
    
    for(int i = 0; i < h->GetNbinsX(); ++i){
      double x, y;
      gr->GetPoint(i, x, y);

      if ( drawEmptyBins || y!=0 ){
	if ( y < 50 )	gr->SetPointEYlow(i, y-fc.CalculateLowerLimit(y,0));
	else            gr->SetPointEYlow(i,sqrt(y));
	if ( y < 30 )   gr->SetPointEYhigh(i,fc.CalculateUpperLimit(y,0)-y);
	else            gr->SetPointEYhigh(i,sqrt(y));
      }

      if(noErrorsXaxis){
        gr->SetPointEXlow(i, 0);
        gr->SetPointEXhigh(i, 0);
      } // Do not use bin width as X-error for points
    }

    return gr;
  }
  //----------------------------------------------------------------------
  TGraph* ShadeBetweenHistograms(TH1* hmin, TH1* hmax)
  {
    int n = hmin->GetNbinsX();
    TGraph* gr = new TGraph(4*n);

    for(int i=1; i<=n; i++)
      {
        double xdown = hmax->GetBinLowEdge(i);
        double xup = hmax->GetBinLowEdge(i+1);
        double ydown = hmin->GetBinContent(i);
        double yup = hmax->GetBinContent(i);

        gr->SetPoint(2*(i-1), xdown, ydown);
        gr->SetPoint(2*(i-1)+1, xup, ydown);
        gr->SetPoint(4*n-2*i, xup, yup);
        gr->SetPoint(4*n-2*i+1, xdown, yup);
      }

      gr->SetFillColor(hmin->GetLineColor() - 9); // In principle this is lighter
      gr->SetLineColor(hmin->GetLineColor() - 9); // In case one does Draw("lf")

      return gr;
  }
  //----------------------------------------------------------------------
  TGraphAsymmErrors * ProfileQuantile(const TH2 * hist,
                                      const std::string & axisName,
                                      const std::string & graphName,
                                      const std::pair<double, double> & quantileDivisions)
  {
    if (hist->GetDimension() != 2)
            throw std::runtime_error(Form("Can't profile a histogram with other than 2 dimensions.  Yours is a %d-D histogram...", hist->GetDimension()));

    const TAxis * axis = nullptr;
    std::function<TH1D*(const char *, Int_t, Int_t, Option_t*)> projectionMethod;
    if (axisName == "x" || axisName == "X")
    {
      axis = hist->GetXaxis();
      projectionMethod = std::bind(&TH2::ProjectionY, hist, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    else if (axisName == "y" || axisName == "Y")
    {
      axis = hist->GetYaxis();
      projectionMethod = std::bind(&TH2::ProjectionX, hist, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    else
      throw std::runtime_error( Form("Can't profile onto unknown axis: '%s'", axisName.c_str()) );


    std::vector<double> points_x, points_y, errors_x_low, errors_x_high, errors_y_low, errors_y_high;

    double quantiles[2] = {0, 0};
    double _quantileDivisions[2] = {quantileDivisions.first, quantileDivisions.second};
    for (int binNum = 1; binNum <= axis->GetNbins(); binNum++)  // remember, 0 is underflow and Nbins+1 is overflow...
    {
      std::unique_ptr<TH1D> proj ( projectionMethod(Form("%s_tmp_%d", hist->GetName(), int(std::time(nullptr))), binNum, binNum, "") );  // randomish name, just to be safe

      double y = proj->GetMean();
      points_x.push_back(axis->GetBinCenter(binNum));
      points_y.push_back(y);

      // for now, make the errors for an empty projection 0
      unsigned int n_qs = 0;
      if (proj->Integral() == 0)
      {
        n_qs = 2;
        quantiles[0] = quantiles[1] = 0;
      }
      else
        n_qs = proj->GetQuantiles(2, quantiles, _quantileDivisions);
      if (n_qs != 2)
        throw std::runtime_error( Form("GetQuantiles() didn't compute all the quantiles in HistoTools.ProfileQuantile().  I requested 2 quantiles, but got %d of them...", n_qs) );

      double binWidth = axis->GetBinWidth(binNum);
      errors_x_low.push_back(binWidth/2.);
      errors_x_high.push_back(binWidth/2.);
      errors_y_low.push_back(y-quantiles[0]);
      errors_y_high.push_back(quantiles[1]-y);
    }

    TGraphAsymmErrors * outGraph = new TGraphAsymmErrors(
            points_x.size(),
            &points_x[0],
            &points_y[0],
            &errors_x_low[0],
            &errors_x_high[0],
            &errors_y_low[0],
            &errors_y_high[0]
    );
    std::string name = (graphName.size()) ? graphName : Form("%s_quantile_errors", hist->GetName());
    outGraph->SetName(name.c_str());

    return outGraph;

  }
} // namespace
