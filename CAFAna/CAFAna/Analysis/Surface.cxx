#include "CAFAna/Analysis/Surface.h"

#include "CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/IFitVar.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/ThreadPool.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/func/IOscCalculator.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TMarker.h"
#include "Minuit2/StackAllocator.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TCollection.h"

#include <iostream>
#include <functional>

namespace ana
{
  //----------------------------------------------------------------------
  Surface::Surface(const IExperiment* expt,
                   osc::IOscCalculatorAdjustable* calc,
                   const IFitVar* xvar, int nbinsx, double xmin, double xmax,
                   const IFitVar* yvar, int nbinsy, double ymin, double ymax,
                   const std::vector<const IFitVar*>& profVars,
                   const std::vector<const ISyst*>& profSysts,
                   const std::map<const IFitVar*, std::vector<double>>& seedPts,
                   const std::vector<SystShifts>& systSeedPts,
                   bool parallel,
                   Fitter::Precision prec)
    : fParallel(parallel), fPrec(prec)
  {
    fHist = ExpandedHistogram(";"+xvar->LatexName()+";"+yvar->LatexName(),
                              nbinsx, xmin, xmax,
                              nbinsy, ymin, ymax);

    for(unsigned int i = 0; i < profVars.size()+profSysts.size(); ++i){
      std::string title;
      if(i < profVars.size())
        title = profVars[i]->LatexName();
      else
        title = profSysts[i-profVars.size()]->LatexName();

      fProfHists.push_back(ExpandedHistogram(title+";"+xvar->LatexName()+";"+yvar->LatexName(),
                                             nbinsx, xmin, xmax,
                                             nbinsy, ymin, ymax));
    }

    std::string progTitle = "Filling surface for "+yvar->ShortName()+" vs "+xvar->ShortName();

    if(!profVars.empty() || !profSysts.empty()){
      progTitle += " (profiling ";

      for(const IFitVar* v: profVars) {
        progTitle += v->ShortName() + ", ";
	fSeedValues.push_back( v->GetValue( calc ) );
      }
      for(const ISyst* s: profSysts)
        progTitle += s->ShortName() + ", ";

      // Always have one superfluous ", " at the end
      progTitle.resize(progTitle.size()-2);
      progTitle += ")";
    }

    FillSurface(progTitle, expt, calc, xvar, yvar, profVars, profSysts, seedPts, systSeedPts);

    // Location of the best minimum found from filled surface
    double minchi = 1e10;
    int minx = nbinsx/2;
    int miny = nbinsy/2;
    for(int x = 1; x < nbinsx+1; ++x){
      for(int y = 1; y < nbinsy+1; ++y){
        const double chi = fHist->GetBinContent(x, y);
        if(chi < minchi){
          minchi = chi;
          minx = x;
          miny = y;
        }
      }
    }

    std::vector<const IFitVar*> allVars = {xvar, yvar};
    allVars.insert(allVars.end(), profVars.begin(), profVars.end());
    Fitter fit(expt, allVars, profSysts);
    fit.SetPrecision(fPrec);
    // Seed from best grid point
    xvar->SetValue(calc, fHist->GetXaxis()->GetBinCenter(minx));
    yvar->SetValue(calc, fHist->GetYaxis()->GetBinCenter(miny));
    for(int i = 0; i < (int)fSeedValues.size(); ++i) profVars[i]->SetValue( calc, fSeedValues[i] );
    SystShifts systSeed = SystShifts::Nominal();
    fMinChi = fit.Fit(calc, systSeed, seedPts);
    fMinX = xvar->GetValue(calc);
    fMinY = yvar->GetValue(calc);

    for(int x = 0; x < nbinsx+2; ++x){
      for(int y = 0; y < nbinsy+2; ++y){
        fHist->SetBinContent(x, y, fHist->GetBinContent(x, y)-fMinChi);
      }
    }

    fHist->SetMinimum(0);
  }
  //---------------------------------------------------------------------
  void Surface::FillSurface(const std::string& progTitle,
                            const IExperiment* expt,
                            osc::IOscCalculatorAdjustable* calc,
                            const IFitVar* xvar, const IFitVar* yvar,
                            const std::vector<const IFitVar*>& profVars,
                            const std::vector<const ISyst*>& profSysts,
                            const std::map<const IFitVar*, std::vector<double>>& seedPts,
                            const std::vector<SystShifts>& systSeedPts)
  {
    if(fParallel && !(profVars.empty() && profSysts.empty())){
      // Minuit calls this at some point, and it creates a static. Which
      // doesn't like happening on multiple threads at once. So do it upfront
      // while we're still single-threaded.
      ROOT::Minuit2::StackAllocatorHolder::Get();
    }

    // Nothing created during surface filling belongs in a
    // directory. Unfortunately the local guards in Spectrum etc are racey when
    // run in parallel. But this should cover the whole lot safely.
    DontAddDirectory guard;

    Progress* prog = 0;
    // Difficult to keep a progress bar properly up to date when threaded
    if(!fParallel) prog = new Progress(progTitle);
    ThreadPool* pool = 0;
    if(fParallel){
      pool = new ThreadPool;
      pool->ShowProgress(progTitle);
    }

    const int Nx = fHist->GetNbinsX();
    const int Ny = fHist->GetNbinsY();

    // Fill bins in "random" order so that the progress bar is accurate
    // straight away instead of possibly being misled by whatever atypical
    // points we start with. This step is a prime which guarantees we get every
    // cell.
    int step = 7919;
    // Very unlikely (Nx or Ny is a multiple of step), but just to be safe.
    if((Nx*Ny)%step == 0) step = 1;

    int bin = 0;
    int neval = 0;

    do{
      const int x = bin%Nx+1;
      const int y = bin/Nx+1;

      const double xv = fHist->GetXaxis()->GetBinCenter(x);
      const double yv = fHist->GetYaxis()->GetBinCenter(y);

      // For parallel running need to set these at point of use otherwise we
      // race.
      if(!fParallel){
        xvar->SetValue(calc, xv);
        yvar->SetValue(calc, yv);
      }

      if(xvar->Penalty(xv, calc) > 1e-10){
	std::cerr << "Warning! " << xvar->ShortName() << " = " << xv
		  << " has penalty of " << xvar->Penalty(xv, calc)
		  << " that could have been applied in surface. "
		  << "This should never happen." << std::endl;
      }
      if(yvar->Penalty(yv, calc) > 1e-10){
        std::cerr << "Warning! " << yvar->ShortName() << " = " << yv
                  << " has penalty of " << yvar->Penalty(yv, calc)
                  << " that could have been applied in surface. "
                  << "This should never happen." << std::endl;
      }

      if(fParallel){
        pool->AddMemberTask(this, &Surface::FillSurfacePoint,
                            expt, calc,
                            xvar, xv, yvar, yv,
                            profVars, profSysts, seedPts, systSeedPts);
      }
      else{
        FillSurfacePoint(expt, calc,
                         xvar, xv, yvar, yv,
                         profVars, profSysts, seedPts, systSeedPts);
        ++neval;
        prog->SetProgress(neval/double(Nx*Ny));
      }

      bin = (bin+step)%(Nx*Ny);
    } while(bin != 0);


    if(fParallel){
      pool->Finish();
      delete pool;
    }
    else{
      prog->Done();
      delete prog;
    }
  }

  /// \brief Helper for Surface::FillSurfacePoint
  ///
  /// The cacheing of the nominal done in PredictionInterp is not
  /// threadsafe. This is an inelegant but pragmatic way of surpressing it.
  class OscCalcNoHash: public osc::IOscCalculatorAdjustable
  {
  public:
    OscCalcNoHash(osc::IOscCalculatorAdjustable* c) : fCalc(c) {}

    osc::IOscCalculatorAdjustable* Copy() const override
    {
      std::cout << "Surface::OscCalcNoHash not copyable." << std::endl;
      abort();
    }
    double P(int a, int b, double E) override {return fCalc->P(a, b, E);}
    /// Marks this calculator "unhashable" so the cacheing won't occur
    TMD5* GetParamsHash() const override {return 0;}

    // It's a shame we have to forward all the getters and setters explicitly,
    // but I don't see how to avoid it. Take away some drudgery with a macro.
#define F(v)\
    void Set##v(double x) override {fCalc->Set##v(x);}\
    double Get##v() const override {return fCalc->Get##v();}
    F(L); F(Rho); F(Dmsq21); F(Dmsq32); F(Th12); F(Th13); F(Th23); F(dCP);
#undef F

  protected:
    osc::IOscCalculatorAdjustable* fCalc;
  };

  //----------------------------------------------------------------------
  void Surface::FillSurfacePoint(const IExperiment* expt,
                                 osc::IOscCalculatorAdjustable* calc,
                                 const IFitVar* xvar, double x,
                                 const IFitVar* yvar, double y,
                                 const std::vector<const IFitVar*>& profVars,
                                 const std::vector<const ISyst*>& profSysts,
                                 const std::map<const IFitVar*, std::vector<double>>& seedPts,
                                 const std::vector<SystShifts>& systSeedPts)
  {
    osc::IOscCalculatorAdjustable* calcNoHash = 0; // specific to parallel mode

    if(fParallel){
      // Need to take our own copy so that we don't get overwritten by someone
      // else's changes.
      calc = calc->Copy();
      xvar->SetValue(calc, x);
      yvar->SetValue(calc, y);

      calcNoHash = new OscCalcNoHash(calc);
    }

    //Make sure that the profiled values of fitvars do not persist between steps.
    for(int i = 0; i < (int)fSeedValues.size(); ++i) profVars[i]->SetValue( calc, fSeedValues[i] );

    double chi;
    if(profVars.empty() && profSysts.empty()){
      chi = expt->ChiSq(fParallel ? calcNoHash : calc);
    }
    else{
      Fitter fitter(expt, profVars, profSysts);
      fitter.SetPrecision(fPrec);
      SystShifts bestSysts;
      chi = fitter.Fit(calc, bestSysts, seedPts, systSeedPts, Fitter::kQuiet);

      for(unsigned int i = 0; i < profVars.size(); ++i){
        fProfHists[i]->Fill(x, y, profVars[i]->GetValue(calc));
      }
      for(unsigned int j = 0; j < profSysts.size(); ++j){
        fProfHists[j+profVars.size()]->Fill(x, y, bestSysts.GetShift(profSysts[j]));
      }
    }

    fHist->Fill(x, y, chi);

    if(fParallel){
      delete calc;
      delete calcNoHash;
    }
  }

  //---------------------------------------------------------------------
  void Surface::EnsureAxes() const
  {
    // Could have a file temporarily open
    DontAddDirectory guard;

    // If this pad has already been drawn in, already has axes
    if(gPad && !gPad->GetListOfPrimitives()->IsEmpty()) return;

    // Old, hackier solution
    /*
    std::cout << gPad->GetListOfPrimitives()->GetEntries() << std::endl;
    // Which pads have we already drawn axes in? Never draw axes in them
    // again. Unfortunately UniqueID() never seems to be set. If that's the
    // case, set it to a random value and hope...
    static std::set<UInt_t> already;
    if(already.count(gPad->GetUniqueID())) return;
    if(gPad->GetUniqueID() == 0) gPad->SetUniqueID(rand());
    already.insert(gPad->GetUniqueID());
    */

    const TAxis* ax = fHist->GetXaxis();
    const TAxis* ay = fHist->GetYaxis();
    const double Nx = ax->GetNbins();
    const double Ny = ay->GetNbins();

    // Axes with limits where the user originally requested them, which we
    // adjusted to be the centres of the first and last bins.
    TH2* axes = new TH2C(UniqueName().c_str(),
                         TString::Format(";%s;%s",
                                         ax->GetTitle(), ay->GetTitle()),
                         Nx-1, ax->GetBinCenter(1), ax->GetBinCenter(Nx),
                         Ny-1, ay->GetBinCenter(1), ay->GetBinCenter(Ny));
    axes->Draw();

    if(fHist){
      // "colz same" will reuse axis's min and max, so set them helpfully here
      axes->SetMinimum(fHist->GetMinimum());
      axes->SetMaximum(fHist->GetMaximum());
    }

    axes->SetTitle(fHist->GetTitle());
    axes->GetXaxis()->SetLabelSize(ax->GetLabelSize());
    axes->GetYaxis()->SetLabelSize(ay->GetLabelSize());
    axes->GetXaxis()->SetLabelOffset(ax->GetLabelOffset());
    axes->GetYaxis()->SetLabelOffset(ay->GetLabelOffset());
    axes->GetXaxis()->CenterTitle();
    axes->GetYaxis()->CenterTitle();
    gPad->Update();
  }

  //---------------------------------------------------------------------
  void Surface::Draw() const
  {
    EnsureAxes();

    fHist->Draw("colz same");

    // colz obliterated them
    gPad->RedrawAxis();

    gPad->Update();
  }

  //---------------------------------------------------------------------
  void Surface::DrawBestFit(Color_t color, Int_t marker) const
  {
    EnsureAxes();

    TMarker* mark = new TMarker(fMinX, fMinY, marker);
    mark->SetMarkerSize(1.5);
    mark->SetMarkerColor(color);
    mark->Draw();
    gPad->Update();
  }

  //----------------------------------------------------------------------
  std::vector<TGraph*> Surface::GetGraphs(TH2* fc, double minchi)
  {
    std::vector<TGraph*> ret;

    if(minchi < 0) minchi = fMinChi;
    std::unique_ptr<TH2> surf((TH2*)fHist->Clone(UniqueName().c_str()));
    surf->Add(fc, -1);

    TVirtualPad* bak = gPad;

    const bool wasbatch = gROOT->IsBatch();
    gROOT->SetBatch(); // User doesn't want to see our temporary canvas
    TCanvas tmp;

    gStyle->SetOptStat(0);

    const double level = minchi-fMinChi;
    surf->SetContour(1, &level);
    surf->Draw("cont list");

    tmp.Update();
    tmp.Paint();

    gROOT->SetBatch(wasbatch);
    gPad = bak;

    // The graphs we need (contained inside TLists, contained inside
    // TObjArrays) are in the list of specials. But we need to be careful about
    // types, because other stuff can get in here too (TDatabasePDG for
    // example).
    TCollection* specs = gROOT->GetListOfSpecials();

    TIter nextSpec(specs);
    while(TObject* spec = nextSpec()){
      if(!spec->InheritsFrom(TObjArray::Class())) continue;
      TObjArray* conts = (TObjArray*)spec;

      if(conts->IsEmpty()) continue;

      if(!conts->At(0)->InheritsFrom(TList::Class())) continue;
      TList* cont = (TList*)conts->At(0);

      TIter nextObj(cont);
      // Contour could be split into multiple pieces
      while(TObject* obj = nextObj()){
        if(!obj->InheritsFrom(TGraph::Class())) continue;

        ret.push_back((TGraph*)obj->Clone(UniqueName().c_str()));
      } // end for obj
    } // end for spec

    return ret;
  }

  //----------------------------------------------------------------------
  void Surface::DrawContour(TH2* fc, Style_t style, Color_t color,
                            double minchi)
  {
    EnsureAxes();

    std::vector<TGraph*> gs = GetGraphs(fc, minchi);

    for(TGraph* g: gs){
      g->SetLineWidth(3);//2);
      g->SetLineStyle(style);
      g->SetLineColor(color);
      g->Draw("l");
    }

    gPad->Update();
  }

  //----------------------------------------------------------------------
  TH2* Surface::ToTH2(double minchi) const
  {
    // Could have a file temporarily open
    DontAddDirectory guard;

    TH2* ret = new TH2F(*fHist);

    if(minchi >= 0){
      for(int x = 0; x < ret->GetNbinsX()+2; ++x){
        for(int y = 0; y < ret->GetNbinsY()+2; ++y){
          ret->SetBinContent(x, y, ret->GetBinContent(x, y)+fMinChi-minchi);
        }
      }
    }

    return ret;
  }

  //----------------------------------------------------------------------
  void Surface::SetTitle(const char* str)
  {
    fHist->SetTitle(str);
  }

  //----------------------------------------------------------------------
  /// Helper function for the gaussian approximation surfaces
  TH2* Flat(double level, const Surface& s)
  {
    TH2* h = s.ToTH2();

    for(int x = 0; x < h->GetNbinsX()+2; ++x)
      for(int y = 0; y < h->GetNbinsY()+2; ++y)
        h->SetBinContent(x, y, level);

    return h;
  }

  // See eg the statistics section of the PDG
  TH2* Gaussian68Percent2D(const Surface& s){return Flat(2.30, s);}
  TH2* Gaussian90Percent2D(const Surface& s){return Flat(4.61, s);}
  TH2* Gaussian95Percent2D(const Surface& s){return Flat(5.99, s);}
  TH2* Gaussian2Sigma2D   (const Surface& s){return Flat(6.18, s);}
  TH2* Gaussian99Percent2D(const Surface& s){return Flat(9.21, s);}
  TH2* Gaussian3Sigma2D   (const Surface& s){return Flat(11.83, s);}

  TH2* Gaussian68Percent1D(const Surface& s){return Flat(1.00, s);}
  TH2* Gaussian90Percent1D(const Surface& s){return Flat(2.71, s);}
  TH2* Gaussian95Percent1D(const Surface& s){return Flat(3.84, s);}
  TH2* Gaussian2Sigma1D   (const Surface& s){return Flat(4.00, s);}
  TH2* Gaussian99Percent1D(const Surface& s){return Flat(6.63, s);}
  TH2* Gaussian3Sigma1D   (const Surface& s){return Flat(9.00, s);}

  //----------------------------------------------------------------------
  void Surface::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();
    TObjString("Surface").Write("type");

    TVectorD v(3);
    v[0] = fMinChi;
    v[1] = fMinX;
    v[2] = fMinY;
    v.Write("minValues");

    fHist->Write("hist");

    TVectorD s(fSeedValues.size(), &fSeedValues[0]);
    s.Write("seeds");

    TDirectory* profDir = dir->mkdir("profHists");
    int idx = 0;
    for(auto it: fProfHists){
      profDir->cd();
      it->Write( TString::Format("hist%d", idx++));
    }

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr< Surface > Surface::LoadFrom(TDirectory* dir)
  {
    DontAddDirectory guard;

    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "Surface");

    const TVectorD v = *(TVectorD*)dir->Get("minValues");
    const TVectorD s = *(TVectorD*)dir->Get("seeds");

    std::unique_ptr<Surface> surf(new Surface);
    surf->fHist = (TH2F*)dir->Get("hist");
    surf->fMinChi = v[0];
    surf->fMinX = v[1];
    surf->fMinY = v[2];

    for(int idx = 0; idx < s.GetNrows(); ++idx){
      surf->fSeedValues.push_back(s[idx]);

      // Search for old "marg" name here too for backwards compatibility
      TH2* h = (TH2*)dir->Get(TString::Format("profHists/hist%d", idx));
      if(h)
        surf->fProfHists.push_back(h);
      else
        surf->fProfHists.push_back((TH2*)dir->Get(TString::Format("margHists/hist%d", idx)));
    }

    return surf;
  }


} // namespace
