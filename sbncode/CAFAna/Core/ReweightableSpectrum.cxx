#include "CAFAna/Core/ReweightableSpectrum.h"

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Core/Utilities.h"

#include "TDirectory.h"
#include "TH2.h"
#include "TObjString.h"

#include <cassert>
#include <iostream>
#include <memory>

namespace ana
{
  //----------------------------------------------------------------------
  ReweightableSpectrum::ReweightableSpectrum(SpectrumLoaderBase& loader,
                                             const HistAxis& recoAxis,
                                             const HistAxis& trueAxis,
                                             const Cut& cut,
                                             const SystShifts& shift,
                                             const Var& wei)
    : ReweightableSpectrum(recoAxis.GetLabels(), recoAxis.GetBinnings(),
                           trueAxis.GetVars()[0])
  {
    assert(trueAxis.NDimensions() == 1);

    fTrueLabel = trueAxis.GetLabels()[0];

    DontAddDirectory guard;

    // Can't use HistCache here because y-axis is not necessarily
    // TrueEnergyBinning. TODO - that should maybe be generalized.

    const std::string name = UniqueName();

    Binning xbins = fBins[0];
    if(fBins.size() > 1){
      int n = 1;
      for(const Binning& b: fBins) n *= b.NBins();
      xbins = Binning::Simple(n, 0, n);
    }

    const Binning ybins = trueAxis.GetBinnings()[0];


    // Ugh combinatorics
    if(xbins.IsSimple() && ybins.IsSimple())
      fHist =  new TH2D(name.c_str(), "",
                        xbins.NBins(), xbins.Min(), xbins.Max(),
                        ybins.NBins(), ybins.Min(), ybins.Max());

    if(xbins.IsSimple() && !ybins.IsSimple())
      fHist =  new TH2D(name.c_str(), "",
                        xbins.NBins(), xbins.Min(), xbins.Max(),
                        ybins.NBins(), &ybins.Edges()[0]);

    if(!xbins.IsSimple() && ybins.IsSimple())
      fHist =  new TH2D(name.c_str(), "",
                        xbins.NBins(), &xbins.Edges()[0],
                        ybins.NBins(), ybins.Min(), ybins.Max());

    if(!xbins.IsSimple() && !ybins.IsSimple())
      fHist =  new TH2D(name.c_str(), "",
                        xbins.NBins(), &xbins.Edges()[0],
                        ybins.NBins(), &ybins.Edges()[0]);

    loader.AddReweightableSpectrum(*this, recoAxis.GetMultiDVar(), cut, shift, wei);
  }

  //----------------------------------------------------------------------
  ReweightableSpectrum::ReweightableSpectrum(const Var& rwVar,
                                             const std::string& xlabel, const std::string& ylabel,
                                             double pot,
                                             int nbinsx, double xmin, double xmax,
                                             int nbinsy, double ymin, double ymax)
    : ReweightableSpectrum(xlabel,
                           Binning::Simple(nbinsx, xmin, xmax),
                           rwVar)
  {
    DontAddDirectory guard;

    fHist = new TH2D(UniqueName().c_str(), "",
                     nbinsx, xmin, xmax, nbinsy, ymin, ymax);

    fTrueLabel = ylabel;

    // Ensure errors get accumulated properly
    fHist->Sumw2();
  }

  //----------------------------------------------------------------------
  ReweightableSpectrum::ReweightableSpectrum(const Var& rwVar, TH2* h,
                                             const std::vector<std::string>& labels,
                                             const std::vector<Binning>& bins,
                                             double pot, double livetime)
    : ReweightableSpectrum(labels, bins, rwVar)
  {
    fPOT = pot;
    fLivetime = livetime;

    if(!h){
      fHist = 0;
      return;
    }

    DontAddDirectory guard;

    const TString className = h->ClassName();

    if(className == "TH2D"){
      // Shortcut if types match
      fHist = new TH2D(*((TH2D*)h));
    }
    else{
      const TAxis* ax = h->GetXaxis();
      const TAxis* ay = h->GetYaxis();
      // Must have both or neither
      assert(bool(ax->GetXbins()->GetArray()) == bool(ay->GetXbins()->GetArray()));

      if(ax->GetXbins()->GetArray()){
        fHist = new TH2D(UniqueName().c_str(), "",
                         ax->GetNbins(), ax->GetXbins()->GetArray(),
                         ay->GetNbins(), ay->GetXbins()->GetArray());
      }
      else{
        fHist = new TH2D(UniqueName().c_str(), "",
                         ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                         ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
      }

      fHist->Add(h);
    }

    fTrueLabel = h->GetYaxis()->GetTitle();
  }

    //----------------------------------------------------------------------
  ReweightableSpectrum::ReweightableSpectrum(const Var& rwVar,
                                             std::unique_ptr<TH2D> h,
                                             const std::vector<std::string>& labels,
                                             const std::vector<Binning>& bins,
                                             double pot, double livetime)
    : ReweightableSpectrum(labels, bins, rwVar)
  {
    fHist = h.release();
    fPOT = pot;
    fLivetime = livetime;

    fTrueLabel = fHist->GetYaxis()->GetTitle();
  }

  //----------------------------------------------------------------------
  ReweightableSpectrum::~ReweightableSpectrum()
  {
    if(fHist && fHist->GetDirectory()){
      static bool once = true;
      if(once){
        once = false;
        std::cerr << "ReweightableSpectrum's fHist is associated with a directory. How did that happen?" << std::endl;
      }
    }

    // Can't use HistCache here because that expects all 2D spectra have the
    // usual true energy binning on their y-axis. In the case that we're
    // actually an OscillatableSpectrum that destructor will indeed take care
    // of invoking HistCache.
    delete fHist;
  }

  //----------------------------------------------------------------------
  ReweightableSpectrum::ReweightableSpectrum(const ReweightableSpectrum& rhs)
    : fRWVar(rhs.fRWVar), fLabels(rhs.fLabels), fBins(rhs.fBins)
  {
    DontAddDirectory guard;

    fHist = new TH2D(*rhs.fHist);

    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;

    assert( rhs.fLoaderCount.empty() ); // Copying with pending loads is unexpected
  }

  //----------------------------------------------------------------------
  ReweightableSpectrum& ReweightableSpectrum::operator=(const ReweightableSpectrum& rhs)
  {
    if(this == &rhs) return *this;

    DontAddDirectory guard;

    fRWVar = rhs.fRWVar;
    fLabels = rhs.fLabels;
    fBins = rhs.fBins;

    delete fHist;
    fHist = new TH2D(*rhs.fHist);
    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;

    assert( fLoaderCount.empty() ); // Copying with pending loads is unexpected

    return *this;
  }

  //----------------------------------------------------------------------
  TH2D* ReweightableSpectrum::ToTH2(double pot) const
  {
    // Could have a file temporarily open
    DontAddDirectory guard;

    TH2D* ret = new TH2D(*fHist);
    if(fPOT){
      ret->Scale(pot/fPOT);
    }
    else{
      // How did it get events with no POT?
      assert(ret->Integral() == 0);
    }


    std::string label;
    for(const std::string& l: fLabels) label += l + " and ";
    label.resize(label.size()-5); // drop the last "and"
    ret->GetXaxis()->SetTitle(label.c_str());
    ret->GetYaxis()->SetTitle(fTrueLabel.c_str());

    return ret;
  }

  //----------------------------------------------------------------------
  void ReweightableSpectrum::Fill(double x, double y, double w)
  {
    fHist->Fill(x, y, w);
  }

  /// Helper for \ref Unweighted
  void ProjectionX(TH2D* from, TH1D* to)
  {
    const int Nx = from->GetNbinsX();
    const int Ny = from->GetNbinsY();

    // Direct access to the bins is faster
    double* fromArr = from->GetArray();
    double* toArr = to->GetArray();

    int fromBin = 0;
    for(int y = 0; y < Ny+2; ++y){
      for(int x = 0; x < Nx+2; ++x){
        // Our loops go over the bins in the order they are internally in
        // 'from', and we do overflows, so we keep up exactly. If you get
        // paranoid, reenable this briefly.
        //        assert(fromBin == from->GetBin(x, y));

        const double z = fromArr[fromBin];
        ++fromBin;

        toArr[x] += z;
      }
    }
  }

  //----------------------------------------------------------------------
  Spectrum ReweightableSpectrum::UnWeighted() const
  {
    DontAddDirectory guard;

    // Create a suitably-sized space for the result
    std::unique_ptr<TH1D> h(HistCache::New("", fHist->GetXaxis()));

    ProjectionX(fHist, h.get());

    return Spectrum(std::move(h), fLabels, fBins, fPOT, fLivetime);
  }

  /// Helper for \ref WeightingVariable
  void ProjectionY(TH2D* from, TH1D* to)
  {
    const int Nx = from->GetNbinsX();
    const int Ny = from->GetNbinsY();

    // Direct access to the bins is faster
    double* fromArr = from->GetArray();
    double* toArr = to->GetArray();

    int fromBin = 0;
    for(int y = 0; y < Ny+2; ++y){
      for(int x = 0; x < Nx+2; ++x){
        // Our loops go over the bins in the order they are internally in
        // 'from', and we do overflows, so we keep up exactly. If you get
        // paranoid, reenable this briefly.
        //        assert(fromBin == from->GetBin(x, y));

        const double z = fromArr[fromBin];
        ++fromBin;

        toArr[y] += z;
      }
    }
  }

  //----------------------------------------------------------------------
  Spectrum ReweightableSpectrum::WeightingVariable() const
  {
    DontAddDirectory guard;

    // Create a suitably-sized space for the result
    std::unique_ptr<TH1D> h(HistCache::New("", fHist->GetYaxis()));

    ProjectionY(fHist, h.get());

    return Spectrum(std::move(h), fLabels, fBins, fPOT, fLivetime);
  }

  //----------------------------------------------------------------------
  Spectrum ReweightableSpectrum::WeightedBy(const TH1* ws) const
  {
    // This function is in the inner loop of oscillation fits, so some
    // optimization has been done.

    DontAddDirectory guard;

    assert(ws->GetNbinsX() == fHist->GetNbinsY());

    TAxis* ax = fHist->GetXaxis();
    TH1D* hRet = HistCache::New("", ax);

    const int X = fHist->GetNbinsX();
    const int Y = fHist->GetNbinsY();

    // Direct access to the bins is faster
    double* retArr = hRet->GetArray();
    double* histArr = fHist->GetArray();

    int bin = 0;
    for(int y = 0; y < Y+2; ++y){
      const double w = ws->GetBinContent(y);
      for(int x = 0; x < X+2; ++x){
        // Our loops go over the bins in the order they are internally in
        // fHist, and we do overflows, so we keep up exactly. If you get
        // paranoid, reenable this briefly.

        // assert(bin == fHist->GetBin(x, y));

        retArr[x] += histArr[bin]*w;
        ++bin;
      }
    }

    // TODO: can this all be more efficient?
    return Spectrum(std::unique_ptr<TH1D>(hRet), fLabels, fBins, fPOT, fLivetime);
  }


  //----------------------------------------------------------------------
  void ReweightableSpectrum::ReweightToTrueSpectrum(const Spectrum& target)
  {
    // This is a big component of what extrapolations do, so it has been
    // optimized for speed

    Ratio corr(target, WeightingVariable());
    std::unique_ptr<TH1D> hcorr(corr.ToTH1());

    assert(hcorr->GetNbinsX() == fHist->GetNbinsY());

    const int X = fHist->GetNbinsX();
    const int Y = fHist->GetNbinsY();

    // Direct access to the bins is faster
    double* histArr = fHist->GetArray();
    double* corrArr = hcorr->GetArray();

    int bin = 0;
    for(int y = 0; y < Y+2; ++y){
      const double w = corrArr[y];
      for(int x = 0; x < X+2; ++x){
        // Our loops go over the bins in the order they are internally in
        // fHist, and we do overflows, so we keep up exactly. If you get
        // paranoid, reenable this briefly.

        // assert(bin == fHist->GetBin(x, y));

        histArr[bin] *= w;
        ++bin;
      }
    }

    TH1D* todel = hcorr.release();
    HistCache::Delete(todel);
  }

  //----------------------------------------------------------------------
  void ReweightableSpectrum::ReweightToRecoSpectrum(const Spectrum& target)
  {
    // This is a big component of what extrapolations do, so it has been
    // optimized for speed

    Ratio corr(target, UnWeighted());
    std::unique_ptr<TH1D> hcorr(corr.ToTH1());

    assert(hcorr->GetNbinsX() == fHist->GetNbinsX());

    const int X = fHist->GetNbinsX();
    const int Y = fHist->GetNbinsY();

    // Direct access to the bins is faster
    double* histArr = fHist->GetArray();
    double* corrArr = hcorr->GetArray();

    int bin = 0;
    for(int y = 0; y < Y+2; ++y){
      for(int x = 0; x < X+2; ++x){
        // Our loops go over the bins in the order they are internally in
        // fHist, and we do overflows, so we keep up exactly. If you get
        // paranoid, reenable this briefly.

        // assert(bin == fHist->GetBin(x, y));

        histArr[bin] *= corrArr[x];
        ++bin;
      }
    }

    TH1D* todel = hcorr.release();
    HistCache::Delete(todel);
  }

  //----------------------------------------------------------------------
  void ReweightableSpectrum::Clear()
  {
    fHist->Reset();
  }

  //----------------------------------------------------------------------
  void ReweightableSpectrum::RemoveLoader(SpectrumLoaderBase* p)
  { fLoaderCount.erase(p); }

  //----------------------------------------------------------------------
  void ReweightableSpectrum::AddLoader(SpectrumLoaderBase* p)
  { fLoaderCount.insert(p); }

  //----------------------------------------------------------------------
  void ReweightableSpectrum::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    TObjString("ReweightableSpectrum").Write("type");

    fHist->GetYaxis()->SetTitle(fTrueLabel.c_str());
    fHist->Write("hist");
    TH1D hPot("", "", 1, 0, 1);
    hPot.Fill(.5, fPOT);
    hPot.Write("pot");
    TH1D hLivetime("", "", 1, 0, 1);
    hLivetime.Fill(.5, fLivetime);
    hLivetime.Write("livetime");

    for(unsigned int i = 0; i < fBins.size(); ++i){
      TObjString(fLabels[i].c_str()).Write(TString::Format("label%d", i).Data());
      fBins[i].SaveTo(dir->mkdir(TString::Format("bins%d", i)));
    }

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<ReweightableSpectrum> ReweightableSpectrum::LoadFrom(TDirectory* dir)
  {
    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "ReweightableSpectrum");

    TH2D* spect = (TH2D*)dir->Get("hist");
    assert(spect);
    TH1* hPot = (TH1*)dir->Get("pot");
    assert(hPot);
    TH1* hLivetime = (TH1*)dir->Get("livetime");
    assert(hLivetime);

    std::vector<std::string> labels;
    std::vector<Binning> bins;

    for(int i = 0; ; ++i){
      TDirectory* subdir = dir->GetDirectory(TString::Format("bins%d", i));
      if(!subdir) break;
      bins.push_back(*Binning::LoadFrom(subdir));
      TObjString* label = (TObjString*)dir->Get(TString::Format("label%d", i));
      labels.push_back(label ? label->GetString().Data() : "");
    }

    if(bins.empty() && labels.empty()){
      // Must be an old file. Make an attempt at backwards compatibility.
      bins.push_back(Binning::FromTAxis(spect->GetXaxis()));
      labels.push_back(spect->GetXaxis()->GetTitle());
    }

    return std::make_unique<ReweightableSpectrum>(kUnweighted,
                                                  spect,
                                                  labels, bins,
                                                  hPot->GetBinContent(1),
                                                  hLivetime->GetBinContent(1));
  }

}
