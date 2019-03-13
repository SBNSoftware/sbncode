#include "CAFAna/Core/Binning.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TObjString.h"
#include "TVectorD.h"

namespace ana
{
  int Binning::fgNextID = 0;

  //----------------------------------------------------------------------
  Binning::Binning()
    : fID(-1)
  {
  }

  //----------------------------------------------------------------------
  Binning Binning::SimpleHelper(int n, double lo, double hi,
                                const std::vector<std::string>& labels)
  {
    assert(labels.empty() || int(labels.size()) == n);

    Binning bins;
    bins.fNBins = n;
    bins.fMin = lo;
    bins.fMax = hi;
    bins.fEdges.resize(n+1);
    for (int i = 0; i <= n; i++)
      bins.fEdges[i] = lo + i*(hi-lo)/n;
    bins.fLabels = labels;
    bins.fIsSimple = true;

    return bins;
  }

  //----------------------------------------------------------------------
  Binning Binning::Simple(int n, double lo, double hi,
                          const std::vector<std::string>& labels)
  {
    Binning bins = SimpleHelper(n, lo, hi, labels);

    auto it = IDMap().find(bins);
    if(it == IDMap().end()){
      bins.fID = fgNextID++;
      IDMap().emplace(bins, bins.fID);
    }
    else{
      bins.fID = it->second;
    }

    return bins;
  }

  //----------------------------------------------------------------------
  Binning Binning::LogUniform(int n, double lo, double hi)
  {
    std::vector<double> edges(n+1);
    const double logSpacing = exp( (log(hi) - log(lo)) / n );
    for (int i = 0; i <= n; ++i) {
      if (i == 0) edges[i] = lo;
      else        edges[i] = logSpacing*edges[i-1];
    }
    return Custom(edges);
  }

  //----------------------------------------------------------------------
  Binning Binning::CustomHelper(const std::vector<double>& edges)
  {
    assert(edges.size() > 1);

    Binning bins;
    bins.fEdges = edges;
    bins.fNBins = edges.size()-1;
    bins.fMin = edges.front();
    bins.fMax = edges.back();
    bins.fIsSimple = false;

    return bins;
  }

  //----------------------------------------------------------------------
  Binning Binning::Custom(const std::vector<double>& edges)
  {
    Binning bins = CustomHelper(edges);

    auto it = IDMap().find(bins);
    if(it == IDMap().end()){
      bins.fID = fgNextID++;
      IDMap().emplace(bins, bins.fID);
    }
    else{
      bins.fID = it->second;
    }

    return bins;
  }

  //----------------------------------------------------------------------
  int Binning::FindBin(float x) const
  {
    // Treat anything outside [fMin, fMax) at Underflow / Overflow
    if (x <  fMin) return 0;               // Underflow
    if (x >= fMax) return fEdges.size();   // Overflow

    // Follow ROOT convention, first bin of histogram is bin 1

    if (this->IsSimple()){
      double binwidth = (fMax - fMin) / fNBins;
      int bin = (x - fMin) / binwidth + 1;
      return bin;
    }

    int bin =
      std::lower_bound(fEdges.begin(), fEdges.end(), x) - fEdges.begin();
    if (x == fEdges[bin]) bin++;
    assert(bin >= 0 && bin < (int)fEdges.size());
    return bin;
  }

  //----------------------------------------------------------------------
  Binning Binning::FromTAxis(const TAxis* ax)
  {
    Binning bins;

    // Evenly spaced binning
    if(!ax->GetXbins()->GetArray()){
      bins = SimpleHelper(ax->GetNbins(), ax->GetXmin(), ax->GetXmax());
    }
    else{
      std::vector<double> edges(ax->GetNbins()+1);
      ax->GetLowEdge(&edges.front());
      edges[ax->GetNbins()] = ax->GetBinUpEdge(ax->GetNbins());

      bins = Binning::Custom(edges);
    }

    auto it = IDMap().find(bins);
    if(it != IDMap().end()){
      bins.fID = it->second;
    }
    else{
      bins.fID = fgNextID++;
      IDMap().emplace(bins, bins.fID);
    }

    return bins;
  }

  // Can we give trueE bin a special ID?
  //----------------------------------------------------------------------
  Binning TrueEnergyBins()
  {
    // Osc P is ~sin^2(1/E). Difference in prob across a bin is ~dP/dE which
    // goes as 1/E^2 times a trigonometric term depending on the parameters but
    // bounded between -1 and +1. So bins should have width ~1/E^2. E_i~1/i
    // achieves that.

    const int kNumTrueEnergyBins = 100;

    // N+1 bin low edges
    std::vector<double> edges(kNumTrueEnergyBins+1);

    const double Emin = .5; // 500 MeV: there's really no events below there

    // How many edges to generate. Allow room for 0-Emin bin
    const double N = kNumTrueEnergyBins-1;
    const double A = N*Emin;

    edges[0] = 0;

    for(int i = 1; i <= N; ++i){
      edges[kNumTrueEnergyBins-i] = A/i;
    }

    edges[kNumTrueEnergyBins] = 120; // Replace the infinity that would be here

    return Binning::Custom(edges);
  }


  //----------------------------------------------------------------------
  Binning TrueLOverTrueEBins()
  {
    // constant binnig

    const int kNumTrueLOverTrueEBins = 2000;
    const double klow = 0.0;
    const double khigh = 5.0;

    return Binning::Simple(kNumTrueLOverTrueEBins, klow, khigh);
  }

  //----------------------------------------------------------------------
  void Binning::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    TObjString("Binning").Write("type");

    TVectorD nminmax(3);

    nminmax[0] = fNBins;
    nminmax[1] = fMin;
    nminmax[2] = fMax;

    nminmax.Write("nminmax");

    TVectorD issimple(1);
    issimple[0] = fIsSimple;
    issimple.Write("issimple");

    TVectorD edges(fEdges.size());
    for(unsigned int i = 0; i < fEdges.size(); ++i)
      edges[i] = fEdges[i];

    edges.Write("edges");

    for(unsigned int i = 0; i < fLabels.size(); ++i)
      TObjString(fLabels[i].c_str()).Write(TString::Format("label%d", i).Data());

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<Binning> Binning::LoadFrom(TDirectory* dir)
  {
    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "Binning");

    TVectorD* vMinMax = (TVectorD*)dir->Get("nminmax");
    assert(vMinMax);

    Binning ret;

    const TVectorD* issimple = (TVectorD*)dir->Get("issimple");
    if((*issimple)[0]){
      ret = Binning::Simple((*vMinMax)[0],
                            (*vMinMax)[1],
                            (*vMinMax)[2]);
    }
    else{
      const TVectorD* vEdges = (TVectorD*)dir->Get("edges");
      std::vector<double> edges(vEdges->GetNrows());
      for(int i = 0; i < vEdges->GetNrows(); ++i) edges[i] = (*vEdges)[i];

      ret = Binning::Custom(edges);
    }

    for(unsigned int i = 0; ; ++i){
      TObjString* s = (TObjString*)dir->Get(TString::Format("label%d", i).Data());
      if(!s) break;
      ret.fLabels.push_back(s->GetString().Data());
    }

    return std::make_unique<Binning>(ret);
  }

  //----------------------------------------------------------------------
  bool Binning::operator==(const Binning& rhs) const
  {
    // NB don't look at ID here or in < because we use these in the maps below
    // that are used to find the IDs in the first place
    if(fIsSimple != rhs.fIsSimple) return false;
    if(fIsSimple){
      return fNBins == rhs.fNBins && fMin == rhs.fMin && fMax == rhs.fMax;
    }
    else{
      return fEdges == rhs.fEdges;
    }
  }

  //----------------------------------------------------------------------
  bool Binning::operator<(const Binning& rhs) const
  {
    if(fIsSimple != rhs.fIsSimple) return fIsSimple < rhs.fIsSimple;
    if(fIsSimple){
      return std::make_tuple(fNBins, fMin, fMax) < std::make_tuple(rhs.fNBins, rhs.fMin, rhs.fMax);
    }
    else{
      return fEdges < rhs.fEdges;
    }
  }

  //----------------------------------------------------------------------
  std::map<Binning, int>& Binning::IDMap()
  {
    static std::map<Binning, int> ret;
    return ret;
  }

}
