#include "CAFAna/Core/Utilities.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/MathUtil.h"

#include "TArrayD.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TObjString.h"
#include "TString.h"
#include "TVector3.h"
#include "TVectorD.h"

#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include "sys/stat.h"
#include "wordexp.h"

namespace ana
{
  double LLPerBinFracSystErr::fgErr = -1;

  //----------------------------------------------------------------------
  IFDHSilent::IFDHSilent()
  {
    const char* s = getenv("IFDH_SILENT");
    fSet = (s && s == std::string("0"));
    if(!fSet) setenv("IFDH_SILENT", "1", 1);
    if(s) std::cout << "IFDH_SILENT=" << s << std::endl;
  }

  //----------------------------------------------------------------------
  IFDHSilent::~IFDHSilent()
  {
    if(!fSet) unsetenv("IFDH_SILENT");
  }

  //----------------------------------------------------------------------
  FloatingExceptionOnNaN::FloatingExceptionOnNaN(bool enable)
  {
    // Don't want any pending FPEs to trigger when we flip exceptions
    // on. Whoever had them off previously had a reason.
    feclearexcept(FE_INVALID);

    fegetexceptflag(&fBackup, FE_INVALID);

#ifndef DARWINBUILD
    if(enable)
      feenableexcept(FE_INVALID);
    else
      fedisableexcept(FE_INVALID);
#else
    std::cerr << "WARNING: CAFAna/Core/Utilities.cxx built on OS X, no feenableexcept available" << std::endl;
#endif
  }

  //----------------------------------------------------------------------
  FloatingExceptionOnNaN::~FloatingExceptionOnNaN()
  {
    fesetexceptflag(&fBackup, FE_INVALID);
  }

  //----------------------------------------------------------------------
  std::string Experiment()
  {
    const char* ret = getenv("EXPERIMENT");

    if(!ret){
      std::cout << "\nERROR: Environment variable $EXPERIMENT not set.\nThis is required for various ifdh/sam functionality.\nYou should probably setup sbndcode or icaruscode, though manually exporting the variable also seems to work." << std::endl;
      exit(1);
    }

    return ret;
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TMatrixD> CalcCovMx(const std::vector<TArrayD*> & binSets, int firstBin, int lastBin)
  {
    if (binSets.size() < 2)
      return std::unique_ptr<TMatrixD>(nullptr);

    if (lastBin < 0)
      lastBin = binSets[0]->GetSize() - 1;  // indexed from 0

    int nBins = lastBin - firstBin + 1;  // firstBin and lastBin are inclusive

    std::vector<double> binMeans(nBins);
    for( const auto & binSet : binSets )
    {
      for ( decltype(lastBin) binIdx = firstBin; binIdx <= lastBin; binIdx++ )
        binMeans[binIdx] += (*binSet)[binIdx];
    }
    for (decltype(lastBin) binIdx = firstBin; binIdx <= lastBin; binIdx++)
      binMeans[binIdx] /= binSets.size();


    auto covmx = std::make_unique<TMatrixD>(nBins, nBins);

    for( unsigned int hist_idx = 0; hist_idx < binSets.size(); ++hist_idx )
    {
      // first calculate the weighted sum of squares of the deviations
      for( decltype(nBins) i = 0; i < nBins; i++ )
      {
        double xi = (*(binSets[hist_idx]))[i];
        for( decltype(nBins) k = i; k < nBins; k++ )
        {
          double xk = (*(binSets[hist_idx]))[k];
          (*covmx)[i][k] += (xi - binMeans[i]) * (xk - binMeans[k]);
          if (i != k)
            (*covmx)[k][i] = (*covmx)[i][k];  // covariance matrices are always symmetric
        }
      }
    } // for (hist_idx)

    // now divide by N-1 to get sample covariance
    (*covmx) *= 1./(binSets.size()-1);

    return covmx;
  }

  //----------------------------------------------------------------------
  double LogLikelihood(double e, double o)
  {
    // http://www.wolframalpha.com/input/?i=d%2Fds+m*(1%2Bs)+-d+%2B+d*ln(d%2F(m*(1%2Bs)))%2Bs%5E2%2FS%5E2%3D0
    // http://www.wolframalpha.com/input/?i=solve+-d%2F(s%2B1)%2Bm%2B2*s%2FS%5E2%3D0+for+s
    const double S = LLPerBinFracSystErr::GetError();
    if(S > 0){
      const double S2 = util::sqr(S);
      const double s = .25*(sqrt(8*o*S2+util::sqr(e*S2-2))-e*S2-2);
      e *= 1+s;
    }

    // With this value, negative expected events and one observed
    // event gives a chisq from this one bin of 182.
    const double minexp = 1e-40; // Don't let expectation go lower than this

    assert(o >= 0);
    if(e < minexp){
      if(!o) return 0;
      e = minexp;
    }

    if(o*1000 > e){
      // This strange form is for numerical stability when e~o
      return 2*o*((e-o)/o + log1p((o-e)/e));
    }
    else{
      // But log1p doesn't like arguments near -1 (observation much smaller
      // than expectation), and it's better to use the usual formula in that
      // case.
      if(o){
        return 2*(e-o + o*log(o/e));
      }
      else{
        return 2*e;
      }
    }
  }

  //----------------------------------------------------------------------
  double LogLikelihood(const TH1* eh, const TH1* oh, bool useOverflow)
  {
    assert(eh->GetNbinsX() == oh->GetNbinsX());

    double chi = 0;

    int bufferBins = useOverflow? 2 : 1;

    for(int i = 0; i < eh->GetNbinsX()+bufferBins; ++i){
      const double e = eh->GetBinContent(i);
      const double o = oh->GetBinContent(i);

      chi += LogLikelihood(e, o);
    }

    return chi;
  }

  //----------------------------------------------------------------------
  double LogLikelihood(const Eigen::ArrayXd& ea, const Eigen::ArrayXd& oa, bool useOverflow)
  {
    assert(ea.size() == oa.size());

    double chi = 0;

    const int bufferBins = useOverflow ? 0 : -1;

    for(int i = 0; i < ea.size()+bufferBins; ++i){
      chi += LogLikelihood(ea[i], oa[i]);
    }

    return chi;
  }

  //----------------------------------------------------------------------
  double Chi2CovMx(const TVectorD* e, const TVectorD* o, const TMatrixD* covmxinv)
  {
    assert (e->GetNrows() == o->GetNrows());

    TVectorD diff = *o - *e;
    return diff * ((*covmxinv) * diff);  // operator* for two TVectorDs is the "dot product" (i.e., v1 * v2 = v1^{trans}v1)
  }

  //----------------------------------------------------------------------
  double Chi2CovMx(const TH1* e, const TH1* o, const TMatrixD* covmxinv)
  {
    TVectorD eVec(e->GetNbinsX());
    TVectorD oVec(o->GetNbinsX());
    for (int bin = 1; bin <= e->GetNbinsX(); bin++)
      eVec[bin-1] = e->GetBinContent(bin);
    for (int bin = 1; bin <= o->GetNbinsX(); bin++)
      oVec[bin-1] = o->GetBinContent(bin);

    return Chi2CovMx(&eVec, &oVec, covmxinv);
  }

  //----------------------------------------------------------------------
  TH2F* ExpandedHistogram(const std::string& title,
                          int nbinsx, double xmin, double xmax, bool xlog,
                          int nbinsy, double ymin, double ymax, bool ylog)
  {
    DontAddDirectory guard;

    if(xlog){xmin = log(xmin); xmax = log(xmax);}
    if(ylog){ymin = log(ymin); ymax = log(ymax);}

    // How wide the bins will be once we're done
    const double xwidth = (xmax-xmin)/(nbinsx-1);
    const double ywidth = (ymax-ymin)/(nbinsy-1);

    // Move the bin edges so that the limits occur at the centres
    xmin -= xwidth/2; ymin -= ywidth/2;
    xmax += xwidth/2; ymax += ywidth/2;

    std::vector<double> xedges(nbinsx+1);
    std::vector<double> yedges(nbinsy+1);

    for(int i = 0; i <= nbinsx; ++i){
      xedges[i] = xmin + (xmax-xmin)*i/double(nbinsx);
      if(xlog) xedges[i] = exp(xedges[i]);
    }
    for(int i = 0; i <= nbinsy; ++i){
      yedges[i] = ymin + (ymax-ymin)*i/double(nbinsy);
      if(ylog) yedges[i] = exp(yedges[i]);
    }

    return new TH2F(UniqueName().c_str(), title.c_str(),
                    nbinsx, &xedges.front(),
                    nbinsy, &yedges.front());
  }


  //----------------------------------------------------------------------
  std::unique_ptr<TMatrixD> SymmMxInverse(const TMatrixD& mx)
  {
    // check if there are any null rows/columns.
    // if there are, they make the matrix singular.
    // we will remove them temporarily,
    // invert the matrix, then put them back afterwards.
    std::set<int> nullRows;
    for (auto row = mx.GetRowLwb(); row <= mx.GetRowUpb(); row++)
    {
      bool rowIsNull = true;
      for (auto col = mx.GetColLwb(); col <= mx.GetColUpb(); col++)
      {
        if (mx[row][col] != 0.)
        {
          rowIsNull = false;
          break;
        }
      }

      if (rowIsNull)
        nullRows.insert(row);
    }

    std::cerr << " Notice: covariance matrix '" << mx.GetName() << "' has " << nullRows.size() << " null rows.\n"
        << "They will be removed before inverting and added back afterwards." << std::endl;

    // create a new matrix for inverting, skipping the null rows
    auto invMx = std::make_unique<TMatrixD>(mx.GetRowLwb(), mx.GetRowUpb() - nullRows.size(),
                                           mx.GetColLwb(), mx.GetColUpb() - nullRows.size());
    unsigned int skippedRows = 0;
    for (auto row = mx.GetRowLwb(); row <= mx.GetRowUpb(); row++)
    {
      if (nullRows.find(row) != nullRows.end())
      {
        skippedRows++;
        continue;
      }
      unsigned int skippedCols = 0;
      for (auto col = mx.GetColLwb(); col <= mx.GetColUpb(); col++)
      {
        // since we assumed the matrix is symmetric,
        // we can just use the null rows list here
        if (nullRows.find(col) != nullRows.end())
        {
          skippedCols++;
          continue;
        }

        (*invMx)[col-skippedCols][row-skippedRows] = (*invMx)[row-skippedRows][col-skippedCols] = mx[row][col];
      }
    }

    invMx->Invert();

    // put back the empty rows if there were any
    if (nullRows.size())
    {
      skippedRows = 0;
      auto retMx = std::make_unique<TMatrixD>(mx.GetRowLwb(), mx.GetRowUpb(),
                                              mx.GetColLwb(), mx.GetColUpb());
      for (auto row = mx.GetRowLwb(); row <= mx.GetRowUpb(); row++)
      {
        if (nullRows.find(row) != nullRows.end())
        {
          skippedRows++;
          continue;
        }

        unsigned int skippedCols = skippedRows;
        for (auto col = row; col <= mx.GetColUpb(); col++)
        {
          if (nullRows.find(col) != nullRows.end())
          {
            skippedCols++;
            continue;
          }

          (*retMx)[col][row] = (*retMx)[row][col] = (*invMx)[row-skippedRows][col-skippedCols];
        }
      }

      return retMx;
    }

    return invMx;
  }

  //----------------------------------------------------------------------
  std::string FindCAFAnaDir()
  {
    return std::string(getenv("MRB_SOURCE"))+"/sbncode/sbncode/CAFAna";
  }

  //----------------------------------------------------------------------
  std::vector<std::string> LoadFileList(const std::string& listfile)
  {
    std::vector<std::string> ret;

    std::ifstream is(listfile);
    if(!is.good()){
      std::cerr << "Can't open file list '" << listfile << "'. Aborting." << std::endl;
      abort();
    }

    while(!is.eof()){
      std::string fname;
      is >> fname;
      if(!fname.empty()) ret.push_back(fname);
    }
    return ret;
  }

  //----------------------------------------------------------------------
  std::map<std::string, std::string> GetCAFMetadata(TDirectory* dir)
  {
    std::map<std::string, std::string> ret;

    TIter next(dir->GetListOfKeys());
    while(TObject* key = next()){
      TObject* obj = dir->Get(key->GetName());
      assert(obj);

      const TString className = obj->ClassName();
      if(className == "TObjString"){
        ret[key->GetName()] = ((TObjString*)obj)->GetString();
      }
      else{
        std::cerr << "Unexpected object " << key->GetName() << " of type " << className << " while looking for metadata. Ignoring" << std::endl;
      }
    }

    return ret;
  }

  //----------------------------------------------------------------------
  void CombineMetadata(std::map<std::string, std::string>& base,
                       const std::map<std::string, std::string>& add,
                       std::set<std::string>& mask)
  {
    for(auto it: add){
      const std::string& key = it.first;

      // Needs special handling anyway, leave it blank.
      if(key == "parents") continue;

      // Accumulate the runs list
      if(key == "runs"){
        const std::string& r1 = base[key];
        const std::string& r2 = it.second;

        assert(!r2.empty());

        if(r1.empty()){
          base[key] = r2;
          continue;
        }

        // "[foo]" + "[bar]"
        std::string sum = r1+&r2[1]; // "[foo]bar]"
        sum[r1.size()-1] = ',';      // "[foo,bar]"
        base[key] = sum;
        continue;
      }

      if(base.find(key) == base.end()){
        // If it's new, add it
        base[key] = it.second;
      }
      else{
        if(key == "simulated.number_of_spills" ||
           key == "event_count" ||
           key == "online.totalevents"){
          // These two fields should be accumulated
          base[key] = TString::Format("%d",
                                      atoi(base[key].c_str()) +
                                      atoi(it.second.c_str())).Data();
        }
        else{
          // If it's a clash, record it
          if(base[key] != it.second) mask.insert(key);
        }
      }
    }
  }


  //----------------------------------------------------------------------
  void WriteCAFMetadata(TDirectory* dir,
                        const std::map<std::string, std::string>& meta)
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    for(auto it: meta){
      TObjString str(it.second.c_str());
      str.Write(it.first.c_str());
    }

    dir->Save();

    tmp->cd();
  }

  //----------------------------------------------------------------------
  bool RunningOnGrid()
  {
    static bool cache;
    static bool cache_set = false;
    if(!cache_set){
      cache = (getenv("_CONDOR_SCRATCH_DIR") != 0);
      cache_set = true;
    }

    return cache;
  }

  //----------------------------------------------------------------------
  size_t Stride(bool allow_default)
  {
    static int cache = -1;

    if(cache < 0){
      char* env = getenv("CAFANA_STRIDE");
      if(env){
        cache = std::atoi(env);
      }
      else{
        if(allow_default){
          cache = 1;
        }
        else{
          std::cout << "Stride() called, but CAFANA_STRIDE is not set (--stride not passed?)" << std::endl;
          abort();
        }
      }
    }

    return cache;
  }

  //----------------------------------------------------------------------
  size_t Offset(bool allow_default)
  {
    static int cache = -1;

    if(cache < 0){
      char* env = getenv("CAFANA_OFFSET");
      if(env){
        cache = std::atoi(env);
      }
      else{
        if(allow_default){
          cache = 0;
        }
        else{
          std::cout << "Offset() called, but CAFANA_OFFSET is not set (--offset not passed?)" << std::endl;
          abort();
        }
      }
    }

    return cache;
  }

  //----------------------------------------------------------------------
  int Limit()
  {
    static int cache = 0;

    if(cache == 0){
      char* env = getenv("CAFANA_LIMIT");
      if(env){
        cache = std::atoi(env);
      }
      else{
        cache = -1;
      }
    }

    return cache;
  }

  //----------------------------------------------------------------------
  size_t JobNumber()
  {
    if(!RunningOnGrid()){
      std::cout << "JobNumber() called, but we are not running on the grid" << std::endl;
      abort();
    }

    return Offset(false);
  }

  //----------------------------------------------------------------------
  size_t NumJobs()
  {
    if(!RunningOnGrid()){
      std::cout << "NumJobs() called, but we are not running on the grid" << std::endl;
      abort();
    }

    return Stride(false);
  }

  //----------------------------------------------------------------------
  FitToFourier::FitToFourier(TH1* h, double xlo, double xhi, int NOsc)
    : fHist(h), fxlo(xlo), fxhi(xhi), fNOsc(NOsc)
  {
  }

  //----------------------------------------------------------------------
  FitToFourier::~FitToFourier()
  {
  }

  //----------------------------------------------------------------------
  double FitToFourier::operator()(double *x, double *par) const
  {
    double x0 = x[0];
    double val = par[0];
    for (int i = 1; i <= fNOsc; i++)
      val += par[2*i-1]*sin(i*M_PI*x0) + par[2*i]*cos(i*M_PI*x0);
    return val;
  }

  //----------------------------------------------------------------------
  TF1* FitToFourier::Fit() const
  {
    //double s[fNOsc] = {0};
    //double c[fNOsc] = {0};

    std::vector<double> s(fNOsc, 0.0);
    std::vector<double> c(fNOsc, 0.0);

    int nBins = 0;
    for(int i = 1; i <= fHist->GetNbinsX(); ++i){
      const double x = M_PI * fHist->GetXaxis()->GetBinCenter(i);
      const double y = fHist->GetBinContent(i);

      if(y == 0) continue;
      ++nBins;

      for(int n = 0; n <= fNOsc; ++n){
        s[n] += y * sin(n*x);
        c[n] += y * cos(n*x);
      }
    }

    for(int n = 0; n <= fNOsc; ++n){
      s[n] *= 2./nBins;
      c[n] *= 2./nBins;
    }

    TF1* f = new TF1(UniqueName().c_str(), this, fxlo, fxhi, 2*fNOsc+1);

    f->SetParameter(0, c[0]/2);
    for(int n = 1; n <= fNOsc; ++n){
      f->SetParameter(n*2-1, s[n]);
      f->SetParameter(n*2,   c[n]);
    }

    // Because ROOT is having problems drawing f if I don't
    double min = fHist->GetMinimum();
    double max = fHist->GetMaximum();
    f->GetYaxis()->SetRangeUser(0.8*min, 1.2*max);
    return f;
  }

  //----------------------------------------------------------------------
  void EnsurePositiveDefinite(TH2* mat)
  {
    // Convert histogram to a proper matrix
    assert(mat->GetNbinsX() == mat->GetNbinsY());
    const int N = mat->GetNbinsX();
    TMatrixD m(N, N);
    for(int i = 0; i < N; ++i)
      for(int j = 0; j < N; ++j)
        m(i, j) = mat->GetBinContent(i+1, j+1);

    // Decompose it
    TVectorD evals;
    TMatrixD evecs = m.EigenVectors(evals);
    TMatrixD evalmat(N, N);
    // Force any negative eigenvalues slightly positive (floating point errors)
    for(int i = 0; i < N; ++i) evalmat(i, i) = std::max(1e-14, evals[i]);

    // Put the original matrix back together
    const TMatrixD evecs_inv(TMatrixD::kTransposed, evecs);
    m = evecs*evalmat*evecs_inv;

    // Decompose again to check for floating point problems
    m.EigenVectors(evals);
    for(int i = 0; i < N; ++i) assert(evals[i] > 0);

    // Copy the new matrix contents back into the histogram
    for(int i = 0; i < N; ++i)
      for(int j = 0; j < N; ++j)
        mat->SetBinContent(i+1, j+1, m(i, j));
  }

  //----------------------------------------------------------------------
  // Note that this does not work for 3D!
  Eigen::ArrayXd GetMaskArray(const Spectrum& s, double xmin, double xmax, double ymin, double ymax)
  {
    if (s.NDimensions() > 2){
      std::cout << "Error: unable to apply a mask in " << s.GetBinnings().size() << " dimensions" << std::endl;
      abort();
    }

    if(s.NDimensions() == 1 && ymax > ymin){
      std::cout << "Error: GetMaskArray(): can't specify y range for 1D spectrum" << std::endl;
      abort();
    }

    const Binning* xbins = &s.GetBinnings()[0];
    const Binning* ybins = (s.NDimensions() == 2) ? &s.GetBinnings()[1] : 0;

    const int Nx = xbins->NBins();
    const int Ny = ybins ? ybins->NBins() : 1;

    // The 1D flattening of 2D binning is pretty confusing. The bins are packed
    // densely, without under/overflow, *except* there is a single underflow at
    // 0 and single overflow at Nx*Ny+1. So we do our calculations as if there
    // were no under/overflow and then add 1 to the output index to account.

    Eigen::ArrayXd ret(Nx*Ny+2);

    // Include underflow and overflow if mask disabled, otherwise exclude
    ret[0] = ret[Nx*Ny+1] = ((xmin < xmax || ymin < ymax) ? 0 : 1);

    for(int i = 0; i < Nx*Ny; ++i){

      const int ix = i / Ny;
      const int iy = i % Ny;

      bool isMask = false;

      if (xmin < xmax){
	if (xbins->Edges()[ix  ] < xmin) isMask = true;
	if (xbins->Edges()[ix+1] > xmax) isMask = true;
      }

      if (ymin < ymax){
	if (ybins->Edges()[iy  ] < ymin) isMask = true;
	if (ybins->Edges()[iy+1] > ymax) isMask = true;
      }

      ret[i+1] = isMask ? 0 : 1;
    }

    return ret;
  }
}
