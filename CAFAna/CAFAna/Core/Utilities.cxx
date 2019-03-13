#include "CAFAna/Core/Utilities.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Ratio.h"

#include "Utilities/func/MathUtil.h"

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

#include <cassert>
#include <cmath>
#include <fstream>
#include "sys/stat.h"
#include "wordexp.h"

namespace ana
{
  double LLPerBinFracSystErr::fgErr = -1;

  //----------------------------------------------------------------------
  std::string UniqueName()
  {
    static int N = 0;
    return TString::Format("cafanauniq%d", N++).Data();
  }

  //----------------------------------------------------------------------
  DontAddDirectory::DontAddDirectory()
  {
    fBackup = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
  }

  //----------------------------------------------------------------------
  DontAddDirectory::~DontAddDirectory()
  {
    TH1::AddDirectory(fBackup);
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

    double chi = 0;

    assert(o >= 0);
    if(e < minexp){
      if(!o) return 0;
      e = minexp;
    }

    chi += 2*(e-o);
    if(o) chi += 2*o*log(o/e);

    return chi;
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
  double LogLikelihoodDerivative(double e, double o, double dedx)
  {
    double ret = 2*dedx;
    if(o) ret -= 2*o*dedx / e;
    return ret;
  }

  //----------------------------------------------------------------------
  double LogLikelihoodDerivative(const TH1D* eh, const TH1D* oh,
                                 const std::vector<double>& dedx)
  {
    const double* ea = eh->GetArray();
    const double* oa = oh->GetArray();

    double ret = 0;
    for(unsigned int i = 0; i < dedx.size(); ++i){
      ret += LogLikelihoodDerivative(ea[i], oa[i], dedx[i]);
    }
    return ret;
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
                          int nbinsx, double xmin, double xmax,
                          int nbinsy, double ymin, double ymax)
  {
    DontAddDirectory guard;

    // How wide the bins will be once we're done
    const double xwidth = (xmax-xmin)/(nbinsx-1);
    const double ywidth = (ymax-ymin)/(nbinsy-1);

    // Move the bin edges so that the limits occur at the centres
    xmin -= xwidth/2; ymin -= ywidth/2;
    xmax += xwidth/2; ymax += ywidth/2;

    return new TH2F(UniqueName().c_str(), title.c_str(),
                    nbinsx, xmin, xmax,
                    nbinsy, ymin, ymax);
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

  // Helper functions for MakeTHND().
  namespace{
    // Eventually the bin parameters will all be unpacked and we just pass them
    // on to the regular constructor.
    template<class T, class... A> T* MakeHist(A... a)
    {
      DontAddDirectory guard;
      return new T(a...);
    }

    // This function consumes bins from the start of the argument list and
    // pushes their translations onto the list of arguments at the end.
    template<class T, class... A> T* MakeHist(const Binning& firstBin,
                                              A... args)
    {
      if(firstBin.IsSimple())
        return MakeHist<T>(args...,
                           firstBin.NBins(), firstBin.Min(), firstBin.Max());
      else
        return MakeHist<T>(args...,
                           firstBin.NBins(), &firstBin.Edges().front());
    }
  }

  // Concrete instantiations. MakeHist() requires us to put the bin arguments
  // first...
  //----------------------------------------------------------------------
  TH1D* MakeTH1D(const char* name, const char* title, const Binning& bins)
  {
    return MakeHist<TH1D>(bins, name, title);
  }

  //----------------------------------------------------------------------
  TH2D* MakeTH2D(const char* name, const char* title,
                 const Binning& binsx,
                 const Binning& binsy)
  {
    return MakeHist<TH2D>(binsx, binsy, name, title);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, ana::EBinType bintype)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(s.ToTH1(exposure, expotype));
    return ToTH2Helper(std::move(h1), binsx, binsy, bintype);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2(const Ratio& r,
             const Binning& binsx, const Binning& binsy)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(r.ToTH1());
    return ToTH2Helper(std::move(h1), binsx, binsy);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx, const Binning& binsy,
		   ana::EBinType bintype)
  {
    // Make sure it's compatible with having been made with this binning
    assert(h1->GetNbinsX() == binsx.NBins()*binsy.NBins());

    TH2* ret = MakeTH2D("", UniqueName().c_str(), binsx, binsy);

    for(int i = 0; i < h1->GetNbinsX(); ++i){
      const double val = h1->GetBinContent(i+1);
      const double err = h1->GetBinError(i+1);

      const int ix = i / binsy.NBins();
      const int iy = i % binsy.NBins();

      ret->SetBinContent(ix+1, iy+1, val);
      ret->SetBinError  (ix+1, iy+1, err);
    }

    if(bintype == ana::EBinType::kBinDensity) ret->Scale(1, "width");

    return ret;
  }

  //----------------------------------------------------------------------

  TH3* ToTH3(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, const Binning& binsz,
	     ana::EBinType bintype)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(s.ToTH1(exposure, expotype));

    return ToTH3Helper(std::move(h1), binsx, binsy, binsz, bintype);
  }

  //----------------------------------------------------------------------

  TH3* ToTH3(const Ratio& r,
             const Binning& binsx, const Binning& binsy, const Binning& binsz)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(r.ToTH1());

    return ToTH3Helper(std::move(h1), binsx, binsy, binsz);
  }

  //----------------------------------------------------------------------
  TH3* ToTH3Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx,
		   const Binning& binsy,
		   const Binning& binsz,
		   ana::EBinType bintype)
  {

    const int nx = binsx.NBins();
    const int ny = binsy.NBins();
    const int nz = binsz.NBins();

    // Make sure it's compatible with having been made with this binning
    assert(h1->GetNbinsX() == nx*ny*nz);

    TH3* ret;

    // If all three axes are simple, we can call a simpler constructor
    if(binsx.IsSimple() && binsy.IsSimple() && binsz.IsSimple()){
      ret = new TH3F(UniqueName().c_str(), "",
                     nx, binsx.Min(), binsx.Max(),
                     ny, binsy.Min(), binsy.Max(),
                     nz, binsz.Min(), binsz.Max());

      if(!binsx.IsSimple() || !binsy.IsSimple() || !binsz.IsSimple()){
        // TH3 doesn't have the constructors for mixed simple and custom
        std::cerr << "ToTH3: one or more axes is custom, but not all three. Applying Simple binning to all three axes" << std::endl;
      }
    }
    else{
      ret = new TH3F(UniqueName().c_str(), "",
                     nx, &binsx.Edges().front(),
                     ny, &binsy.Edges().front(),
                     nz, &binsz.Edges().front());
    }

    for(int i = 0; i < h1->GetNbinsX(); ++i){
      const double val = h1->GetBinContent(i+1);
      const double err = h1->GetBinError(i+1);

      const int nynz = ny*nz;
      const int nmodnynz = i%nynz;
      const int ix = i/nynz;
      const int iy = nmodnynz/nz;
      const int iz = i%nz;

      ret->SetBinContent(ix+1, iy+1, iz+1, val);
      ret->SetBinError  (ix+1, iy+1, iz+1, err);
    }

    if(bintype == ana::EBinType::kBinDensity) ret->Scale(1, "width");

    return ret;

  }

  //----------------------------------------------------------------------
  std::vector<std::string> Wildcard(const std::string& wildcardString)
  {
    // Expand environment variables and wildcards like the shell would
    wordexp_t p;
    const int status = wordexp(wildcardString.c_str(), &p, WRDE_SHOWERR);

    if(status != 0){
      std::cerr << "Wildcard string '" << wildcardString
                << "' returned error " << status << " from wordexp()."
                << std::endl;
      return {};
    }

    std::vector<std::string> fileList;

    for(unsigned int i = 0; i < p.we_wordc; ++i){
      // Check the file exists before adding it
      struct stat sb;
      if(stat(p.we_wordv[i], &sb) == 0)
        fileList.push_back(p.we_wordv[i]);
    }

    wordfree(&p);

    return fileList;
  }

  //----------------------------------------------------------------------
  std::string FindCAFAnaDir()
  {
    const char* pub = getenv("SRT_PUBLIC_CONTEXT");
    const char* priv = getenv("SRT_PRIVATE_CONTEXT");

    if(priv && priv != std::string(".")){
      const std::string ret = std::string(priv)+"/CAFAna/";
      struct stat junk;
      if(stat(ret.c_str(), &junk) == 0) return ret;
    }
    assert(pub);
    return std::string(pub)+"/CAFAna/";
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
    return (getenv("_CONDOR_SCRATCH_DIR") != 0);
  }


  //----------------------------------------------------------------------
  bool AlmostEqual(double a, double b)
  {
    if(a == 0 && b == 0) return true;

    return fabs(a-b)/std::max(a, b) < .0001; // allow 0.01% error
  }


  //----------------------------------------------------------------------
  std::string pnfs2xrootd(std::string loc, bool unauth)
  {
    static bool first = true;
    static bool onsite = false;

    if (first && unauth) {
      first = false;
      char chostname[255];
      gethostname(chostname, 255);
      std::string hostname = chostname;

      if ( hostname.find("fnal.gov") != std::string::npos ){
        onsite = true;
        std::cout << "Using unauthenticated xrootd access (port 1095) while on-site, hostname: " << hostname << std::endl;
      }
      else {
        onsite = false;
        std::cout << "Using authenticated xrootd access (port 1094) access while off-site, hostname: " << hostname << std::endl;
      }
    }

    if(loc.rfind("/pnfs/", 0) == 0){ // ie begins with
      if ( onsite && unauth )
        loc = std::string("root://fndcagpvm01.fnal.gov:1095//pnfs/fnal.gov/usr/")+&loc.c_str()[6];
      else
        loc = std::string("root://fndcagpvm01.fnal.gov:1094//pnfs/fnal.gov/usr/")+&loc.c_str()[6];
    }
    return loc;
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
    double s[fNOsc] = {0};
    double c[fNOsc] = {0};
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
  TH1* GetMaskHist(const Spectrum& s, double xmin, double xmax, double ymin, double ymax)
  {
    if (s.GetBinnings().size() > 2){
      std::cout << "Error: unable to apply a mask in " << s.GetBinnings().size() << " dimensions" << std::endl;
      abort();
    }

    // The exposure isn't important here
    TH1* fMaskND  = s.ToTHX(s.POT());
    TH1D* fMask1D = s.ToTH1(s.POT());

    int ybins = fMaskND->GetNbinsY();

    for(int i = 0; i < fMask1D->GetNbinsX()+2; ++i){

      int ix = i / ybins;
      int iy = i % ybins;

      bool isMask = false;

      if (xmin < xmax){
	if (fMaskND->GetXaxis()->GetBinLowEdge(ix+1) < xmin) isMask=true;
	if (fMaskND->GetXaxis()->GetBinUpEdge(ix+1) > xmax) isMask=true;
      }

      if (ymin < ymax){
	if (fMaskND->GetYaxis()->GetBinLowEdge(iy+1) < ymin) isMask=true;
	if (fMaskND->GetYaxis()->GetBinUpEdge(iy+1) > ymax) isMask=true;
      }

      if (isMask) fMask1D->SetBinContent(i+1, 0);
      else fMask1D->SetBinContent(i+1, 1);

    }
    return fMask1D;
  }
}
