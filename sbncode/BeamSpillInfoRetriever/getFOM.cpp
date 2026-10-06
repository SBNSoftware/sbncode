/**
 * @file sbncode/BeamSpillInfoRetriever/getFOM.cpp
 * @brief Figure of Merit for BNB Spills using SBND information adapted from MicroBooNE FOM
 * @author Max Dubnowski (maxdub@sas.upenn.edu or `@Max Dubnowski` on SBN Slack)
 */
#include "sbncode/BeamSpillInfoRetriever/getFOM.h"
#include <math.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <tuple>
#include "TH1D.h"
#include "TFitResult.h"
#include <vector>


using namespace std;


namespace sbn {
  
  namespace {

    /// Value used by the retrievers for a device that could not be read.
    constexpr double MissingValue = -999.;

    /// A device reading is usable if it is finite and not the "missing" marker.
    bool isValid(double v) { return std::isfinite(v) && v != MissingValue; }

    /// Number of wires per plane in the multiwire chambers.
    constexpr std::size_t NWires = 48;

    /// Straight-line extrapolation to the target centre through two BPMs.
    /// Positions in mm, z in m: the slope is in mm/m, i.e. mrad, which is what
    /// calcFOM() expects for the angle (no atan: the slope already *is* the angle).
    std::pair<double, double> extrapolate
      (double delta0, double z0, double delta1, double z1, double ztarget)
    {
      double const ang = (delta1 - delta0) / (z1 - z0);
      double const pos = delta0 + ang * (ztarget - z0);
      return { ang, pos };
    }

  } // local namespace


  namespace {

    // Z positions of the monitors [m]
    constexpr double vp873_zpos = 191.153656;
    constexpr double hp875_zpos = 202.116104;
    constexpr double vp875_zpos = 202.3193205;
    constexpr double hptg1_zpos = 204.833267;
    constexpr double hptg2_zpos = 205.240662;
    constexpr double vptg2_zpos = 205.036835;
    constexpr double target_center_zpos = 206.870895;

    // multiwire sigma -> sigma at the target (quadratic polynomials)
    constexpr double p875x[] = {0.431857, 0.158077, 0.00303551};
    constexpr double p875y[] = {0.279128, 0.337048, 0};
    constexpr double p876x[] = {0.166172, 0.30999, -0.00630299};
    constexpr double p876y[] = {0.13425, 0.580862, 0};

    constexpr double smallSigmaX = 0.5, largeSigmaX = 10, smallSigmaY = 0.3, largeSigmaY = 10;
    constexpr double maxChi2X = 20, maxChi2Y = 20;
    bool inWindow(double sx, double sy)
      { return sx > smallSigmaX && sx < largeSigmaX && sy > smallSigmaY && sy < largeSigmaY; }

    /// A database (online-fitted) multiwire sigma above this [mm] means the
    /// chamber saw no beam: the online fit went to noise (beam: 1-2.5 mm, empty
    /// profile: ~7 mm). Such a width is not used.
    constexpr double EmptyMWSigma = 4.0;

    double poly(double const* p, double s) { return p[0] + p[1]*s + p[2]*s*s; }

    /**
     * Autotune calibration times [s, UTC] of the two horizontal target BPMs
     * (from IFBeam, Z. Pavlovic; same table as sbnana's getBNBFoM.cxx).
     * Autotune steers the beam on HP875 and on the target BPM calibrated most
     * recently, so that one reads the beam at its setpoint; the offset of the
     * other one may be stale (HPTG1 before 2023-03-03: +1.5 mm, i.e. the beam
     * projected 2.6 mm off-centre). New autotune calibrations must be added here.
     */
    constexpr unsigned long HPTG1CalibTimes[] = { 1420092000, 1576014670, 1606510357, 1677887110, 1711030691 };
    constexpr unsigned long HPTG2CalibTimes[] = { 1420092000, 1574190444, 1576014670, 1588603117, 1605706657,
                                                  1606510357, 1608593857, 1668802510 };

    template <std::size_t N>
    unsigned long lastCalibration(unsigned long const (&times)[N], unsigned long t) {
      unsigned long last = 0;
      for (unsigned long c: times) if (c <= t) last = c;
      return last;
    }

    /// Fills position, angle and the related status bits of `state`.
    void fillGeometry(BNBSpillInfo const& spill, BNBBeamState& state) {
      using namespace fomstatus;
      // ---- horizontal: HP875 and the primary target BPM, offsets subtracted.
      // The other target BPM is not used as a fallback: its offset is not the
      // one autotune steers on, so the position would be biased (mm level).
      bool const useHPTG2 = hptg2IsPrimary(spill.spill_time_s);
      if (useHPTG2) state.status |= HPTG2Primary;
      bool const okHP875 = isValid(spill.HP875) && isValid(spill.HP875Offset);
      bool const okHPTG1 = isValid(spill.HPTG1) && isValid(spill.HPTG1Offset);
      bool const okHPTG2 = isValid(spill.HPTG2) && isValid(spill.HPTG2Offset);
      bool const okPrimary = useHPTG2? okHPTG2: okHPTG1;
      if (!okHP875 || !okPrimary) {
        state.status |= NoHBPM;
        if (okHP875 && (useHPTG2? okHPTG1: okHPTG2)) state.status |= FallbackRemoved;
      }
      else {
        double const delta_hp875 = spill.HP875 - spill.HP875Offset;
        std::tie(state.hang, state.hpos) = useHPTG2
          ? extrapolate(delta_hp875, hp875_zpos, spill.HPTG2 - spill.HPTG2Offset, hptg2_zpos, target_center_zpos)
          : extrapolate(delta_hp875, hp875_zpos, spill.HPTG1 - spill.HPTG1Offset, hptg1_zpos, target_center_zpos);
      }

      // ---- vertical: VP875 and VP873, offsets subtracted (VPTG2 is not used
      // as a fallback for the same reason).
      bool const okVP875 = isValid(spill.VP875) && isValid(spill.VP875Offset);
      bool const okVP873 = isValid(spill.VP873) && isValid(spill.VP873Offset);
      bool const okVPTG2 = isValid(spill.VPTG2) && isValid(spill.VPTG2Offset);
      if (!okVP875 || !okVP873) {
        state.status |= NoVBPM;
        if (okVP875 && okVPTG2) state.status |= FallbackRemoved;
      }
      else {
        double const delta_vp875 = spill.VP875 - spill.VP875Offset;
        std::tie(state.vang, state.vpos) =
          extrapolate(delta_vp875, vp875_zpos, spill.VP873 - spill.VP873Offset, vp873_zpos, target_center_zpos);
      }
    }

    /// Width from the multiwire profiles fitted here (M876, then M875).
    /// The target multiwire (MMBTBB) is not used: its profiles are noise.
    bool fittedWidth(BNBSpillInfo const& spill, double& tgtsx, double& tgtsy, int& source) {
      struct MWDevice_t { std::vector<int> const* data; double const* px; double const* py; int source; };
      MWDevice_t const mwDevices[] = {
        { &spill.M876BB, p876x, p876y, BNBBeamState::FitM876 },
        { &spill.M875BB, p875x, p875y, BNBBeamState::FitM875 },
      };
      for (auto const& dev: mwDevices) {
        // each profile is 48 horizontal wires followed by 48 vertical ones
        if (dev.data->size() < 2*NWires) continue;
        std::vector<double> const mw(dev.data->begin(), dev.data->begin() + 2*NWires);
        double xx, yy, sx, sy, chi2x, chi2y;
        bool const fitOK = processBNBprofile(&mw[0], xx, sx, chi2x)
                         & processBNBprofile(&mw[NWires], yy, sy, chi2y);
        if (!fitOK) continue;
        sx = poly(dev.px, sx);
        sy = poly(dev.py, sy);
        if (inWindow(sx, sy) && chi2x < maxChi2X && chi2y < maxChi2Y) {
          tgtsx = sx; tgtsy = sy; source = dev.source;
          return true;
        }
      }
      return false;
    }

    /// Width from the online ("pre-fit", database) multiwire widths.
    /// Sets `empty` if the device chosen saw no beam (no width returned then).
    bool databaseWidth(BNBSpillInfo const& spill, double& tgtsx, double& tgtsy, int& source, bool& empty) {
      struct DBWidth_t { double hs, vs; double const* px; double const* py; int source; };
      DBWidth_t const dbWidths[] = {
        { spill.M876HS, spill.M876VS, p876x, p876y, BNBBeamState::DatabaseM876 },
        { spill.M875HS, spill.M875VS, p875x, p875y, BNBBeamState::DatabaseM875 },
      };
      empty = false;
      for (auto const& w: dbWidths) {
        if (!isValid(w.hs) || !isValid(w.vs)) continue;
        double const sx = poly(w.px, w.hs);
        double const sy = poly(w.py, w.vs);
        if (!inWindow(sx, sy)) continue;
        if (w.hs > EmptyMWSigma || w.vs > EmptyMWSigma) { empty = true; return false; }
        tgtsx = sx; tgtsy = sy; source = w.source;
        return true;
      }
      return false;
    }

  } // local namespace


  bool hptg2IsPrimary(unsigned long spill_time_s) {
    return lastCalibration(HPTG2CalibTimes, spill_time_s) >= lastCalibration(HPTG1CalibTimes, spill_time_s);
  }


  double beamIntensity(BNBSpillInfo const& spill, bool* tor875used) {
    // TOR860 sometimes reads ~1e-3 of the beam TOR875 sees (TOR860 drop-out):
    // then TOR875 is the intensity. A missing toroid is -999.
    bool const ok860 = isValid(spill.TOR860) && spill.TOR860 > 0.;
    bool const ok875 = isValid(spill.TOR875) && spill.TOR875 > 0.;
    bool const dropout = ok860 && ok875 && spill.TOR860 < TOR860DropoutFraction * spill.TOR875;
    if (tor875used) *tor875used = dropout;
    if (ok860 && !dropout) return spill.TOR860;
    if (ok875) return spill.TOR875;
    return -1.;
  }


  BNBBeamState getBNBBeamState(BNBSpillInfo const& spill) {
    using namespace fomstatus;
    BNBBeamState state;
    bool tor875used = false;
    state.tor = beamIntensity(spill, &tor875used);
    if (tor875used) state.status |= TOR875Used;
    if (state.tor <= 0.) state.status |= NoTOR;
    fillGeometry(spill, state);

    // width: fitted profile, else database width, else nominal
    bool empty = false;
    if (!fittedWidth(spill, state.sx, state.sy, state.widthSource)
      && !databaseWidth(spill, state.sx, state.sy, state.widthSource, empty))
    {
      state.sx = state.sy = MissingValue;
      state.widthSource = BNBBeamState::Nominal;
      state.status |= NoMWWidth;
      if (empty) state.status |= MWEmpty;
    }
    return state;
  }


  double computeFOM(BNBBeamState const& state, bool nominalWidth) {
    if (!state.hasFOM()) return MissingValue;
    bool const measured = !nominalWidth && state.widthSource != BNBBeamState::Nominal;
    return measured
      ? 1-pow(10, sbn::calcFOM(state.hpos, state.hang, state.vpos, state.vang, state.tor, state.sx, state.sy))
      : 1-pow(10, sbn::calcFOM(state.hpos, state.hang, state.vpos, state.vang, state.tor));
  }


  std::tuple<float, float, float> getBNBqualityFOM(BNBSpillInfo const& spill)
  {
    BNBBeamState const state = getBNBBeamState(spill);
    if (state.status & fomstatus::NoTOR)  return {-1, -1, -1};
    if (state.status & fomstatus::NoHBPM) return {2, 2, 2};
    if (state.status & fomstatus::NoVBPM) return {3, 3, 3};

    // FOM with the width fitted from the multiwire profiles
    double fom = MissingValue, prefitfom = MissingValue;
    double sx = MissingValue, sy = MissingValue;
    int source = BNBBeamState::Nominal;
    if (fittedWidth(spill, sx, sy, source))
      fom = 1-pow(10, sbn::calcFOM(state.hpos, state.hang, state.vpos, state.vang, state.tor, sx, sy));

    // "pre-fit" FOM with the widths fitted online (database M876/M875)
    bool empty = false;
    if (databaseWidth(spill, sx, sy, source, empty))
      prefitfom = 1-pow(10, sbn::calcFOM(state.hpos, state.hang, state.vpos, state.vang, state.tor, sx, sy));

    // FOM with the nominal beam width (scale factors 1)
    double const noMWfom = computeFOM(state, true);
    return {fom, prefitfom, noMWfom};
  }


/**
    * @brief Extracts statistics from multiwire monitor data.
    * @param mwdata pointer to multiwire data (48 channels horizontal or vertical)
    * @param[out] x mean position from the fit [mm]
    * @param[out] sx &sigma; from the fit [mm]
    * @param[out] chi2 &chi;&sup2;/NDF for the Gaussian fit
    * @return whether the fit succeeded
    *
    * This function takes multiwire data (`mwdata`),
    * finds the min and max,
    * finds the first and last bin where amplitude is greater than 20%,
    * fits the peak between first and last bin with Gaussian (assuming 2% relative errors)
    * and returns the parameters of the fit.
    * The 48 wires cover 24 mm (0.5 mm pitch).
    */
  bool processBNBprofile(const double* mwdata, double &x, double& sx, double& chi2)
  {
    x = sx = chi2 = 99999;
    // values' sign is inverted
    double minx = std::min(-*std::max_element(mwdata, mwdata + NWires), 0.0);
    double maxx = std::max(-*std::min_element(mwdata, mwdata + NWires), 0.0);
    int first_x = -1; int last_x = -1;
    // local histogram, not registered in gDirectory (no name clashes, no leak)
    TH1D hProf("hProfMW","",NWires,-12.0,12.0);
    hProf.SetDirectory(nullptr);
    double const threshold = (maxx-minx)*0.2; // 20% of the range
    double const error = (maxx-minx)*0.02; // 2% of the range
    for (unsigned int i=0;i<NWires;i++) {
      hProf.SetBinContent(i+1,-mwdata[i]-minx);
      if (-mwdata[i]-minx    > threshold && first_x==-1) first_x=i;
      if (-mwdata[i]-minx    > threshold)                last_x=i+1;
      hProf.SetBinError(i+1,error);
    }
    if (hProf.GetSumOfWeights() <= 0 || first_x < 0) return false;

    TFitResultPtr const fit = hProf.Fit("gaus","QNS","",-12+first_x*0.5,-12+last_x*0.5);
    if (!fit.Get() || fit->Status() != 0 || fit->Ndf() <= 0) return false;
    x   = fit->Parameter(1);
    sx  = fit->Parameter(2);
    chi2= fit->Chi2() / fit->Ndf();
    return true;
  }


  double calcFOM(double horpos, double horang, double verpos, double verang, double ppp, double tgtsx, double tgtsy)
  {
    ppp /= 1e12; //converts to 10^12 POT
    

    //code from MiniBooNE AnalysisFramework with the addition of scaling the beam profile to match tgtsx, tgtsy
    //form DQ_BeamLine_twiss_init.F
    double bx  =  4.68;
    double ax  =  0.0389;
    double gx  = (1+ax*ax)/bx;
    double nx  =  0.0958;
    double npx = -0.0286;
    double by  = 59.12;
    double ay  =  2.4159;
    double gy  = (1+ay*ay)/by;
    double ny  =  0.4577;
    double npy = -0.0271;
    //from DQ_BeamLine_make_tgt_fom2.F
    double ex = 0.1775E-06 + 0.1827E-07*ppp;
    double ey = 0.1382E-06 + 0.2608E-08*ppp;
    double dp = 0.4485E-03 + 0.6100E-04*ppp;
    double tex = ex;
    double tey = ey;
    double tdp = dp;
    //from DQ_BeamLine_beam_init.F
    double sigma1[6][6]={{0}};
    double centroid1[6]={0};
    centroid1[0] = horpos;
    centroid1[1] = horang;
    centroid1[2] = verpos;
    centroid1[3] = verang;
    centroid1[4] = 0.0;
    centroid1[5] = 0.0;
    sigma1[5][5] =  tdp*tdp;
    sigma1[0][0] =  tex*bx+ nx*nx *tdp*tdp;
    sigma1[0][1] = -tex*ax+ nx*npx*tdp*tdp;
    sigma1[1][1] =  tex*gx+npx*npx*tdp*tdp;
    sigma1[0][5] =  nx*tdp*tdp;
    sigma1[1][5] =  npx*tdp*tdp;
    sigma1[1][0] =  sigma1[0][1];
    sigma1[5][0] =  sigma1[0][5];
    sigma1[5][1] =  sigma1[1][5];
    sigma1[2][2] =  tey*by+ny*ny*tdp*tdp;
    sigma1[2][3] = -tey*ay+ny*npy*tdp*tdp;
    sigma1[3][3] =  tey*gy+npy*npy*tdp*tdp;
    sigma1[2][5] =  ny*tdp*tdp;
    sigma1[3][5] =  npy*tdp*tdp;
    sigma1[3][2] =  sigma1[2][3];
    sigma1[5][2] =  sigma1[2][5];
    sigma1[5][3] =  sigma1[3][5];

    
    double begtocnt[6][6]={
      { 0.65954,  0.43311,  0.00321,  0.10786, 0.00000,  1.97230},
      { 0.13047,  1.60192,  0.00034,  0.00512, 0.00000,  1.96723},
      {-0.00287, -0.03677, -0.35277, -4.68056, 0.00000,  0.68525},
      {-0.00089, -0.00430, -0.17722, -5.18616, 0.00000,  0.32300},
      {-0.00104,  0.00232, -0.00001, -0.00224, 1.00000, -0.00450},
      { 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  1.00000}
    };
    double cnttoups[6][6]={{0}};
    double cnttodns[6][6]={{0}};
    double identity[6][6]={{0}};
    double begtoups[6][6]={{0}};
    double begtodns[6][6]={{0}};
    for (int i=0;i<6;i++) {
      for (int j=0;j<6;j++) {
	if (i==j) {
	  cnttoups[i][j] = 1.0;
	  cnttodns[i][j] = 1.0;
	  identity[i][j] = 1.0;
	} else {
	  cnttoups[i][j] = 0.0;
	  cnttodns[i][j] = 0.0;
	  identity[i][j] = 0.0;
	}
      }
    }
    cnttoups[0][1] = -0.35710;
    cnttoups[2][3] = -0.35710;
    cnttodns[0][1] = +0.35710;
    cnttodns[2][3] = +0.35710;
    for (int i=0;i<6;i++) {
      for (int j=0;j<6;j++) {
	for (int k=0;k<6;k++) {
	  begtoups[i][k] = begtoups[i][k] + cnttoups[i][j]*begtocnt[j][k];
	  begtodns[i][k] = begtodns[i][k] + cnttodns[i][j]*begtocnt[j][k];
	}
      }
    }
    //swim to upstream of target
    double cx, cy, sx, sy, rho;
    sbn::swimBNB(centroid1,sigma1,
		 cnttoups, begtoups,
		 cx, cy, sx, sy, rho);

    double scalex, scaley;
    if( tgtsx !=-999 && tgtsy!=-999){
        scalex=tgtsx/sx;
        scaley=tgtsy/sy;
    }
    else{
        scalex=1;
        scaley=1;
    }
    double fom_a=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    //swim to center of target
    sbn::swimBNB(centroid1,sigma1,
		 identity, begtocnt,
		 cx, cy, sx, sy, rho);
    double fom_b=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    //swim to downstream of target
    sbn::swimBNB(centroid1,sigma1,
		 cnttodns, begtodns,
		 cx, cy, sx, sy, rho);
    double fom_c=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    // add a guard for double precision

    if(fom_a <= -10000. || fom_b <= -10000. || fom_c <= -10000) return -10000.;
    double fom2=fom_a*0.6347 +
      fom_b*0.2812 +
      fom_c*0.0841;
    return fom2;
  }
  
  void swimBNB(const double centroid1[6], const double sigma1[6][6],
	       const double xferc[6][6], const double xfers[6][6],
	       double &cx, double& cy, double &sx, double &sy, double &rho)
  {
    //centroid
    double centroid2[6]={0};
    for (int i=0;i<6;i++) {
      for (int j=0;j<6;j++) {
	centroid2[i] = centroid2[i] + xferc[i][j]*centroid1[j];
      }
    }
    cx = centroid2[0];
    cy = centroid2[2];
    //sigma
    double sigma2[6][6]={{0}};
    for (int i = 0; i<6;i++) {
      for (int j = 0;j<6;j++) {
	for (int k = 0;k<6;k++) {
	  for (int m = 0;m<6;m++) {
	    sigma2[i][m] = sigma2[i][m] + xfers[i][j]*sigma1[j][k]*xfers[m][k];
	  }
	}
      }
    }
    //get beam sigma
    sx  = sqrt(sigma2[0][0])*1000.0;
    sy  = sqrt(sigma2[2][2])*1000.0;
    rho = sigma2[0][2]/sqrt(sigma2[0][0]*sigma2[2][2]);
  }
  
  
  
  double func_intbivar(const double cx, const double cy, const double sx, const double sy, const double rho )
  {
    //integrate beam overlap with target cylinder
    double x0  =  cx;
    double y0  =  cy;
    double dbin = 0.1;
    double dx = dbin;
    double dy = dbin;
    double r    = 4.75;
    double rr   = r*r;
    double rho2 = rho*rho;
    double xmin = -r;
    double ymin = -r;
    int imax = round((2.0*r)/dx);
    int jmax = round((2.0*r)/dy);
    double sum =  0.0;
    double x = xmin;
    for (int i=0;i<=imax;i++) {
      double y = ymin;
      for (int j=0;j<=jmax;j++) {
	if ( (x*x+y*y)<rr ) {
	  //Looking up Bivariate Normal Distribution, tx and ty should be x-x0 not x+x0; doesn't change much for a symmetric spill
	  double tx = (x-x0)/sx;
	  double ty = (y-y0)/sy;
	  //Code used in the MicroBooNE code, possible error in the formula
	  //double tx = (x+x0)/sx;
	  //double ty = (y+y0)/sy;
	  double z = tx*tx - 2.0*rho*tx*ty + ty*ty;
	  double t = exp(-z/(2.0*(1.0-rho2)));
	  sum = sum + t;
	}
	y = y + dy;
      }
      x = x + dx;
    }
    sum = sum*dx*dy/(2.0*M_PI*sx*sy*sqrt(1.0-rho2));


    // add a guard for double precision
    if(sum >= 1.) return -10000.;
    return log10(1-sum);
  }
  
}
