// \file  TrackMomentumCalculator.cxx
//
// \author sowjanyag@phys.ksu.edu

#include "TrackMomentumCalculator.h"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Math/Functor.h"

#include <array>
#include <cassert>

using std::cout;
using std::endl;

namespace {

  constexpr auto range_gramper_cm()
  {
    std::array<float, 29> Range_grampercm{{
      9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1,
      1.063E2,  1.725E2, 2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3,
      2.297E3,  4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
      4.326E4,  5.768E4, 7.734E4, 1.060E5, 1.307E5}};
    for (float& value : Range_grampercm) {
      value /= 1.396; // convert to cm
    }
    return Range_grampercm;
  }

  constexpr auto Range_grampercm = range_gramper_cm();
  constexpr std::array<float, 29> KE_MeV{{
    10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
    400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
    20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000}};
  TGraph const KEvsR{29, Range_grampercm.data(), KE_MeV.data()};
  TSpline3 const KEvsR_spline3{"KEvsRS", &KEvsR};


  TVector3 const basex{1, 0, 0};
  TVector3 const basey{0, 1, 0};
  TVector3 const basez{0, 0, 1};
  constexpr float kcal{0.0024};

  class FcnWrapper {
  public:
    explicit FcnWrapper(std::vector<double>&& xmeas,
                        std::vector<double>&& ymeas,
                        std::vector<double>&& eymeas)
      : xmeas_{xmeas}
      , ymeas_{ymeas}
      , eymeas_{eymeas}
    {}

    double
    my_mcs_chi2(double const* x) const
    {
      double result = 0.0;

      double p = x[0];
      double theta0 = x[1];

      auto const n = xmeas_.size();
      assert(n == ymeas_.size());
      assert(n == eymeas_.size());

      for (std::size_t i = 0; i < n; ++i) {
        double const xx = xmeas_[i];
        double const yy = ymeas_[i];
        double const ey = eymeas_[i];

        if (std::abs(ey) < std::numeric_limits<double>::epsilon()) {
          std::cout << " Zero denominator in my_mcs_chi2 ! " << std::endl;
          return -1;
        }

        constexpr double rad_length{14.0};
        double const l0 = xx / rad_length;
        double res1 = 0.0;

        if (xx > 0 && p > 0)
          res1 = (13.6 / p) * std::sqrt(l0) * (1.0 + 0.038 * std::log(l0));

        res1 = std::sqrt(res1 * res1 + theta0 * theta0);

        double const diff = yy - res1;
        result += cet::square(diff / ey);
      }

      result += 2.0 / (4.6) * theta0; // *std::log( 1.0/14.0 );

      if (std::isnan(result) || std::isinf(result)) {
        MF_LOG_DEBUG("TrackMomentumCalculator") << " Is nan in my_mcs_chi2 ! ";
        return -1;
      }

      return result;
    }

  private:
    std::vector<double> const xmeas_;
    std::vector<double> const ymeas_;
    std::vector<double> const eymeas_;
  };

}

namespace trkf {

  TrackMomentumCalculator::TrackMomentumCalculator(double const min,
                                                   double const max)
    : minLength{min}
    , maxLength{max}
  {
    for (int i = 1; i <= n_steps; i++) {
      steps.push_back(steps_size * i);
    }
  }

  double
  TrackMomentumCalculator::GetTrackMomentum(double trkrange, int pdg) const
  {
    /* Muon range-momentum tables from CSDA (Argon density = 1.4 g/cm^3)
       website:
       http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf

       CSDA table values:
       float Range_grampercm[30] = {9.833E-1, 1.786E0, 3.321E0,
       6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
       2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
       4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
       4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5}; float KE_MeV[30] = {10, 14,
       20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000,
       4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000,
       200000, 300000, 400000};

       Functions below are obtained by fitting polynomial fits to KE_MeV vs
       Range (cm) graph. A better fit was obtained by splitting the graph into
       two: Below range<=200cm,a polynomial of power 4 was a good fit; above
       200cm, a polynomial of power 6 was a good fit

       Fit errors for future purposes:
       Below 200cm, Forpoly4 fit: p0 err=1.38533;p1 err=0.209626; p2
       err=0.00650077; p3 err=6.42207E-5; p4 err=1.94893E-7; Above 200cm,
       Forpoly6 fit: p0 err=5.24743;p1 err=0.0176229; p2 err=1.6263E-5; p3
       err=5.9155E-9; p4 err=9.71709E-13; p5 err=7.22381E-17;p6
       err=1.9709E-21;*/

    ///////////////////////////////////////////////////////////////////////////
    //*********For muon, the calculations are valid up to 1.91E4 cm range
    //corresponding to a Muon KE of 40 GeV**********//
    ///////////////////////////////////////////////////////////////////////////

    /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
      website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html

      CSDA values:
      double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
      350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
      1500, 2000, 2500, 3000, 4000, 5000};

      double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
      2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
      9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
      2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
      7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};

      Functions below are obtained by fitting power and polynomial fits to
      KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
      graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
      polynomial of power 6 was a good fit

      Fit errors for future purposes:
      For power function fit: a=0.388873; and b=0.00347075
      Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
      err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
      err=2.96958E-17;*/

    ///////////////////////////////////////////////////////////////////////////
    //*********For proton, the calculations are valid up to 3.022E3 cm range
    //corresponding to a Muon KE of 5 GeV**********//
    ///////////////////////////////////////////////////////////////////////////

    if (trkrange < 0 || std::isnan(trkrange)) {
      mf::LogError("TrackMomentumCalculator")
        << "Invalid track range " << trkrange << " return -1";
      return -1.;
    }

    double KE, Momentum, M;
    constexpr double Muon_M = 105.7, Proton_M = 938.272;

    if (abs(pdg) == 13) {
      M = Muon_M;
      KE = KEvsR_spline3.Eval(trkrange);
    } else if (abs(pdg) == 2212) {
      M = Proton_M;
      if (trkrange > 0 && trkrange <= 80)
        KE = 29.9317 * std::pow(trkrange, 0.586304);
      else if (trkrange > 80 && trkrange <= 3.022E3)
        KE =
          149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
          (4.34587E-6 * trkrange * trkrange * trkrange) +
          (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
          (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
          (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *
           trkrange);
      else
        KE = -999;
    } else
      KE = -999;

    if (KE < 0)
      Momentum = -999;
    else
      Momentum = std::sqrt((KE * KE) + (2 * M * KE));

    Momentum = Momentum / 1000;

    return Momentum;
  }

  // Momentum measurement via Multiple Coulomb Scattering (MCS)

  // author: Leonidas N. Kalousis (July 2015)

  // email: kalousis@vt.edu

  double
  TrackMomentumCalculator::GetMomentumMultiScatterLLHD(
    const art::Ptr<recob::Track>& trk)
  {
    std::vector<float> recoX;
    std::vector<float> recoY;
    std::vector<float> recoZ;

    int n_points = trk->NumberTrajectoryPoints();

    for (int i = 0; i < n_points; i++) {
      auto const& pos = trk->LocationAtPoint(i);
      recoX.push_back(pos.X());
      recoY.push_back(pos.Y());
      recoZ.push_back(pos.Z());
    }

    if (recoX.size() < 2)
      return -1.0;

    if (!plotRecoTracks_(recoX, recoY, recoZ))
      return -1.0;

    constexpr double seg_size{10.};

    auto const segments = getSegTracks_(recoX, recoY, recoZ, seg_size);
    if (!segments.has_value())
      return -1.0;

    auto const seg_steps = segments->x.size();
    if (seg_steps < 2)
      return -1;

    double const recoL = segments->L.at(seg_steps - 1);
    if (recoL < minLength || recoL > maxLength)
      return -1;

    std::vector<float> dEi;
    std::vector<float> dEj;
    std::vector<float> dthij;
    std::vector<float> ind;
    if (getDeltaThetaij_(dEi, dEj, dthij, ind, *segments, seg_size) != 0)
      return -1.0;

    double logL = 1e+16;
    double bf = -666.0; // double errs = -666.0;

    int const start1{};
    int const end1{750};
    int const start2{};
    int const end2{}; // 800.0;

    for (int k = start1; k <= end1; ++k) {
      double const p_test = 0.001 + k * 0.01;

      for (int l = start2; l <= end2; ++l) {
        double const res_test = 2.0; // 0.001+l*1.0;
        double const fv = my_mcs_llhd(dEi, dEj, dthij, ind, p_test, res_test);

        if (fv < logL) {
          bf = p_test;
          logL = fv;
          // errs = res_test;
        }
      }
    }
    return bf;
  }

  TVector3
  TrackMomentumCalculator::GetMultiScatterStartingPoint(
    const art::Ptr<recob::Track>& trk)
  {
    double const LLHDp = GetMuMultiScatterLLHD3(trk, true);
    double const LLHDm = GetMuMultiScatterLLHD3(trk, false);

    if (LLHDp != -1 && LLHDm != -1 && LLHDp > LLHDm) {
      int const n_points = trk->NumberTrajectoryPoints();
      return trk->LocationAtPoint<TVector3>(n_points - 1);
    } else {
      return trk->LocationAtPoint<TVector3>(0);
    }

    // Should never get here
    return TVector3{};
  }

  double
  TrackMomentumCalculator::GetMuMultiScatterLLHD3(
    art::Ptr<recob::Track> const& trk,
    bool const dir)
  {
    std::vector<float> recoX;
    std::vector<float> recoY;
    std::vector<float> recoZ;

    int const n_points = trk->NumberTrajectoryPoints();
    for (int i = 0; i < n_points; ++i) {
      auto const index = dir ? i : n_points - 1 - i;
      auto const& pos = trk->LocationAtPoint(index);
      recoX.push_back(pos.X());
      recoY.push_back(pos.Y());
      recoZ.push_back(pos.Z());
    }

    if (recoX.size() < 2)
      return -1.0;

    if (!plotRecoTracks_(recoX, recoY, recoZ))
      return -1.0;

    constexpr double seg_size{5.0};
    auto const segments = getSegTracks_(recoX, recoY, recoZ, seg_size);
    if (!segments.has_value())
      return -1.0;

    auto const seg_steps = segments->x.size();
    if (seg_steps < 2)
      return -1;

    double const recoL = segments->L.at(seg_steps - 1);
    if (recoL < 15.0 || recoL > maxLength)
      return -1;

    std::vector<float> dEi;
    std::vector<float> dEj;
    std::vector<float> dthij;
    std::vector<float> ind;
    if (getDeltaThetaij_(dEi, dEj, dthij, ind, *segments, seg_size) != 0)
      return -1.0;

    double const p_range = recoL * kcal;
    double const logL = my_mcs_llhd(dEi, dEj, dthij, ind, p_range, 5.65);

    return logL;
  }

  int
  TrackMomentumCalculator::getDeltaThetaij_(std::vector<float>& ei,
                                            std::vector<float>& ej,
                                            std::vector<float>& th,
                                            std::vector<float>& ind,
                                            Segments const& segments,
                                            double const thick) const
  {
    int const a1 = segments.x.size();
    int const a2 = segments.y.size();
    int const a3 = segments.z.size();

    if (a1 != a2 || a1 != a3) {
      std::cout << " ( Get thij ) Error ! " << std::endl;
      return -1.0;
    }

    auto const& segnx = segments.nx;
    auto const& segny = segments.ny;
    auto const& segnz = segments.nz;
    auto const& segL = segments.L;

    int tot = a1 - 1;
    double thick1 = thick + 0.13;

    for (int i = 0; i < tot; i++) {
      double const dx = segnx.at(i);
      double const dy = segny.at(i);
      double const dz = segnz.at(i);

      TVector3 const vec_z{dx, dy, dz};
      TVector3 vec_x;
      TVector3 vec_y;

      double const switcher = basex.Dot(vec_z);
      if (std::abs(switcher) <= 0.995) {
        vec_y = vec_z.Cross(basex).Unit();
        vec_x = vec_y.Cross(vec_z);
      }
      else {
        // cout << " It switched ! Isn't this lovely !!! " << endl; //
        vec_y = basez.Cross(vec_z).Unit();
        vec_x = vec_y.Cross(vec_z);
      }

      TVector3 const Rx{vec_x.Dot(basex), vec_x.Dot(basey), vec_x.Dot(basez)};
      TVector3 const Ry{vec_y.Dot(basex), vec_y.Dot(basey), vec_y.Dot(basez)};
      TVector3 const Rz{vec_z.Dot(basex), vec_z.Dot(basey), vec_z.Dot(basez)};

      double const refL = segL.at(i);

      for (int j = i; j < tot; j++) {
        double const L1 = segL.at(j);
        double const L2 = segL.at(j + 1);

        double const dz1 = L1 - refL;
        double const dz2 = L2 - refL;

        if (dz1 <= thick1 && dz2 > thick1) {
          double const here_dx = segnx.at(j);
          double const here_dy = segny.at(j);
          double const here_dz = segnz.at(j);

          TVector3 const here_vec{here_dx, here_dy, here_dz};
          TVector3 const rot_here{
            Rx.Dot(here_vec), Ry.Dot(here_vec), Rz.Dot(here_vec)};

          double const scx = rot_here.X();
          double const scy = rot_here.Y();
          double const scz = rot_here.Z();

          double const azy = find_angle(scz, scy);
          double const azx = find_angle(scz, scx);

          constexpr double ULim = 10000.0;
          constexpr double LLim = -10000.0;

          double const cL = kcal;
          double const Li = segL.at(i);
          double const Lj = segL.at(j);

          if (azy <= ULim && azy >= LLim) {
            ei.push_back(Li * cL);
            ej.push_back(Lj * cL);
            th.push_back(azy);
            ind.push_back(2);
          }

          if (azx <= ULim && azx >= LLim) {
            ei.push_back(Li * cL);
            ej.push_back(Lj * cL);
            th.push_back(azx);
            ind.push_back(1);
          }

          break; // of course !
        }
      }
    }

    return 0;
  }

  double
  TrackMomentumCalculator::GetMomentumMultiScatterChi2(
    const art::Ptr<recob::Track>& trk)
  {
    std::vector<float> recoX;
    std::vector<float> recoY;
    std::vector<float> recoZ;

    int n_points = trk->NumberTrajectoryPoints();

    for (int i = 0; i < n_points; i++) {
      auto const& pos = trk->LocationAtPoint(i);
      recoX.push_back(pos.X());
      recoY.push_back(pos.Y());
      recoZ.push_back(pos.Z());
    }

    if (recoX.size() < 2)
      return -1.0;

    if (!plotRecoTracks_(recoX, recoY, recoZ))
      return -1.0;

    double const seg_size{steps_size};
    auto const segments = getSegTracks_(recoX, recoY, recoZ, seg_size);
    if (!segments.has_value())
      return -1.0;

    auto const seg_steps = segments->x.size();
    if (seg_steps < 2)
      return -1;

    double const recoL = segments->L.at(seg_steps - 1);
    if (recoL < minLength || recoL > maxLength)
      return -1;

    double ymax = -999.0;
    double ymin = +999.0;

    std::vector<double> xmeas;
    std::vector<double> ymeas;
    std::vector<double> eymeas;
    xmeas.reserve(n_steps);
    ymeas.reserve(n_steps);
    eymeas.reserve(n_steps);
    for (int j = 0; j < n_steps; j++) {
      double const trial = steps.at(j);
      auto const [mean, rms, rmse] = getDeltaThetaRMS_(*segments, trial);

      if (std::isnan(mean) || std::isinf(mean)) {
        mf::LogDebug("TrackMomentumCalculator") << "Returned mean is either nan or infinity.";
        continue;
      }
      if (std::isnan(rms) || std::isinf(rms)) {
        mf::LogDebug("TrackMomentumCalculator") << "Returned rms is either nan or infinity.";
        continue;
      }
      if (std::isnan(rmse) || std::isinf(rmse)) {
        mf::LogDebug("TrackMomentumCalculator") << "Returned rmse is either nan or infinity.";
        continue;
      }

      xmeas.push_back(trial); // Is this what is intended?
      ymeas.push_back(rms);
      eymeas.push_back(std::sqrt(cet::sum_of_squares(rmse, 0.05 * rms))); // <--- conservative syst. error to fix chi^{2} behaviour !!!

      if (ymin > rms)
        ymin = rms;
      if (ymax < rms)
        ymax = rms;
    }

    assert(xmeas.size() == ymeas.size());
    assert(xmeas.size() == eymeas.size());
    if (xmeas.empty()) {
      return -1.0;
    }

    TGraphErrors gr_meas{n_steps, xmeas.data(), ymeas.data(), nullptr, eymeas.data()};

    gr_meas.SetTitle("(#Delta#theta)_{rms} versus material thickness; Material "
                     "thickness in cm; (#Delta#theta)_{rms} in mrad");

    gr_meas.SetLineColor(kBlack);
    gr_meas.SetMarkerColor(kBlack);
    gr_meas.SetMarkerStyle(20);
    gr_meas.SetMarkerSize(1.2);

    gr_meas.GetXaxis()->SetLimits(steps.at(0) - steps.at(0),
                                  steps.at(n_steps - 1) + steps.at(0));
    gr_meas.SetMinimum(0.0);
    gr_meas.SetMaximum(1.80 * ymax);

    ROOT::Minuit2::Minuit2Minimizer mP{};
    FcnWrapper const wrapper{move(xmeas), move(ymeas), move(eymeas)};
    ROOT::Math::Functor FCA([&wrapper](double const* xs) { return wrapper.my_mcs_chi2(xs); }, 2);

    mP.SetFunction(FCA);
    mP.SetLimitedVariable(0, "p_{MCS}", 1.0, 0.01, 0.001, 7.5);
    mP.SetLimitedVariable(1, "#delta#theta", 0.0, 1.0, 0.0, 45.0);
    mP.SetMaxFunctionCalls(1.E9);
    mP.SetMaxIterations(1.E9);
    mP.SetTolerance(0.01);
    mP.SetStrategy(2);
    mP.SetErrorDef(1.0);

    bool const mstatus = mP.Minimize();

    mP.Hesse();

    const double* pars = mP.X();
    const double* erpars = mP.Errors();

    double const deltap = (recoL * kcal) / 2.0;

    double const p_mcs = pars[0] + deltap;
    double const p_mcs_e [[maybe_unused]] = erpars[0];
    return mstatus ? p_mcs : -1.0;
  }

  bool
  TrackMomentumCalculator::plotRecoTracks_(std::vector<float> const& xxx,
                                           std::vector<float> const& yyy,
                                           std::vector<float> const& zzz)
  {
    auto const n = xxx.size();
    auto const y_size = yyy.size();
    auto const z_size = zzz.size();

    if (n != y_size || n != z_size) {
      cout << " ( Get reco tacks ) Error ! " << endl;
      return false;
    }

    // Here, we perform a const-cast to float* because, sadly,
    // TPolyLine3D requires a pointer to a non-const object.  We will
    // trust that ROOT does not mess around with the underlying data.
    auto xs = const_cast<float*>(xxx.data());
    auto ys = const_cast<float*>(yyy.data());
    auto zs = const_cast<float*>(zzz.data());

    auto const narrowed_size =
      static_cast<int>(n); // ROOT requires a signed integral type
    delete gr_reco_xyz;
    gr_reco_xyz = new TPolyLine3D{narrowed_size, zs, xs, ys};
    gr_reco_yz = TGraph{narrowed_size, zzz.data(), yyy.data()};
    gr_reco_xz = TGraph{narrowed_size, zzz.data(), xxx.data()};
    gr_reco_xy = TGraph{narrowed_size, xxx.data(), yyy.data()};

    return true;
  }

  std::optional<TrackMomentumCalculator::Segments>
  TrackMomentumCalculator::getSegTracks_(std::vector<float> const& xxx,
                                         std::vector<float> const& yyy,
                                         std::vector<float> const& zzz,
                                         double const seg_size)
  {
    double stag = 0.0;

    int a1 = xxx.size();
    int a2 = yyy.size();
    int a3 = zzz.size();

    if ((a1 != a2) || (a1 != a3) || (a2 != a3)) {
      cout << " ( Digitize reco tacks ) Error ! " << endl;
      return std::nullopt;
    }

    int const stopper = seg_stop / seg_size;

    std::vector<float> segx, segnx;
    std::vector<float> segy, segny;
    std::vector<float> segz, segnz;
    std::vector<float> segL;

    int ntot = 0;

    n_seg = 0;

    double x0{};
    double y0{};
    double z0{};

    double x00 = xxx.at(0);
    double y00 = yyy.at(0);
    double z00 = zzz.at(0);

    int indC = 0;

    std::vector<float> vx;
    std::vector<float> vy;
    std::vector<float> vz;

    for (int i = 0; i < a1; i++) {
      x0 = xxx.at(i);
      y0 = yyy.at(i);
      z0 = zzz.at(i);

      double const RR0 =
        std::sqrt(cet::sum_of_squares(x00 - x0, y00 - y0, z00 - z0));

      if (RR0 >= stag) {
        segx.push_back(x0);
        segy.push_back(y0);
        segz.push_back(z0);

        segL.push_back(stag);

        x_seg[n_seg] = x0;
        y_seg[n_seg] = y0;
        z_seg[n_seg] = z0;

        n_seg++;

        vx.push_back(x0);
        vy.push_back(y0);
        vz.push_back(z0);

        ntot++;

        indC = i + 1;

        break;
      }
    }

    for (int i = indC; i < a1 - 1; i++) {
      double const x1 = xxx.at(i);
      double const y1 = yyy.at(i);
      double const z1 = zzz.at(i);

      double const dr1 = std::sqrt(cet::sum_of_squares(x1 - x0, y1 - y0, z1 - z0));

      double const x2 = xxx.at(i + 1);
      double const y2 = yyy.at(i + 1);
      double const z2 = zzz.at(i + 1);

      double const dr2 = std::sqrt(cet::sum_of_squares(x2 - x0, y2 - y0, z2 - z0));

      if (dr1 < seg_size) {
        vx.push_back(x1);
        vy.push_back(y1);
        vz.push_back(z1);

        ntot++;
      }

      if (dr1 <= seg_size && dr2 > seg_size) {
        double const dx = x2 - x1;
        double const dy = y2 - y1;
        double const dz = z2 - z1;
        double const dr = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (dr == 0) {
          cout << " ( Zero ) Error ! " << endl;
          return std::nullopt;
        }

        double const beta =
          2.0 * ((x1 - x0) * dx + (y1 - y0) * dy + (z1 - z0) * dz) / (dr * dr);

        double const gamma = (dr1 * dr1 - seg_size * seg_size) / (dr * dr);
        double const delta = beta * beta - 4.0 * gamma;

        if (delta < 0.0) {
          cout << " ( Discriminant ) Error ! " << endl;
          return std::nullopt;
        }

        double const lysi1 = (-beta + std::sqrt(delta)) / 2.0;
        double const t = lysi1;

        double const xp = x1 + t * dx;
        double const yp = y1 + t * dy;
        double const zp = z1 + t * dz;

        segx.push_back(xp);
        segy.push_back(yp);
        segz.push_back(zp);

        segL.push_back(1.0 * n_seg * 1.0 * seg_size + stag);

        x_seg[n_seg] = xp;
        y_seg[n_seg] = yp;
        z_seg[n_seg] = zp;
        n_seg++;

        x0 = xp;
        y0 = yp;
        z0 = zp;

        vx.push_back(x0);
        vy.push_back(y0);
        vz.push_back(z0);

        ntot++;

        auto const na = vx.size();
        double sumx = 0.0;
        double sumy = 0.0;
        double sumz = 0.0;

        for (std::size_t i = 0; i < na; ++i) {
          sumx += vx.at(i);
          sumy += vy.at(i);
          sumz += vz.at(i);
        }

        sumx /= na;
        sumy /= na;
        sumz /= na;

        std::vector<double> mx;
        std::vector<double> my;
        std::vector<double> mz;

        TMatrixDSym m(3);

        for (std::size_t i = 0; i < na; ++i) {
          double const xxw1 = vx.at(i);
          double const yyw1 = vy.at(i);
          double const zzw1 = vz.at(i);

          mx.push_back(xxw1 - sumx);
          my.push_back(yyw1 - sumy);
          mz.push_back(zzw1 - sumz);

          double const xxw0 = mx.at(i);
          double const yyw0 = my.at(i);
          double const zzw0 = mz.at(i);

          m(0, 0) += xxw0 * xxw0 / na;
          m(0, 1) += xxw0 * yyw0 / na;
          m(0, 2) += xxw0 * zzw0 / na;

          m(1, 0) += yyw0 * xxw0 / na;
          m(1, 1) += yyw0 * yyw0 / na;
          m(1, 2) += yyw0 * zzw0 / na;

          m(2, 0) += zzw0 * xxw0 / na;
          m(2, 1) += zzw0 * yyw0 / na;
          m(2, 2) += zzw0 * zzw0 / na;
        }

        TMatrixDSymEigen me(m);

        TVectorD eigenval = me.GetEigenValues();
        TMatrixD eigenvec = me.GetEigenVectors();

        double max1 = -666.0;

        double ind1 = 0;

        for (int i = 0; i < 3; ++i) {
          double const p1 = eigenval(i);

          if (p1 > max1) {
            max1 = p1;
            ind1 = i;
          }
        }

        double ax = eigenvec(0, ind1);
        double ay = eigenvec(1, ind1);
        double az = eigenvec(2, ind1);

        if (n_seg > 1) {
          if (segx.at(n_seg - 1) - segx.at(n_seg - 2) > 0)
            ax = std::abs(ax);
          else
            ax = -1.0 * std::abs(ax);

          if (segy.at(n_seg - 1) - segy.at(n_seg - 2) > 0)
            ay = std::abs(ay);
          else
            ay = -1.0 * std::abs(ay);

          if (segz.at(n_seg - 1) - segz.at(n_seg - 2) > 0)
            az = std::abs(az);
          else
            az = -1.0 * std::abs(az);

          segnx.push_back(ax);
          segny.push_back(ay);
          segnz.push_back(az);
        }
        else
          return std::nullopt;

        ntot = 0;

        vx.clear();
        vy.clear();
        vz.clear();

        vx.push_back(x0);
        vy.push_back(y0);
        vz.push_back(z0);

        ntot++;
      }
      else if (dr1 > seg_size) {
        double const dx = x1 - x0;
        double const dy = y1 - y0;
        double const dz = z1 - z0;
        double const dr = std::sqrt(cet::sum_of_squares(dx, dy, dz));

        if (dr == 0) {
          cout << " ( Zero ) Error ! " << endl;
          return std::nullopt;
        }

        double const t = seg_size / dr;
        double const xp = x0 + t * dx;
        double const yp = y0 + t * dy;
        double const zp = z0 + t * dz;

        segx.push_back(xp);
        segy.push_back(yp);
        segz.push_back(zp);
        segL.push_back(1.0 * n_seg * 1.0 * seg_size + stag);

        x_seg[n_seg] = xp;
        y_seg[n_seg] = yp;
        z_seg[n_seg] = zp;
        n_seg++;

        x0 = xp;
        y0 = yp;
        z0 = zp;

        i = (i - 1);

        vx.push_back(x0);
        vy.push_back(y0);
        vz.push_back(z0);

        ntot++;

        double na = vx.size();

        double sumx = 0.0;
        double sumy = 0.0;
        double sumz = 0.0;

        for (std::size_t i = 0; i < na; ++i) {
          sumx += vx.at(i);
          sumy += vy.at(i);
          sumz += vz.at(i);
        }

        sumx /= na;
        sumy /= na;
        sumz /= na;

        std::vector<double> mx;
        std::vector<double> my;
        std::vector<double> mz;

        TMatrixDSym m(3);

        for (int i = 0; i < na; ++i) {
          double const xxw1 = vx.at(i);
          double const yyw1 = vy.at(i);
          double const zzw1 = vz.at(i);

          mx.push_back(xxw1 - sumx);
          my.push_back(yyw1 - sumy);
          mz.push_back(zzw1 - sumz);

          double const xxw0 = mx.at(i);
          double const yyw0 = my.at(i);
          double const zzw0 = mz.at(i);

          m(0, 0) += xxw0 * xxw0 / na;
          m(0, 1) += xxw0 * yyw0 / na;
          m(0, 2) += xxw0 * zzw0 / na;

          m(1, 0) += yyw0 * xxw0 / na;
          m(1, 1) += yyw0 * yyw0 / na;
          m(1, 2) += yyw0 * zzw0 / na;

          m(2, 0) += zzw0 * xxw0 / na;
          m(2, 1) += zzw0 * yyw0 / na;
          m(2, 2) += zzw0 * zzw0 / na;
        }

        TMatrixDSymEigen me(m);

        TVectorD eigenval = me.GetEigenValues();
        TMatrixD eigenvec = me.GetEigenVectors();

        double max1 = -666.0;
        double ind1 = 0;

        for (int i = 0; i < 3; ++i) {
          double const p1 = eigenval(i);

          if (p1 > max1) {
            max1 = p1;
            ind1 = i;
          }
        }

        double ax = eigenvec(0, ind1);
        double ay = eigenvec(1, ind1);
        double az = eigenvec(2, ind1);

        if (n_seg > 1) {
          if (segx.at(n_seg - 1) - segx.at(n_seg - 2) > 0)
            ax = std::abs(ax);
          else
            ax = -1.0 * std::abs(ax);

          if (segy.at(n_seg - 1) - segy.at(n_seg - 2) > 0)
            ay = std::abs(ay);
          else
            ay = -1.0 * std::abs(ay);

          if (segz.at(n_seg - 1) - segz.at(n_seg - 2) > 0)
            az = std::abs(az);
          else
            az = -1.0 * std::abs(az);

          segnx.push_back(ax);
          segny.push_back(ay);
          segnz.push_back(az);
        }

        else
          return std::nullopt;

        ntot = 0;

        vx.clear();
        vy.clear();
        vz.clear();

        vx.push_back(x0);
        vy.push_back(y0);
        vz.push_back(z0);

        ntot++;
      }

      if (n_seg >= (stopper + 1.0) && seg_stop != -1)
        break;
    }

    delete gr_seg_xyz;
    gr_seg_xyz = new TPolyLine3D{n_seg, z_seg, x_seg, y_seg};
    gr_seg_yz = TGraph{n_seg, z_seg, y_seg};
    gr_seg_xz = TGraph{n_seg, z_seg, x_seg};
    gr_seg_xy = TGraph{n_seg, x_seg, y_seg};

    return std::make_optional<Segments>(Segments{segx, segnx, segy, segny, segz, segnz, segL});
  }

  std::tuple<double, double, double>
  TrackMomentumCalculator::getDeltaThetaRMS_(Segments const& segments,
                                             double const thick) const
  {
    auto const& segnx = segments.nx;
    auto const& segny = segments.ny;
    auto const& segnz = segments.nz;
    auto const& segL = segments.L;

    int const a1 = segnx.size();
    int const a2 = segny.size();
    int const a3 = segnz.size();

    if (a1 != a2 || a1 != a3) {
      cout << " ( Get RMS ) Error ! " << endl;
      return std::make_tuple(0., 0., 0.);
    }

    int const tot = a1 - 1;

    double const thick1 = thick + 0.13;

    std::vector<float> buf0;

    for (int i = 0; i < tot; i++) {
      double const dx = segnx.at(i);
      double const dy = segny.at(i);
      double const dz = segnz.at(i);

      TVector3 const vec_z{dx, dy, dz};
      TVector3 vec_x;
      TVector3 vec_y;

      double const switcher = basex.Dot(vec_z);

      if (switcher <= 0.995) {
        vec_y = vec_z.Cross(basex).Unit();
        vec_x = vec_y.Cross(vec_z);
      }
      else {
        // cout << " It switched ! Isn't this lovely !!! " << endl;
        vec_y = basez.Cross(vec_z).Unit();
        vec_x = vec_y.Cross(vec_z);
      }

      TVector3 const Rx{vec_x.Dot(basex), vec_x.Dot(basey), vec_x.Dot(basez)};
      TVector3 const Ry{vec_y.Dot(basex), vec_y.Dot(basey), vec_y.Dot(basez)};
      TVector3 const Rz{vec_z.Dot(basex), vec_z.Dot(basey), vec_z.Dot(basez)};

      double const refL = segL.at(i);

      for (int j = i; j < tot; j++) {
        double const L1 = segL.at(j);
        double const L2 = segL.at(j + 1);

        double const dz1 = L1 - refL;
        double const dz2 = L2 - refL;

        if (dz1 <= thick1 && dz2 > thick1) {
          double const here_dx = segnx.at(j);
          double const here_dy = segny.at(j);
          double const here_dz = segnz.at(j);

          TVector3 const here_vec{here_dx, here_dy, here_dz};
          TVector3 const rot_here{Rx.Dot(here_vec), Ry.Dot(here_vec), Rz.Dot(here_vec)};

          double const scx = rot_here.X();
          double const scy = rot_here.Y();
          double const scz = rot_here.Z();

          double azy = find_angle(scz, scy);
          azy *= 1.0;

          double const azx = find_angle(scz, scx);

          double const ULim = 10000.0;
          double const LLim = -10000.0;

          if (azx <= ULim && azx >= LLim) {
            buf0.push_back(azx);
          }

          break; // of course !
        }
      }
    }

    int const nmeas = buf0.size();
    double nnn = 0.0;

    double mean = 0.0;
    double rms = 0.0;
    double rmse = 0.0;

    for (int i = 0; i < nmeas; i++) {
      mean += buf0.at(i);
      nnn++;
    }
    mean = mean / nnn;

    for (int i = 0; i < nmeas; i++)
      rms += ((buf0.at(i)) * (buf0.at(i)));

    rms = rms / (nnn);
    rms = std::sqrt(rms);
    rmse = rms / std::sqrt(2.0 * tot);

    double rms1 = rms;

    rms = 0.0;

    double ntot1 = 0.0;
    double const lev1 = 2.50;

    for (int i = 0; i < nmeas; i++) {
      double const amp = buf0.at(i);
      if (amp < (mean + lev1 * rms1) && amp > (mean - lev1 * rms1)) {
        ++ntot1;
        rms += amp * amp;
      }
    }

    rms = rms / (ntot1);
    rms = std::sqrt(rms);
    rmse = rms / std::sqrt(2.0 * ntot1);
    return std::make_tuple(mean, rms, rmse);
  }

  double
  TrackMomentumCalculator::find_angle(double vz, double vy) const
  {
    double thetayz = -999.0;

    if (vz > 0 && vy > 0) {
      double ratio = std::abs(vy / vz);
      thetayz = std::atan(ratio);
    }

    else if (vz < 0 && vy > 0) {
      double ratio = std::abs(vy / vz);
      thetayz = std::atan(ratio);
      thetayz = TMath::Pi() - thetayz;
    }

    else if (vz < 0 && vy < 0) {
      double ratio = std::abs(vy / vz);
      thetayz = std::atan(ratio);
      thetayz = thetayz + TMath::Pi();
    }

    else if (vz > 0 && vy < 0) {
      double ratio = std::abs(vy / vz);
      thetayz = std::atan(ratio);
      thetayz = 2.0 * TMath::Pi() - thetayz;
    }

    else if (vz == 0 && vy > 0) {
      thetayz = TMath::Pi() / 2.0;
    }

    else if (vz == 0 && vy < 0) {
      thetayz = 3.0 * TMath::Pi() / 2.0;
    }

    if (thetayz > TMath::Pi()) {
      thetayz = thetayz - 2.0 * TMath::Pi();
    }

    return 1000.0 * thetayz;
  }

  double
  TrackMomentumCalculator::my_g(double xx, double Q, double s) const
  {
    if (s == 0.) {
      cout << " Error : The code tries to divide by zero ! " << endl;
      return 0.;
    }

    double const arg = (xx - Q) / s;
    double const result =
        -0.5 * std::log(2.0 * TMath::Pi()) - std::log(s) - 0.5 * arg * arg;

    if (std::isnan(result) || std::isinf(result)) {
      cout << " Is nan ! my_g ! " << -std::log(s) << ", " << s << endl;
    }

    return result;
  }

  double
  TrackMomentumCalculator::my_mcs_llhd(std::vector<float> const& dEi,
                                       std::vector<float> const& dEj,
                                       std::vector<float> const& dthij,
                                       std::vector<float> const& ind,
                                       double const x0,
                                       double const x1) const
  {
    double p = x0;
    double theta0x = x1;

    double result = 0.0;
    double nnn1 = dEi.size();

    double red_length = (10.0) / 14.0;
    double addth = 0;

    for (int i = 0; i < nnn1; i++) {
      double Ei = p - dEi.at(i);
      double Ej = p - dEj.at(i);

      if (Ei > 0 && Ej < 0)
        addth = 3.14 * 1000.0;

      Ei = std::abs(Ei);
      Ej = std::abs(Ej);

      double tH0 = (13.6 / std::sqrt(Ei * Ej)) *
                   (1.0 + 0.038 * std::log(red_length)) * std::sqrt(red_length);

      double rms = -1.0;

      if (ind.at(i) == 1) {
        rms = std::sqrt(tH0 * tH0 + pow(theta0x, 2.0));

        double const DT = dthij.at(i) + addth;
        double const prob = my_g(DT, 0.0, rms);

        result = result - 2.0 * prob;
      }
    }

    if (std::isnan(result) || std::isinf(result)) {
      std::cout << " Is nan ! my_mcs_llhd ( 1 ) ! " << std::endl;
    }
    return result;
  }

} // namespace track
