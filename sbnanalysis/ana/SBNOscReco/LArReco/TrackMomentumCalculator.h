/// \file  TrackMomentumCalculator.h
//  \author  sowjanyag@phys.ksu.edu

#ifndef TrackMomentumCalculator_H
#define TrackMomentumCalculator_H

#include "canvas/Persistency/Common/Ptr.h"
#include "larcorealg/CoreUtils/quiet_Math_Functor.h"
#include "lardataobj/RecoBase/Track.h"

#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TPolyLine3D.h"
#include "TSpline.h"
#include "TVector3.h"

#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserParameters.h"

#include <cmath>
#include <iostream>
#include <optional>
#include <vector>
#include <tuple>

namespace trkf {

  class TrackMomentumCalculator {
  public:
    TrackMomentumCalculator(double minLength = 100.0,
                            double maxLength = 1350.0);

    double GetTrackMomentum(double trkrange, int pdg) const;
    double GetMomentumMultiScatterChi2(art::Ptr<recob::Track> const& trk);
    double GetMomentumMultiScatterLLHD(art::Ptr<recob::Track> const& trk);
    double GetMuMultiScatterLLHD3(art::Ptr<recob::Track> const& trk, bool dir);
    TVector3 GetMultiScatterStartingPoint(art::Ptr<recob::Track> const& trk);

  private:
    bool plotRecoTracks_(std::vector<float> const& xxx,
                         std::vector<float> const& yyy,
                         std::vector<float> const& zzz);

    struct Segments {
      std::vector<float> x, nx;
      std::vector<float> y, ny;
      std::vector<float> z, nz;
      std::vector<float> L;
    };

    std::optional<Segments> getSegTracks_(std::vector<float> const& xxx,
                                          std::vector<float> const& yyy,
                                          std::vector<float> const& zzz,
                                          double seg_size);

    std::tuple<double, double, double> getDeltaThetaRMS_(Segments const& segments,
                                                         double thick) const;

    int getDeltaThetaij_(std::vector<float>& ei,
                         std::vector<float>& ej,
                         std::vector<float>& th,
                         std::vector<float>& ind,
                         Segments const& segments,
                         double thick) const;

    double my_g(double xx, double Q, double s) const;

    double my_mcs_llhd(std::vector<float> const& dEi,
                       std::vector<float> const& dEj,
                       std::vector<float> const& dthij,
                       std::vector<float> const& ind,
                       double x0, double x1) const;

    float seg_stop{-1.};
    int n_seg{};

    float x_seg[50000];
    float y_seg[50000];
    float z_seg[50000];

    double find_angle(double vz, double vy) const;

    float steps_size{10.};
    int n_steps{6};
    std::vector<float> steps;

    double minLength;
    double maxLength;

    // The following are objects that are created but not drawn or
    // saved.  This class should consider accepting a "debug"
    // parameter where if it is specified, then the graphs will be
    // created; otherwise, their creation is unnecessary and impedes
    // efficiency.
    //
    // N.B. TPolyLine3D objects are owned by ROOT, and we thus refer
    // to them by pointer.  It is important that 'delete' is not
    // called on the TPolyLine3D pointers during destruction of a
    // TrackMomentumCalculator object.
    TPolyLine3D* gr_reco_xyz{nullptr};
    TGraph gr_reco_xy{};
    TGraph gr_reco_yz{};
    TGraph gr_reco_xz{};

    TPolyLine3D* gr_seg_xyz{nullptr};
    TGraph gr_seg_xy{};
    TGraph gr_seg_yz{};
    TGraph gr_seg_xz{};

  };

} // namespace trkf

#endif // TrackMomentumCalculator_H
