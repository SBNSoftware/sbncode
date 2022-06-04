#ifndef TRAJECTORYMCSFITTER_H
#define TRAJECTORYMCSFITTER_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/RecoObjects/TrackState.h"

namespace trkf::sbn {
  /**
   * @file  larreco/RecoAlg/TrajectoryMCSFitter.h
   * @class trkf::TrajectoryMCSFitter
   *
   * @brief Class for Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory
   *
   * Class for Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory.
   *
   * Inputs are: a Track or Trajectory, and various fit parameters (pIdHypothesis, minNumSegments, segmentLength, pMin, pMax, pStep, angResol)
   *
   * Outputs are: a recob::MCSFitResult, containing:
   *   resulting momentum, momentum uncertainty, and best likelihood value (both for fwd and bwd fit);
   *   vector of segment (radiation) lengths, vector of scattering angles, and PID hypothesis used in the fit.
   *
   * For configuration options see TrajectoryMCSFitter#Config
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */
  class TrajectoryMCSFitter {
    // 
  public:
    //
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> pIdHypothesis {
        Name("pIdHypothesis"),
	Comment("Default particle Id Hypothesis to be used in the fit when not specified."),
	13
      };
      fhicl::Atom<int> minNumSegments {
        Name("minNumSegments"),
	Comment("Minimum number of segments the track is split into."),
	3
      };
      fhicl::Atom<double> segmentLength {
        Name("segmentLength"),
	Comment("Nominal length of track segments used in the fit."),
	14.
      };
      fhicl::Atom<int> minHitsPerSegment {
        Name("minHitsPerSegment"),
	Comment("Exclude segments with less hits than this value."),
	2
      };
      fhicl::Atom<int> nElossSteps {
        Name("nElossSteps"),
	Comment("Number of steps for computing energy loss uptream to current segment."),
	10
      };
      fhicl::Atom<int> eLossMode {
        Name("eLossMode"),
	Comment("Default is MPV Landau. Choose 1 for MIP (constant); 2 for Bethe-Bloch."),
	0
      };
      fhicl::Atom<double> pMin {
        Name("pMin"),
	Comment("Minimum momentum value in likelihood scan."),
	0.01
      };
      fhicl::Atom<double> pMax {
        Name("pMax"),
	Comment("Maximum momentum value in likelihood scan."),
	7.50  
      };
      fhicl::Atom<double> pStep {
        Name("pStep"),
	Comment("Step in momentum value in likelihood scan."),
	0.01
      };
      fhicl::Atom<double> angResol {
        Name("angResol"),
	Comment("Angular resolution parameter used in modified Highland formula. Unit is mrad."),
	3.0
      };
    };
    using Parameters = fhicl::Table<Config>;
    //
    TrajectoryMCSFitter(int pIdHyp, int minNSegs, double segLen, int minHitsPerSegment, int nElossSteps, int eLossMode, double pMin, double pMax, double pStep, double angResol){
      pIdHyp_ = pIdHyp;
      minNSegs_ = minNSegs;
      segLen_ = segLen;
      minHitsPerSegment_ = minHitsPerSegment;
      nElossSteps_ = nElossSteps;
      eLossMode_ = eLossMode;
      pMin_ = pMin;
      pMax_ = pMax;
      pStep_ = pStep;
      angResol_ = angResol;
    }
    explicit TrajectoryMCSFitter(const Parameters & p)
      : TrajectoryMCSFitter(p().pIdHypothesis(),p().minNumSegments(),p().segmentLength(),p().minHitsPerSegment(),p().nElossSteps(),p().eLossMode(),p().pMin(),p().pMax(),p().pStep(),p().angResol()) {}
    //
    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj, bool momDepConst = true) const { return fitMcs(traj,pIdHyp_,momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Track& track,          bool momDepConst = true) const { return fitMcs(track,pIdHyp_,momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Trajectory& traj,      bool momDepConst = true) const { return fitMcs(traj,pIdHyp_,momDepConst); }
    //
    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst = true) const;
    recob::MCSFitResult fitMcs(const recob::Track& track,          int pid, bool momDepConst = true) const { return fitMcs(track.Trajectory(),pid,momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Trajectory& traj,      int pid, bool momDepConst = true) const {
      recob::TrackTrajectory::Flags_t flags(traj.NPoints());
      const recob::TrackTrajectory tt(traj,std::move(flags));
      return fitMcs(tt,pid,momDepConst);
    }
    //
    void breakTrajInSegments(const recob::TrackTrajectory& traj, std::vector<size_t>& breakpoints, std::vector<float>& segradlengths, std::vector<float>& cumseglens) const;
    void linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, recob::tracking::Vector_t& pcdir) const;
    double mcsLikelihood(double p, double theta0x, std::vector<float>& dthetaij, std::vector<float>& seg_nradl, std::vector<float>& cumLen, bool fwd, bool momDepConst, int pid) const;
    //
    struct ScanResult {
      public:
        ScanResult(double ap, double apUnc, double alogL) : p(ap), pUnc(apUnc), logL(alogL) {}
        double p, pUnc, logL;
    };
    //
    const ScanResult doLikelihoodScan(std::vector<float>& dtheta, std::vector<float>& seg_nradlengths, std::vector<float>& cumLen, bool fwdFit, bool momDepConst, int pid) const;
    //
    inline double MomentumDependentConstant(const double p) const {
      //these are from https://arxiv.org/abs/1703.06187
      constexpr double a = 0.1049;
      constexpr double c = 11.0038;
      return (a/(p*p)) + c;
    }
    double mass(int pid) const {
      if (abs(pid)==13)   { return mumass; }
      if (abs(pid)==211)  { return pimass; }
      if (abs(pid)==321)  { return kmass;  }
      if (abs(pid)==2212) { return pmass;  }
      return util::kBogusD;
    }
    double energyLossBetheBloch(const double mass,const double e2) const;
    double energyLossLandau(const double mass2,const double E2, const double x) const;
    //
    double GetE(const double initial_E, const double length_travelled, const double mass) const;
    //
  private:
    int    pIdHyp_;
    int    minNSegs_;
    double segLen_;
    int    minHitsPerSegment_;
    int    nElossSteps_;
    int    eLossMode_;
    double pMin_;
    double pMax_;
    double pStep_;
    double angResol_;
  };
}

#endif
