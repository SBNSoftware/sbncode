#include "TrajectoryMCSFitter.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

using namespace std;
using namespace trkf::sbn;
using namespace recob::tracking;

recob::MCSFitResult TrajectoryMCSFitter::fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst) const {
  //
  // Break the trajectory in segments of length approximately equal to segLen_
  //
  vector<size_t> breakpoints;
  vector<float> segradlengths;
  vector<float> cumseglens;
  vector<bool> breakpointsgood;
  breakTrajInSegments(traj, breakpoints, segradlengths, cumseglens, breakpointsgood);
  //
  // Fit segment directions, and get 3D angles between them
  //
  if (segradlengths.size()<2) return recob::MCSFitResult();
  vector<float> dtheta;
  Vector_t pcdir0;
  Vector_t pcdir1;
  for (unsigned int p = 0; p<segradlengths.size(); p++) {
    linearRegression(traj, breakpoints[p], breakpoints[p+1], pcdir1);
    if (p>0) {
      if (segradlengths[p]<-100. || segradlengths[p-1]<-100.) {
        dtheta.push_back(-999.);
      } 
      if (!breakpointsgood[p] || !breakpointsgood[p-1]) {
        dtheta.push_back(-999.);
      }
      else { 
        const double cosval = pcdir0.X()*pcdir1.X()+pcdir0.Y()*pcdir1.Y()+pcdir0.Z()*pcdir1.Z();
        //assert(std::abs(cosval)<=1);
        //units are mrad
        double dt = 1000.*acos(cosval);//should we try to use expansion for small angles?
        dtheta.push_back(dt);
      }
    }
    pcdir0 = pcdir1;
  }
  //
  // Perform likelihood scan in forward and backward directions
  //
  vector<float> cumLenFwd;
  vector<float> cumLenBwd;
  for (unsigned int i = 0; i<cumseglens.size()-2; i++) {
    cumLenFwd.push_back(cumseglens[i]);
    cumLenBwd.push_back(cumseglens.back()-cumseglens[i+2]);
  }
  const ScanResult fwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenFwd, true,  momDepConst, pid);
  const ScanResult bwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenBwd, false, momDepConst, pid);
  // std::cout << "fwdResult.p=" << fwdResult.p << " fwdResult.pUnc=" << fwdResult.pUnc << " fwdResult.logL=" << fwdResult.logL << std::endl;
  return recob::MCSFitResult(pid,
                            fwdResult.p,fwdResult.pUnc,fwdResult.logL,
                            bwdResult.p,bwdResult.pUnc,bwdResult.logL,
                            segradlengths,dtheta);
}

void TrajectoryMCSFitter::breakTrajInSegments(const recob::TrackTrajectory& traj, vector<size_t>& breakpoints, vector<float>& segradlengths, vector<float>& cumseglens, vector<bool>& breakpointsgood) const {
  // set fiducial volume
  std::vector<geo::BoxBoundedGeo> fiducialVolumes;
  fiducialVolumes = setFiducialVolumes();
  // set excluded volumes
  std::vector<geo::BoxBoundedGeo> excludeVolumes;
  excludeVolumes = setExcludeVolumes();

  const double trajlen = traj.Length();
  const int nseg = std::max(minNSegs_,int(trajlen/segLen_));
  const double thisSegLen = trajlen/double(nseg);
  // std::cout << "track with length=" << trajlen << " broken in nseg=" << nseg << " of length=" << thisSegLen << " where segLen_=" << segLen_ << std::endl;
  //
  constexpr double lar_radl_inv = 1./14.0;
  cumseglens.push_back(0.);//first segment has zero cumulative length from previous segments
  double thislen = 0.;
  auto nextValid=traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  auto pos0 = traj.LocationAtPoint(nextValid);
  bool pos0good = isInVolume(fiducialVolumes, pos0) && !isInVolume(excludeVolumes, pos0);
  breakpointsgood.push_back(pos0good);
  nextValid = traj.NextValidPoint(nextValid+1);
  int npoints = 0;
  while (nextValid!=recob::TrackTrajectory::InvalidIndex) {
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += ( (pos1-pos0).R() );
    pos0=pos1;
    npoints++;
    if (thislen>=thisSegLen) {
      bool pos1good = isInVolume(fiducialVolumes, pos1) && !isInVolume(excludeVolumes, pos1);
      breakpoints.push_back(nextValid);
      breakpointsgood.push_back(pos1good);
      if (npoints>=minHitsPerSegment_) segradlengths.push_back(thislen*lar_radl_inv);
      else segradlengths.push_back(-999.);
      cumseglens.push_back(cumseglens.back()+thislen);
      thislen = 0.;
      npoints = 0;
    }
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //then add last segment
  if (thislen>0.) {
    auto endpointpos = traj.LocationAtPoint(nextValid);
    bool endpointposgood = isInVolume(fiducialVolumes, endpointpos) && !isInVolume(excludeVolumes, endpointpos);
    breakpoints.push_back(traj.LastValidPoint()+1);
    breakpointsgood.push_back(endpointposgood);
    segradlengths.push_back(thislen*lar_radl_inv);
    cumseglens.push_back(cumseglens.back()+thislen);
  }
  return;
}

const TrajectoryMCSFitter::ScanResult TrajectoryMCSFitter::doLikelihoodScan(std::vector<float>& dtheta, std::vector<float>& seg_nradlengths, std::vector<float>& cumLen, bool fwdFit, bool momDepConst, int pid) const {
  int    best_idx  = -1;
  double best_logL = std::numeric_limits<double>::max();
  double best_p    = -1.0;
  std::vector<float> vlogL;
  for (double p_test = pMin_; p_test <= pMax_; p_test+=pStep_) {
    double logL = mcsLikelihood(p_test, angResol_, dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst, pid);
    if (logL < best_logL) {
      best_p    = p_test;
      best_logL = logL;
      best_idx  = vlogL.size();
    }
    vlogL.push_back(logL);
  }
  //
  //uncertainty from left side scan
  double lunc = -1.0;
  if (best_idx>0) {
    for (int j=best_idx-1;j>=0;j--) {
      double dLL = vlogL[j]-vlogL[best_idx];
      if ( dLL<0.5 ) {
        lunc = (best_idx-j)*pStep_;
      } else break;
    }
  }
  //uncertainty from right side scan
  double runc = -1.0;
  if (best_idx<int(vlogL.size()-1)) {  
    for (unsigned int j=best_idx+1;j<vlogL.size();j++) {
      double dLL = vlogL[j]-vlogL[best_idx];
      if ( dLL<0.5 ) {
        runc = (j-best_idx)*pStep_;
      } else break;
    }
  }
  return ScanResult(best_p, std::max(lunc,runc), best_logL);
}

void TrajectoryMCSFitter::linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
  //
  int npoints = 0;
  geo::vect::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  while (nextValid<lastPoint) {
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    nextValid = traj.NextValidPoint(nextValid+1);
    npoints++;
  }
  const auto avgpos = middlePointCalc.middlePoint();
  const double norm = 1./double(npoints);
  //
  //assert(npoints>0);
  //
  TMatrixDSym m(3);
  nextValid = firstPoint;
  while (nextValid<lastPoint) {
    const auto p = traj.LocationAtPoint(nextValid);
    const double xxw0 = p.X()-avgpos.X();
    const double yyw0 = p.Y()-avgpos.Y();
    const double zzw0 = p.Z()-avgpos.Z();
    m(0, 0) += xxw0*xxw0*norm;
    m(0, 1) += xxw0*yyw0*norm;
    m(0, 2) += xxw0*zzw0*norm;
    m(1, 0) += yyw0*xxw0*norm;
    m(1, 1) += yyw0*yyw0*norm;
    m(1, 2) += yyw0*zzw0*norm;
    m(2, 0) += zzw0*xxw0*norm;
    m(2, 1) += zzw0*yyw0*norm;
    m(2, 2) += zzw0*zzw0*norm;
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i=1; i<3; ++i) {
    if (eigenval(i)>maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  } 
  //
  pcdir = Vector_t(eigenvec(0, maxevalidx),eigenvec(1, maxevalidx),eigenvec(2, maxevalidx));
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir)<0.) pcdir*=-1.;
  //
}

double TrajectoryMCSFitter::mcsLikelihood(double p, double theta0x, std::vector<float>& dthetaij, std::vector<float>& seg_nradl, std::vector<float>& cumLen, bool fwd, bool momDepConst, int pid) const {
  //
  const int beg  = (fwd ? 0 : (dthetaij.size()-1));
  const int end  = (fwd ? dthetaij.size() : -1);
  const int incr = (fwd ? +1 : -1);
  //
  // bool print = false;//(p>1.999 && p<2.001);
  //
  const double m = mass(pid);
  const double m2 = m*m;
  const double Etot = sqrt(p*p + m2);//Initial energy
  double Eij2 = 0.;
  //
  double const fixedterm = 0.5 * std::log( 2.0 * M_PI );
  double result = 0;
  for (int i = beg; i != end; i+=incr ) {
    if (dthetaij[i]<0) {
      continue;
    }
    //
    if (eLossMode_==1) {
      // ELoss mode: MIP (constant)
      constexpr double kcal = 0.002105;
      const double Eij = Etot - kcal*cumLen[i];//energy at this segment
      Eij2 = Eij*Eij;
    } else {
      // Non constant energy loss distribution
      const double Eij = GetE(Etot,cumLen[i],m);
      Eij2 = Eij*Eij;
    }
    //
    if ( Eij2 <= m2 ) {
      result = std::numeric_limits<double>::max();
      break;
    }
    const double pij = sqrt(Eij2 - m2);//momentum at this segment
    const double beta = sqrt( 1. - ((m2)/(pij*pij + m2)) );
    constexpr double tuned_HL_term1 = 11.0038; // https://arxiv.org/abs/1703.06187
    constexpr double HL_term2 = 0.038;
    const double tH0 = ( (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) / (pij*beta) ) * ( 1.0 + HL_term2 * std::log( seg_nradl[i] ) ) * sqrt( seg_nradl[i] );
    const double rms = sqrt( 2.0*( tH0 * tH0 + theta0x * theta0x ) );
    if (rms==0.0) {
      std::cout << " Error : RMS cannot be zero ! " << std::endl;
      return std::numeric_limits<double>::max();
    } 
    const double arg = dthetaij[i]/rms;
    result += ( std::log( rms ) + 0.5 * arg * arg + fixedterm);
    // if (print && fwd==true) cout << "TrajectoryMCSFitter pij=" << pij << " dthetaij[i]=" << dthetaij[i] << " tH0=" << tH0 << " rms=" << rms << " prob=" << ( std::log( rms ) + 0.5 * arg * arg + fixedterm) << " const=" << (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) << " beta=" << beta << " red_length=" << seg_nradl[i] << endl;
  }
  return result;
}

double TrajectoryMCSFitter::energyLossLandau(const double mass2,const double e2, const double x) const {
  //
  // eq. (33.11) in http://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (except density correction is ignored)
  //
  if (x<=0.) return 0.;
  constexpr double Iinv2 = 1./(188.E-6*188.E-6);
  constexpr double matConst = 1.4*18./40.;//density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  constexpr double j = 0.200;
  //
  const double beta2 = (e2-mass2)/e2;
  const double gamma2 = 1./(1.0 - beta2);
  const double epsilon = 0.5*kappa*x*matConst/beta2;
  //
  return 0.001*epsilon*( log(2.*me*beta2*gamma2*epsilon*Iinv2) + j - beta2 );
}
//
double TrajectoryMCSFitter::energyLossBetheBloch(const double mass,const double e2) const {
  // stolen, mostly, from GFMaterialEffects.
  constexpr double Iinv = 1./188.E-6;
  constexpr double matConst = 1.4*18./40.;//density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  //
  const double beta2 = (e2-mass*mass)/e2;
  const double gamma2 = 1./(1.0 - beta2);
  const double massRatio = me/mass;
  const double argument = (2.*me*gamma2*beta2*Iinv) * std::sqrt(1+2*std::sqrt(gamma2)*massRatio + massRatio*massRatio);
  //
  double dedx = kappa*matConst/beta2;
  //
  if (mass==0.0) return(0.0);
  if (argument <= exp(beta2)) {
      dedx = 0.;
  } else{
    dedx *= (log(argument)-beta2)*1.E-3; // Bethe-Bloch, converted to GeV/cm
    if (dedx<0.) dedx = 0.;
  }
  return dedx;
}
//
double TrajectoryMCSFitter::GetE(const double initial_E, const double length_travelled, const double m) const {
  //
  const double step_size = length_travelled / nElossSteps_;
  //
  double current_E = initial_E;
  const double m2 = m*m;
  //
  for (auto i = 0; i < nElossSteps_; ++i) {
    if (eLossMode_==2) {
      double dedx = energyLossBetheBloch(m,current_E);
      current_E -= (dedx * step_size);
    } else {
      // MPV of Landau energy loss distribution
      current_E -= energyLossLandau(m2,current_E*current_E,step_size);
    }
    if ( current_E <= m ) {
      // std::cout<<"WARNING: current_E less than mu mass. it is "<<current_E<<std::endl;
      return 0.;
    }
  }
  return current_E;
}

std::vector<geo::BoxBoundedGeo> TrajectoryMCSFitter::setFiducialVolumes() const {
  std::vector<geo::BoxBoundedGeo> fiducialVolumes;
  std::vector<std::vector<geo::BoxBoundedGeo>> TPCVolumes;
  std::vector<geo::BoxBoundedGeo> ActiveVolumes;
  
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  for (auto const &cryoID: geometry->Iterate<geo::CryostatID>()) {
    std::vector<geo::BoxBoundedGeo> thisTPCVolumes;
    for (auto const& tpc: geometry->Iterate<geo::TPCGeo>(cryoID)) {
      thisTPCVolumes.push_back(tpc.ActiveBoundingBox());
    }
    TPCVolumes.push_back(std::move(thisTPCVolumes));
  }
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: TPCVolumes) {
    double xMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double xMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double yMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double yMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double zMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();
    double zMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();
    ActiveVolumes.emplace_back(xMin, xMax, yMin, yMax, zMin, zMax);
  }
  double fidInsetMinX = fiducialVolumeInsets_[0];
  double fidInsetMaxX = fiducialVolumeInsets_[1];
  double fidInsetMinY = fiducialVolumeInsets_[2];
  double fidInsetMaxY = fiducialVolumeInsets_[3];
  double fidInsetMinZ = fiducialVolumeInsets_[4];
  double fidInsetMaxZ = fiducialVolumeInsets_[5];

  if (fiducialVolumeInsets_.size() != 6) {
    std::cout << "Error: fiducialVolumeInsets vector must have length of 6, not fiducializing" << std::endl;
    fidInsetMinX = 0.0;
    fidInsetMaxX = 0.0;
    fidInsetMinY = 0.0;
    fidInsetMaxY = 0.0;
    fidInsetMinZ = 0.0;
    fidInsetMaxZ = 0.0;
  }
  for (const geo::BoxBoundedGeo &AV: ActiveVolumes) {
    fiducialVolumes.emplace_back(AV.MinX() + fidInsetMinX, AV.MaxX() - fidInsetMaxX,
                                 AV.MinY() + fidInsetMinY, AV.MaxY() - fidInsetMaxY,
                                 AV.MinZ() + fidInsetMinZ, AV.MaxZ() - fidInsetMaxZ);
  }
  return fiducialVolumes;
}

std::vector<geo::BoxBoundedGeo> TrajectoryMCSFitter::setExcludeVolumes() const {
  std::vector<geo::BoxBoundedGeo> excludeVolumes;
  if (excludeVolumes_.size()%6 != 0) {
    std::cout << "Error: excludeVolumes vector must have length multiple of 6, not excluding any regions" << std::endl;
    excludeVolumes.emplace_back(-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0);
    return excludeVolumes;
  }
  for (unsigned int i=0; i<excludeVolumes_.size()/6; i++) {
    excludeVolumes.emplace_back(excludeVolumes_[6*i+0], excludeVolumes_[6*i+1], excludeVolumes_[6*i+2], excludeVolumes_[6*i+3], excludeVolumes_[6*i+4], excludeVolumes_[6*i+5]);
  }
  return excludeVolumes;
}

bool TrajectoryMCSFitter::isInVolume(const std::vector<geo::BoxBoundedGeo> &volumes, const geo::Point_t &point) const {
  for (const geo::BoxBoundedGeo &volume: volumes) {
    if (point.X()>=volume.MinX() && point.X()<=volume.MaxX() &&
        point.Y()>=volume.MinY() && point.Y()<=volume.MaxY() &&
        point.Z()>=volume.MinZ() && point.Z()<=volume.MaxZ()) {
      return true;
    }
  }
  return false;
}