
#include <TLorentzVector.h>

// Mass values from GENIE v3.0.6
constexpr double NEUTRON_MASS = 0.93956541; // GeV
constexpr double PROTON_MASS = 0.93827208; // GeV
constexpr double MUON_MASS = 0.10565837; // GeV

constexpr double EB_N = 0.02478; // GeV/nucleon in 40Ar
constexpr double EB_P = 0.02122; // GeV/nucleon in 40Ar


// compute and return the four-vector for the final state hadronic system given a list of 3-momenta 
// assuming that they are all for protons
TLorentzVector GetHadVec(const std::vector<TVector3>& momenta) {

  TLorentzVector outvec(0,0,0,0); // initialized to all 0's
  if(momenta.empty()) { // this could happen if the event is not signal, e.g. no visible protons
    return outvec;
  }

  // loop over proton 3-momenta
  for(TVector3 const& vec : momenta) {
    if(isinf(vec.Mag())) {
      std::cout << "mom vec with inf mag!"
                << "\n\t(px,py,pz) = (" << vec.X()
                << ", " << vec.Y() << ", "
                << vec.Z() << ")" << std::endl;
      continue;
    }
    // mass from knock-out protons already existed in initial state so doesn't contribute to E_h
    double binde = 0.008595261375; //binding energy per nucleon for 40-Ar in GeV
    if(outvec.M() == 0){
      outvec += TLorentzVector(vec,sqrt(vec.Mag2()+pow(PROTON_MASS,2))+binde); 
      if(outvec.M() == 0)std::cout << "ERROR: outvec is still empty after initial proton" <<std::endl;
    }
    else{
      outvec += TLorentzVector(vec,sqrt(vec.Mag2()+pow(PROTON_MASS,2))-PROTON_MASS+binde);
      if(outvec.M() == 0)std::cout << "ERROR: outvec has become empty adding in knockout protons!" <<std::endl;
    }
  }

  // some sanity checks
  if(outvec.M()==0) {
   std::cout << "WARNING in KinematicReco::GetHadVec W is empty!" <<std::endl;
  }  

  if(isnan(outvec.M()) || isinf(outvec.M())) {
    std::cout << "ERROR in GetHadVec: W is NaN or is INF!" <<std::endl;
    for(TVector3 const& vec : momenta)          
       std::cout << "\t vec.Mag() = " << vec.Mag() << std::endl;
  }

  /*if(outvec.M()<0){
    std::cout << "WARNING in GetHadVec: negative invariant mass W"
              << "\n\tp_h = (" << outvec.X() << "," << outvec.Y()
              << "," << outvec.Z() << "," << outvec.T() << ")"
              << " -> W = " << outvec.M() << " GeV"
              << "\n\t with " << momenta.size() << " protons"
              << std::endl;
  }*/

  return outvec;

}

// compute and return the four-vector for the final state hadronic system given a list of 3-momenta 
// assuming that they are all for protons but doing a more careful correction for proton mass than above
// the option wfix is whether to start by fixing W to be Mp (wfix=true) or obtain W from
// phad + dispersion rel. , correct it for (N-1)Mp, then back-correct Ehad (default)
TLorentzVector GetHadVecW(const std::vector<TVector3>& momenta, bool wfix=false) {

  TLorentzVector outvec(0,0,0,0); // initialized to all 0's
  if(momenta.empty()) { // this could happen if the event is not signal, e.g. no visible protons
    return outvec;
  }

  // loop over proton 3-momenta
  for(TVector3 const& vec : momenta) {
    if(isinf(vec.Mag())) {
      std::cout << "mom vec with inf mag!"
                << "\n\t(px,py,pz) = (" << vec.X()
                << ", " << vec.Y() << ", "
                << vec.Z() << ")" << std::endl;
      continue;
    }
    outvec += TLorentzVector(vec,sqrt(vec.Mag2()+pow(PROTON_MASS,2)));
  }

  //spectator proton mass correction
  double w = 0.;
  if(wfix) 
    w = PROTON_MASS;
  else 
    w = outvec.M()-(momenta.size()-1)*PROTON_MASS; 

  outvec.SetE(sqrt(pow(outvec.P(),2)+w*w));

  // some sanity checks
  if(outvec.M()==0) {
   std::cout << "WARNING in KinematicReco::GetHadVec W is empty!" <<std::endl;
  }

  if(isnan(outvec.M()) || isinf(outvec.M())) {
    std::cout << "ERROR in GetHadVec: W is NaN or is INF!" <<std::endl;
    for(TVector3 const& vec : momenta)
       std::cout << "\t vec.Mag() = " << vec.Mag() << std::endl;
  }

  if(outvec.M()<0){
    std::cout << "WARNING in GetHadVec: negative invariant mass W"
              << "\n\tp_h = (" << outvec.X() << "," << outvec.Y()
              << "," << outvec.Z() << "," << outvec.T() << ")"
              << " -> W = " << outvec.M() << " GeV"
              << "\n\t with " << momenta.size() << " protons"
              << std::endl;
  }

  return outvec;

}

// compute and return the neutrino energy [GeV] with the kinematic method
// taking the final state hadronic system 4-vector as input.
// N.B. this assumes nu_mu CC interactions (not anti-nu_mu) where the 
// neutron mass would need to be swapped out for the proton mass.
double KinRecoEnu(const TLorentzVector& hadvec, const bool cc=false, const bool rhc=false) {

  double enu = pow(hadvec.P(),2) - pow(hadvec.E()-NEUTRON_MASS,2) + pow(MUON_MASS,2);
  enu *= 0.5/(hadvec.Pz() + NEUTRON_MASS - hadvec.E());

  if(isnan(enu)) std::cout << "KinRecoEnu gives NaN!" << std::endl;
  if(isinf(enu)) std::cout << "KinRecoEnu gives inf!" << std::endl;
  if(isnan(enu) || isinf(enu)) {
    std::cout << "\tW = " << hadvec.M()
              << "\n\tE_h = " << hadvec.E()
              << "\n\tP_h,z = " << hadvec.Pz() << std::endl;
  }

  return enu;
}

double KinRecoEnu(const double& el, const double& eh) {

  return el + eh - NEUTRON_MASS;
}

// brief: return the outgoing lepton total energy [GeV] from the kinematic method
// @param hadvec : FS hadronic 4-momentum
// @param cc : whether the interaction is CC (not NC)
// @param rhc : whether the beam is RHC and anti-neutrino dominated (not FHC/neutrino dominated)
double KinRecoEl(const TLorentzVector& hadvec, const bool cc=false, const bool rhc=false) {

  double mnuc = NEUTRON_MASS, ml =0.;
  if(cc) ml = MUON_MASS;
  if(rhc) mnuc = PROTON_MASS;

  double el = hadvec.M2() - mnuc*mnuc - ml*ml
              + 2*(mnuc-hadvec.E()) * (hadvec.E()-hadvec.Pz());
  el *= 0.5/(hadvec.E() - hadvec.Pz() - mnuc);

  if(isnan(el)) std::cout << "KinRecoE gives NaN!" << std::endl;
  if(isinf(el)) std::cout << "KinRecoE gives inf!" << std::endl;
  if(isnan(el) || isinf(emu)) {
    std::cout << "\tW = " << hadvec.M()
              << "\n\tE_h = " << hadvec.E()
              << "\n\tP_h,z = " << hadvec.Pz() << std::endl;
  }

  return el;
}

// compute and return the muon momentum [GeV/c] with the kinematic method
// taking the final state hadronic system 4-vector as input.
// N.B. this assumes nu_mu CC interactions (not anti-nu_mu) where the 
// neutron mass would need to be swapped out for the proton mass.
double KinRecoPl(const TLorentzVector& hadvec, const bool cc=false, const bool rhc=false) {

  double pl = 0., ml=0.;
  if(cc) ml = MUON_MASS;

  const double el = KinRecoEl(hadvec,cc,rhc);
  if(abs(el) < ml) 
    pl = -1. * sqrt(ml*ml - el*el); //pure imaginary -> negative

  else if(el<0)
    pl = -1*sqrt(el*el - ml*ml);

  else
    pl = sqrt(el*el - ml*ml);

  if(isnan(pl)){
    std::cout << "ERROR IN KinRecoPl: Pl is NaN! El = " << el << std::endl;
  }
  return pl;
}

// compute and return the muon polar angle [rad] with the kinematic method
// taking the final state hadronic system 4-vector as input
// N.B. this assumes nu_mu CC interactions (not anti-nu_mu) where the 
// neutron mass would need to be swapped out for the proton mass
double KinRecoCosTheta(const TLorentzVector& hadvec) {

  double elep = KinRecoE(hadvec);
  double den = elep*elep-pow(MUON_MASS,2);
  if(den<0)
      den = -1.*sqrt(-1.*den); //pure imaginary -> negative
  else
      den = sqrt(den);
  return (elep+hadvec.E()-NEUTRON_MASS-hadvec.Pz())/den;
}


// compute and return the muon polar angle [rad] with the kinematic method
// taking the final state hadronic system 4-vector as input
// N.B. this assumes nu_mu CC interactions (not anti-nu_mu) where the 
// neutron mass would need to be swapped out for the proton mass
double KinRecoTheta(const TLorentzVector& hadvec) {

  double costheta = KinRecoCosTheta(hadvec);
  if(costheta>1.000000000001 || costheta<-1.000000000001) {
      return -DBL_MAX;
  }

  return acos(costheta);
}
