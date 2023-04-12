#include "AnisotropicThreeBodyDecay.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

#include "TLorentzVector.h"
#include "TRandom3.h"

/*Implementation of HNL Three Body Decays Anisotropies
  Valid as long as the HNL is decaying into a neutrino and identical final-state charged leptons
  See arXiv:2104.05719 for more details @LuisPelegrina */


/* Constants used for the calculations of C_i Coefficients given the electro-weak mixing angle
   See arXiv:2104.05719 (Table 3) for more details */

double evgen::ldm::AnThreeBD::GL(bool CC)
{
  double const sw2 = Constants::Instance().sin2thetaW;
  if(CC)  return 0.5 * (1 + 2 * sw2); //When the decay can happen by CC gL changes
  else  return -0.5 * (1 - 2 * sw2);
}

double evgen::ldm::AnThreeBD::GR()
{
  double const sw2 = Constants::Instance().sin2thetaW;
  return sw2;
}


/* Coefficients of Lorentz-Invariant Objects given couplings gL and gR (Majorana HNL)
   See arXiv:2104.05719 (Table 3) for more details
   Valid as long as the HNL is decaying into a neutrino and identical final-state charged leptons */
void evgen::ldm::AnThreeBD::CMgLgR(double (&CM)[6], double Ue4, double Umu4, double Ut4, int LeptonPDG)
{
  double U4[3] = {Ue4, Umu4, Ut4};
  
  for (int i=0;i<6;i++) CM[i] = 0;

  double gl;
  double gr=GR();
  //Calculate the contribution of each mixing matrix element and sum them 
  for (int i=0; i<3; i++) {
    if (LeptonPDG == 11 + 2*i) gl = GL(true);// Take into account the CC contribution 
    else gl=GL(false); 
    CM[0] += 128 * U4[i] * gr * gl; //C1
    CM[1] += 64 * U4[i] * (gl * gl + gr * gr);//C4
    CM[2] += 64 * U4[i] * (gl * gl + gr * gr);//C5
    CM[3] += 0;//C8
    CM[4] += 64 * U4[i] * (gr * gr - gl * gl);//C9
    CM[5] += -64 * U4[i] * (gr * gr - gl * gl);//C10
  }
    
  return;
}

/*Coefficients of Lorentz-Invariant Objects given couplings gL and gR (Dirac HNL)
  See arXiv:2104.05719 (Table 3) for more details
  Valid as long as the HNL is decaying into a neutrino and identical final-state charged leptons */
void evgen::ldm::AnThreeBD::CDNgLgR(double (&CD)[6], double Ue4, double Umu4, double Ut4, int LeptonPDG)
{
  double U4[3]={Ue4, Umu4, Ut4};
  
  for (int i=0;i<6;i++) CD[i] = 0;

  double gl;
  double gr=GR();
  //Calculate the contribution of each mixing matrix element and sum them 
  for (int i=0;i<3;i++) {
    if (LeptonPDG == 11 + 2*i) gl = GL(true);// Take into account the CC contribution 
    else gl = GL(false); 
    
    CD[0] += 64 * U4[i] * gr * gl; //C1
    CD[1] += 64 * U4[i] * (gl * gl);//C4
    CD[2] += 64 * U4[i] * (gr * gr);//C5
    CD[3] += -64 * U4[i] * gr * gl;//C8
    CD[4] += -64 * U4[i] * (gl * gl);//C9
    CD[5] += -64 * U4[i] * (gr * gr);//C10
  }
  
  return;
}


/*Coefficients of Lorentz-Invariant Objects given couplings gL and gR (Dirac bar{HNL} (AntiHNL))
  See arXiv:2104.05719 (Table 3) for more details
  Valid as long as the HNL is decaying into a neutrino and identical final-state charged leptons */
void evgen::ldm::AnThreeBD::CDAgLgR(double (&CD)[6], double Ue4, double Umu4, double Ut4, int LeptonPDG)
{
  double U4[3]={Ue4, Umu4, Ut4};

  for (int i=0;i<6;i++) CD[i] = 0;

  
  double gl;
  double gr=GR();
  //Calculate the contribution of each mixing matrix element and sum them 
  for (int i=0;i<3;i++) {
    if (LeptonPDG == 11 + 2*i) gl = GL(true);// Take into account the CC contribution 
    else gl=GL(false); 
    
    CD[0] += 64 * U4[i] * gr * gl; //C1
    CD[1] += 64 * U4[i] * (gr * gr);//C4
    CD[2] += 64 * U4[i] * (gl * gl);//C5
    CD[3] += 64 * U4[i] * gr * gl;//C8
    CD[4] += 64 * U4[i] * (gr * gr);//C9
    CD[5] += 64 * U4[i] * (gl * gl);//C10
  }
  
  return;
}

/* Calculate Lorentz-Invariants K1, K4, K5, K8, K9, K10 (the only ones required for N->l_alphal_alphanu decay) given the kinematical quantities (z_{ll}, z_{\nu m}, \cos\theta_{ll}, and \gamma_{ll}),
   as well as the masses (N, daughter charged leptons) and the N polarization [-1, 1]
   See arXiv:2104.05719 for more details */
void evgen::ldm::AnThreeBD::KinDep(double (&K)[6], double Pol, double z2_num, double z2_ll, double ct_ll, double gam_ll, double m_HNL, double m_p, double m_m)
{
  //Get the invariant masses from z2_ll and z2_num @LuisPelegrina
  double m2_ll = z2_ll * m_HNL * m_HNL; 
  double m2_num = z2_num * m_HNL * m_HNL;
  
  K[0] = 0.5 * m_m * m_p * (m_HNL * m_HNL - m2_ll); //K1
  K[1] = 0.25 * (m2_num - m_m * m_m) * (m_HNL * m_HNL + m_p * m_p - m2_num);//K4
  K[2] = 0.25 * (m2_ll + m2_num - m_p * m_p) * (m_HNL * m_HNL + m_m * m_m - m2_ll - m2_num);//K5

  //Calculate some needed variables for the spin dependent Lorentz-Invariants @LuisPelegrina
  double E_m = (m2_ll + m2_num - m_p * m_p) / (2.0 * m_HNL);
  double p_m = sqrt(E_m * E_m - m_m * m_m);
  double ct_num = (E_m * (m_HNL * m_HNL - m2_ll) - m_HNL * (m2_num - m_m * m_m)) / ((m_HNL * m_HNL - m2_ll) * p_m);
  double st_num = sqrt(1.0 - ct_num * ct_num);
  double st_ll = sqrt(1.0 - ct_ll * ct_ll);
  
  K[3] = 0.5 * Pol * m_m * m_p * (m_HNL * m_HNL - m2_ll) * ct_ll;//K8
  K[4] = 0.5 * Pol * m_HNL * (m2_num - m_m * m_m) * (p_m * (cos(gam_ll) * st_ll * st_num - ct_ll * ct_num) - (m_HNL * m_HNL - m2_ll) / (2 * m_HNL) * ct_ll);//K9
  K[5] = 0.5 * Pol * m_HNL * p_m * (m_HNL * m_HNL + m_m * m_m - m2_num - m2_ll) * (ct_ll * ct_num - st_num * cos(gam_ll) * st_ll);//K10
  
  return;
}

/*Check whether or not z2_num is allowed by calculating its maximum and minimum with z2_ll
  This expresion was derived from  https://pdg.lbl.gov/2018/reviews/rpp2018-rev-kinematics.pdf Eq. (47.23a) and (47.23b) */
bool evgen::ldm::AnThreeBD::IsDalitzAllowed(double z2_ll, double s, double d,double z2_num)
{
  double st = 2.0 * (1.0 - z2_ll) * sqrt((z2_ll - d*d) * (z2_ll - s*s)) / (4.0 * z2_ll);
  double bt = (d*d * z2_ll + 2.0 * s * d + z2_ll * (2.0 - 2.0 * z2_ll + s*s)) / (4.0 * z2_ll);
  
  if(((bt - st) >  z2_num) | (z2_num > (bt + st))) return false;
   else return true;
  
}

/*Obtain the matrix-element-squared of the N decay
  Requires the N/charged-lepton masses, gL and gR, N polarization, thekinematical quantities (z_{ll}, z_{\nu m}, \cos\theta_{ll}, and \gamma_{ll}),
  and whether the HNL it is Dirac or Majorana
  Coefficients of Lorentz-Invariants are also required
  See arXiv:2104.05719 for more details */
double evgen::ldm::AnThreeBD::MSqDM(double (&C)[6], double z2_num, double z2_ll,  double ct_ll, double gam_ll, double m_HNL, double m_p, double m_m, double Pol)
{
  //Check if the kinematical quantities given are allowed inside the Dalitz
  bool DalitzAllowed = false;
  double s = (m_m + m_p) / m_HNL;
  double d = (m_m - m_p) / m_HNL;
  DalitzAllowed = evgen::ldm::AnThreeBD::IsDalitzAllowed(z2_ll, s, d, z2_num);
  if (!DalitzAllowed) return 0;

  //Calculate the  Lorentz-Invariants K[i]
  double K[6]={0, 0, 0, 0, 0, 0};
  evgen::ldm::AnThreeBD::KinDep(K, Pol, z2_num, z2_ll, ct_ll, gam_ll, m_HNL, m_p, m_m);

  //Multiply the Lorentz-Invariants by its Coefficients and sum them
  double MSq = 0; 
  for (int i=0; i<6; i++) MSq = MSq + C[i] * K[i];

  //return the matrix-element-squared
  return MSq;
}


/*Return the total maximum matrix-element-squared for a given M_HNL,Lepton_PDG, Polarization and mixing matrix element
  Currently optimized for electrons and muons for m_HNl<388 MeV
  The values returned by this function were calculated using 10e7 iterations of the monte Carlo Algorith used when no value is found. If the theory regarding anysotropic decays is changed this values must be calculated again
  This function can be improved, but the maximun is not trivial to estimate, as the allowed Dalitz region restricts the use of some algorithms */
double evgen::ldm::AnThreeBD::MaxMSqDM(double m_HNL, int  LeptonPDG, double Ue4, double Umu4, double Ut4, double m_p, double m_m, double Pol, double (&C)[6], bool Majorana)
{
  double Ce4 = 0;
  double Cmu4 = 0;
  double Ct4 = 0;

  /*If the pair of leptons are electrons or muons and m_HNl<388 MeV dont calculate the maximum
    And access a previously calculated value to save computation time */
  if (Majorana) {
    if(LeptonPDG == 13) { 
      if(m_HNL <= 0.22) {
	Cmu4 = 0.0024;
	Ce4 = 0.00023;
	Ct4 = 0.00023;
      }else if(m_HNL <= 0.24) {
	Cmu4 = 0.0081;
	Ce4 = 0.00083;
	Ct4 = 0.00083;
      }else if(m_HNL <= 0.3) {
	Cmu4 = 0.04;
	Ce4 = 0.0045;
	Ct4 = 0.0045;
      }else if(m_HNL <= 0.389) {
	Cmu4 = 0.135;
	Ce4 = 0.018;
	Ct4 = 0.018;
      }
      
    }else if(LeptonPDG == 11) {
      if(m_HNL <= 0.025) {
	Cmu4 = 4.6 * pow(10,-7);
	Ce4 = 2.9 * pow(10,-6);
	Ct4 = 4.6 * pow(10,-7);
      }else if(m_HNL <= 0.07) {
	Cmu4 = 2.8 * pow(10,-5);
	Ce4 = 0.0002;
	Ct4 = 2.8 * pow(10,-5);
      }else if(m_HNL <= 0.12) {
	Cmu4 = 0.00025;
	Ce4 = 0.0016;
	Ct4 = 0.00025;
      }else if(m_HNL <= 0.185) {
	Cmu4 = 0.0015;
	Ce4 = 0.009;
	Ct4 = 0.0015;
      }else if(m_HNL <= 0.26) {
	Cmu4 = 0.006;
	Ce4 = 0.034;
	Ct4 = 0.006;
      }else if(m_HNL <= 0.390) {
	Cmu4 = 0.027;
	Ce4 = 0.17;
	Ct4 = 0.027;
      }
    }
    
  }else if(!Majorana) {
    if(LeptonPDG == 13) {
	if(m_HNL <= 0.22) {
	  Cmu4 = 0.00142;
	  Ce4 = 0.00023;
	  Ct4 = 0.00023;
	}else if(m_HNL <= 0.24) {
	  Cmu4 = 0.005;
	  Ce4 = 0.00083;
	  Ct4 = 0.00083;
	}else if(m_HNL <= 0.3) {
	  Cmu4 = 0.022;
	  Ce4 = 0.0045;
	  Ct4 = 0.0045;
	}else if(m_HNL <= 0.389) {
	  Cmu4 = 0.077;
	  Ce4 = 0.018;
	  Ct4 = 0.018;
	}
	
    }else if(LeptonPDG == 11) {
      if(m_HNL <= 0.025) {
	Cmu4 = 4.1 * pow(10,-7);
	Ce4 = 1.8 * pow(10,-6);
	Ct4 = 4.1 * pow(10,-7);
      }else if(m_HNL <= 0.07) {
	Cmu4 = 2.6 * pow(10,-5);
	Ce4 = 0.00013;
	Ct4 = 2.6 * pow(10,-5);
      }else if(m_HNL <= 0.12) {
	Cmu4 = 0.00023;
	Ce4 = 0.0011;
	Ct4 = 0.00023;
      }else if(m_HNL <= 0.185) {
	Cmu4 = 0.0013;
	Ce4 = 0.006;
	Ct4 = 0.0013;
      }else if(m_HNL <= 0.26) {
	Cmu4 = 0.005;
	Ce4 = 0.023;
	Ct4 = 0.005;
      }else if(m_HNL <= 0.390) {
	Cmu4 = 0.025;
	Ce4 = 0.11;
	Ct4 = 0.025;
      }
    }
  }

  /*If the HNL mass is greater than 388 MeV or the pair of leptons produce are tau calculate the maximum
    Searching for it randomly, pretty slow aproach, can be improved */
  if((Ce4 == 0)|(Cmu4 == 0)|(Ct4 == 0)) {
    double z2_ll,z2_num,ct_ll,gam_ll;
    //Use 10e6 as a good aproximation of the maximum, if more precision is required the number can be increased
    int N_iter = pow(10,6);
    double MSqMax = 0;
    TRandom3 *rand = new TRandom3(0);
    for (int i=0;i<N_iter;i++) {	  
      double MSq = 0;
      do {
	//Generate a position in the Dalitz Plot randomly
	double r1 = rand->Rndm();
	double r2 = rand->Rndm();
	double r3 = rand->Rndm();
	double r4 = rand->Rndm();
	      
	z2_ll = ((m_HNL * m_HNL - (m_m + m_p) * (m_m + m_p)) * r1 + (m_m + m_p) * (m_m + m_p)) / (m_HNL * m_HNL);
	z2_num = (((m_HNL - m_p) * (m_HNL - m_p) - m_m * m_m) * r2 + m_m * m_m) / (m_HNL * m_HNL);
	ct_ll = r3 * 2 - 1;
	gam_ll = 2 * TMath::Pi() * r4;
	MSq = MSqDM(C, z2_num, z2_ll, ct_ll, gam_ll, m_HNL, m_p, m_m, Pol);
      }while(MSq == 0);
      
      //If the generated MSq is greater than the previus maximum save it as a new maximum    
      if (MSq > MSqMax) MSqMax = MSq;
    }
    
    //Let the user know that the maximum is not being read for memory and is being calculated inestead
    std::cout << "The maximum calculation is slow" <<std::endl;
    
    return MSqMax;
  }
  
  return Ce4 * Ue4 + Cmu4 * Umu4 + Ct4 * Ut4;
}

/*Given the parameters of the final-state (invariant masses and angles)
  determine the rest-frame four-vectors of the outgoing neutrino and charged-lepton pair
  Parameters:
  z2ll = m_{\ell\ell}^2/m_N^2: reduced invariant mass of the charged lepton pair
  z2num = m_{\nu m}^2/m_N^2: reduced invariant mass of the neutrino/negatively-charged-lepton
  ctll: cosine of the angle between the charged-lepton-pair (sum of the four-vectors) and the z-axis
  gamll: rotation angle about the direction of the charged-lepton pair
  masses: [mN, mm, mp] the HNL and daughter charged-lepton masses (mm: negatively-charged, mp: positively-charged) */
void evgen::ldm::AnThreeBD::RF4vecs(TLorentzVector &pnu, TLorentzVector &pm, TLorentzVector &pp, double z2_num, double z2_ll, double ct_ll, double gam_ll, double m_HNL, double m_p, double m_m)
{
  double PnuRF[4], PmRF[4], PpRF[4];
  
  double m2_ll = z2_ll * m_HNL * m_HNL;
  double m2_num = z2_num * m_HNL * m_HNL;
  
  TRandom3 *rand1 = new TRandom3(0);
  double r1 = rand1->Rndm();

  //Calculate some useful kinematical cuantities
  double E_m = (m2_ll + m2_num - m_p * m_p) / (2.0 * m_HNL);
  double E_nu = 0.5 * m_HNL * (1.0 - z2_ll);
  
  double p_m = sqrt(E_m * E_m - m_m * m_m);
  double ct_num = (E_m * (m_HNL * m_HNL - m2_ll) - m_HNL * (m2_num - m_m * m_m)) / ((m_HNL * m_HNL - m2_ll) * p_m);
  double st_num = sqrt(1.0 - ct_num * ct_num);
  double st_ll = sqrt(1.0 - ct_ll * ct_ll);

  //Select the phi angle randomly as it doesnt affect the anisotropies
  double phi = r1 * 2.0 * TMath::Pi();
  phi = 1;
    
  PnuRF[0] = E_nu;
  PnuRF[1] = -E_nu * st_ll * sin(phi);
  PnuRF[2] = -E_nu * st_ll * cos(phi);
  PnuRF[3] = -E_nu * ct_ll;//Neutrino four momentum
  pnu.SetXYZT(PnuRF[1], PnuRF[2], PnuRF[3], PnuRF[0]);
  
  PmRF[0] = E_m;
  PmRF[1] = p_m * (st_num * cos(phi) * sin(gam_ll) - sin(phi) * (st_num * cos(gam_ll) * ct_ll + ct_num * st_ll));
  PmRF[2] = p_m * (-cos(phi) * (st_num * cos(gam_ll) * ct_ll + ct_num * st_ll) - sin(phi) * st_num * sin(gam_ll));
  PmRF[3] = p_m * (-ct_num * ct_ll + st_num * cos(gam_ll) * st_ll);//Negatively charged lepton momentum
  pm.SetXYZT(PmRF[1], PmRF[2], PmRF[3], PmRF[0]);
  
  PpRF[0] = m_HNL - PnuRF[0] - PmRF[0];
  PpRF[1] = -PnuRF[1] - PmRF[1];
  PpRF[2] = -PnuRF[2] - PmRF[2];
  PpRF[3] = -PnuRF[3] - PmRF[3];//positively charged lepton momentum
  pp.SetXYZT(PpRF[1], PpRF[2], PpRF[3], PpRF[0]);
      
  return;
 }


/* Generate a rest-frame HNL anisotropic decay into electron/positron pairs or muon/antimuon pairs (tau not implemented)
  Quantities needed:
  masses: [mN, mm, mp] the HNL and daughter charged-lepton masses (mm: negatively-charged, mp: positively-charged)
  Ue4, Umu4 and Ut4: Mixing matrix element squared
  Lepton PDG: The Pdg of the lepton pair 13=Muon, 11=Electron
  Majorana: Wether the HNL is a Majorana Particle or not
  AntiHNL: Wether the HNL is a AntiHNL or a HNL (Only relevant in Dirac case)
  Pol: Polarization of the HNL
  See arXiv:2104.05719 for more details */
void evgen::ldm::AnThreeBD::AnisotropicThreeBodyDist(TLorentzVector &pnu, TLorentzVector &pm, TLorentzVector& pp, double m_HNL, double Ue4, double Umu4, double Ut4, int LeptonPDG, bool Majorana, bool AntiHNL, double Pol)
{
  //By default select the electron as the lepton pair
  double m_p = evgen::ldm::Constants::Instance().elec_mass;
  double m_m = evgen::ldm::Constants::Instance().elec_mass;

  TRandom3 *rand = new TRandom3(0);

  //Select wether the final leptons are electrons or muons
  if(LeptonPDG == 13) {
    m_p = evgen::ldm::Constants::Instance().muon_mass;
    m_m = evgen::ldm::Constants::Instance().muon_mass;
  } else if(LeptonPDG == 11) {
    m_p = evgen::ldm::Constants::Instance().elec_mass;
    m_m = evgen::ldm::Constants::Instance().elec_mass;
  }


  /*Get the coefficients of Lorentz-Invariant Objects for Majorana or Dirac HNL
    See arXiv:2104.05719 for more details */
  double C[6]={0, 0, 0, 0, 0, 0};
  if (Majorana) evgen::ldm::AnThreeBD::CMgLgR(C, Ue4, Umu4, Ut4, LeptonPDG);
  else {
    if (AntiHNL) evgen::ldm::AnThreeBD::CDAgLgR(C, Ue4, Umu4, Ut4, LeptonPDG);
    else evgen::ldm::AnThreeBD::CDNgLgR(C, Ue4, Umu4, Ut4, LeptonPDG);
  }
  
  double z2_num = 0;
  double z2_ll = 0;
  double ct_ll = 0;
  double gam_ll = 0;

  //Get the Maximum of the Matrix-element squared and use rejection sampling to get the anisotropic final states distribution
  double MSqMax = evgen::ldm::AnThreeBD::MaxMSqDM(m_HNL, LeptonPDG, Ue4, Umu4, Ut4, m_p, m_m, Pol, C, Majorana);
  
  double MSq = 0;
  bool IsW = false;
  do {
    do {
      /*Calculate randomly kinematicals quantities of the decay
            z2ll = m_{\ell\ell}^2/m_N^2: reduced invariant mass of the charged lepton pair
            z2num = m_{\nu m}^2/m_N^2: reduced invariant mass of the neutrino/negatively-charged-lepton
            ctll: cosine of the angle between the charged-lepton-pair (sum of the four-vectors) and the z-axis
            gamll: rotation angle about the direction of the charged-lepton pair
      */
      double r1 = rand->Rndm();
      double r2 = rand->Rndm();
      double r3 = rand->Rndm();
      double r4 = rand->Rndm();
      
      z2_ll = ((m_HNL * m_HNL - (m_m + m_p) * (m_m + m_p)) * r1 + (m_m + m_p) * (m_m + m_p))/(m_HNL * m_HNL);
      z2_num = (((m_HNL - m_p) * (m_HNL - m_p) - m_m * m_m) * r2 + m_m * m_m) / (m_HNL * m_HNL);
      ct_ll = r3 * 2 - 1;
      gam_ll = 2 * TMath::Pi() * r4;
      MSq = evgen::ldm::AnThreeBD::MSqDM(C, z2_num, z2_ll, ct_ll, gam_ll, m_HNL, m_p, m_m, Pol);
    }while(MSq==0);//Check if the randomly selected kinematical quantities  where allowed in the Dalitz
    double r5 = rand->Rndm();
    if(r5 < MSq/MSqMax) IsW = true;
  }while(!IsW); //Use rejection-sampling to select each final state given its probability to happen
  
  //Change from the generated kinematical quantities to 4 Vectors that can be used by the generator
  evgen::ldm::AnThreeBD::RF4vecs(pnu, pm, pp, z2_num, z2_ll, ct_ll, gam_ll, m_HNL, m_p, m_m);
  
  return;
}   
