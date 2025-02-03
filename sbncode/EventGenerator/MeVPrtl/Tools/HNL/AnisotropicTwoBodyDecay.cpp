#include "AnisotropicTwoBodyDecay.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

//Implementation of HNL Two Body Decays Anisotropies (N -> m+l-)
//Valid as long as the HNL is decaying into a charged lepton and a charged pseudoscalar meson
//See arXiv:1905.00284 for more details @LuisPelegrina

/*
Fuction to calculate the momentum of the daugther particles in a 2 body decay
  parent_mass: The mass of the parent particle in the decay
  childA_mass, childB_mass: Mass of the child particles in the decay
*/
double evgen::ldm::AnTwoBD::two_body_momentum(double parent_mass, double childA_mass, double childB_mass) {
  return sqrt(pow(parent_mass * parent_mass + childA_mass * childA_mass - childB_mass * childB_mass, 2) -
              4 * parent_mass * parent_mass * childA_mass * childA_mass) / (2 * parent_mass);
}

//Helper function to calculate  a random 3 dimensional vector
TVector3 evgen::ldm::AnTwoBD::RandomUnitVector(TRandom3 *rand) {
  // In order to pick a random point on a sphere -- pick a random value of _costh_, __not__ theta
  // b.c. d\Omega = d\phi dcos\theta, i.e. d\Omega != d\phi d\theta
  double costheta = rand->Uniform(-1, 1);
  double sintheta = sqrt(1. - costheta * costheta);
  double phi = rand->Uniform(0, 2*M_PI);
  return TVector3(sintheta * cos(phi), sintheta * sin(phi), costheta);
}

//Helper function to calculate the lambda term present in the I_1_plus, I_1_minus e I_1^+ functions
double evgen::ldm::AnTwoBD::lambda(double a, double b, double c) {
  return a*a+b*b+c*c-2*a*b-2*b*c-2*c*a;
}

//Function to calculate the I_1^+ function that characterizes the angular anisotropies. See arXiv:1905.00284 for more details
double evgen::ldm::AnTwoBD::I_1_plus(double x,double y, double cos_theta) {
  double sqrt_lambda = sqrt(evgen::ldm::AnTwoBD::lambda(1,x,y));
  return sqrt_lambda*( (1-x)*(1-x) - y * (1+x) + (x - 1)* sqrt_lambda* cos_theta    )/(4*M_PI);
}

//Function to calculate the I_1^- function that characterizes the angular anisotropies. See arXiv:1905.00284  for more details
double evgen::ldm::AnTwoBD::I_1_minus(double x,double y, double cos_theta) {
  double sqrt_lambda = sqrt(evgen::ldm::AnTwoBD::lambda(1,x,y));
  return sqrt_lambda*( (1-x)*(1-x) - y * (1+x) - (x - 1)* sqrt_lambda* cos_theta    )/(4*M_PI);
}


/*
  Function to apply a rotation to the calculated vector to match the HNL Pol direction. 
  Formulas in arXiv:1905.00284 are defined with an HNL with the spin in the z direction. 
  A rotation must be applied to the decay products to change the z direction to the beam direction (or HNL direction) 
  in order for the results to be consistent, as per the reference arXiv:1905.00284
  p_HNL: HNL momentum vector (3D)
  p_l: Lepton momentum vector (3D)
  p_m: Meson momentum vector (3D)
*/
void evgen::ldm::AnTwoBD::RotateToHNLPolFrame(TVector3 p_HNL, TVector3 &p_l, TVector3 &p_m) {
  //Define the HNL direction
  TVector3 HNL_dir = p_HNL.Unit();
  
  //Define the direction of the polarized HNL in which formulas are calculated (In this case z)
  TVector3 z(0, 0, 1);
  
  // Calculate the axis of rotation (cross product of z and hnl direction)
  TVector3 rotation_axis = z.Cross(HNL_dir);

  // Calculate the angle of rotation 
  double angle = z.Angle(HNL_dir);
  
  //Apply the rotation
  if (rotation_axis.Mag() != 0) {
    p_l.Rotate(angle, rotation_axis);
    p_m.Rotate(angle, rotation_axis);
  } else if(angle == M_PI) {
    p_l = -p_l;
    p_m = -p_m;
  }

  return;
}

/*
  Funtcion to calculate the maximum Integral. This Integral is defined as the weighted sum of I+ and I-. 
  It must be calculated to apply the rejection-acceptance method
*/
double evgen::ldm::AnTwoBD::GetMaxIntegral(int LeptonPDG, double Pol, double eps_m, double eps_l) {
  double max_I;
  if(LeptonPDG > 0) {
    if(Pol > 0) {
      max_I = (1 - Pol) * evgen::ldm::AnTwoBD::I_1_minus(eps_l*eps_l, eps_m*eps_m, -1)/2 + (1 + Pol) * evgen::ldm::AnTwoBD::I_1_plus(eps_l*eps_l, eps_m*eps_m, -1)/2;
    } else {
      max_I = (1 - Pol) * evgen::ldm::AnTwoBD::I_1_minus(eps_l*eps_l, eps_m*eps_m, 1)/2 + (1 + Pol) * evgen::ldm::AnTwoBD::I_1_plus(eps_l*eps_l, eps_m*eps_m, 1)/2;
    }
  } else {
    if(Pol > 0) {
      max_I = (1 - Pol) * evgen::ldm::AnTwoBD::I_1_plus(eps_l*eps_l, eps_m*eps_m, 1)/2 + (1 + Pol) * evgen::ldm::AnTwoBD::I_1_minus(eps_l*eps_l, eps_m*eps_m, 1)/2;
    } else {
      max_I = (1 - Pol) * evgen::ldm::AnTwoBD::I_1_plus(eps_l*eps_l, eps_m*eps_m, -1)/2 + (1 + Pol) * evgen::ldm::AnTwoBD::I_1_minus(eps_l*eps_l, eps_m*eps_m, -1)/2;
    }
  }
  return max_I;
}


/* Generate a rest-frame HNL anisotropic decay into a lepton meson final state
  Quantities needed:
  masses: [m_HNL, m_l, m_m] the HNL and daughter lepton (m_l) and meson (m_m) mass
  Lepton PDG: The Pdg of the lepton 13=Muon, 11=Electron
  Meson PDG: The Pdg of the meson 211=Pion, 321=Kaon
  Pol: Polarization of the HNL
  See arXiv:1905.00284 for more details */
void evgen::ldm::AnTwoBD::AnisotropicTwoBodyDist( TLorentzVector pHNL, TLorentzVector &pl, TLorentzVector &pm, double m_HNL, int LeptonPDG, int MesonPDG, double Pol)
{
  //Get the lepton mass, by default choose the lepton to be a muon
  double m_l = evgen::ldm::Constants::Instance().muon_mass;
  if (abs(LeptonPDG) == 13) {
    m_l = evgen::ldm::Constants::Instance().muon_mass;
  } else if (abs(LeptonPDG) == 11) {
    m_l = evgen::ldm::Constants::Instance().elec_mass;
  }

  //Get the messon mass, by default choose the meson to be a pion
  double m_m = evgen::ldm::Constants::Instance().piplus_mass;
  if (abs(MesonPDG) == 211) {
    m_m = evgen::ldm::Constants::Instance().piplus_mass;
  } else if (abs(MesonPDG) == 321) {
    m_m = evgen::ldm::Constants::Instance().kplus_mass;
  }

  //Initialize the seed to apply rejection-acceptance method. Use 0 to automatically compute the seed via a TUUID object.
  TRandom3 *rand = new TRandom3(0);

  /*
    Declare variables to use inside the loop and for the calculations   
    p_mag: Momentum of the child particles in the two body decay
    eps_m: Quotient of dividing the lepton mass and the HNL mass
    eps_m: Quotient of dividing the meson mass and the HNL mass
    I: value of the integral for each iteration
    max_I: Maximum value of the integral, needed for the acceptance-rejection method
    l_dir: Direction of the lepton in the decay
  */
  double p_mag = evgen::ldm::AnTwoBD::two_body_momentum(m_HNL, m_l, m_m);
  double eps_l = m_l/m_HNL;
  double eps_m = m_m/m_HNL;
  double I;
  double max_I = evgen::ldm::AnTwoBD::GetMaxIntegral(LeptonPDG, Pol, eps_m, eps_l);
  TVector3 l_dir;

  //Use a rejection-acceptance method to select final states following the desired anisotropic distribution
  do {
    //Get a random direction for the lepton in the decay
    l_dir = evgen::ldm::AnTwoBD::RandomUnitVector(rand);
    double cos_theta = l_dir.Z();
    //Calculate the integral value in the direction chosen as a weigthed sum of the + and - polariztaion integral as given in arXiv:1905.00284
    if(LeptonPDG > 0) {
      I =(1 -Pol) * evgen::ldm::AnTwoBD::I_1_minus(eps_l*eps_l, eps_m*eps_m, cos_theta)/2 + (1 + Pol) * evgen::ldm::AnTwoBD::I_1_plus(eps_l*eps_l, eps_m*eps_m, cos_theta)/2;
    } else {
      I = (1 -Pol) * evgen::ldm::AnTwoBD::I_1_plus(eps_l*eps_l, eps_m*eps_m, cos_theta)/2 + (1 + Pol) * evgen::ldm::AnTwoBD::I_1_minus(eps_l*eps_l, eps_m*eps_m, cos_theta)/2;
    }
    
    //If the integral is bigger than the maximum possible Integral give out an error
    if(I > max_I)  {
      std::cerr << "VERY VERY BAD!!!! Integral in Two Body anisotropic N->ml decay is bigger than maximum!!!\n";
      std::cout << "VERY VERY BAD!!!! Integral in Two Body anisotropic N->ml decay is bigger than maximum!!!\n";
    }
  
  //Reject or accept the state
  } while( rand->Uniform(0, 1) > I/max_I);

  //Calculate the daughters three-momentum and prepare the HNL momentum for calculations
  TVector3 pl_RF = TVector3(p_mag * l_dir);
  TVector3 pm_RF = TVector3(-p_mag * l_dir);
  TVector3 pHNL_RF = TVector3(pHNL.X(), pHNL.Y(), pHNL.Z());
  
  //Rotate the leptons to match the HNL Polarization direction frame.
  RotateToHNLPolFrame(pHNL_RF, pl_RF, pm_RF);

  //Set the momentum following the anisotropic distribution
  pl.SetXYZT(pl_RF.X(), pl_RF.Y(), pl_RF.Z(), sqrt(pl_RF.Mag()*pl_RF.Mag() + m_l*m_l));
  pm.SetXYZT(pm_RF.X(), pm_RF.Y(), pm_RF.Z(), sqrt(pm_RF.Mag()*pm_RF.Mag() + m_m*m_m));
  
  return;
}   
