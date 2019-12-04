///
/// \file    EarthModel.cxx
/// \brief   Provide interface to different Earth models
/// \author  messier@indiana.edu
/// \version $Id: EarthModel.cxx,v 1.2 2012-07-01 22:21:59 gsdavies Exp $
///
#include "OscLib/EarthModel.h"
#include "CAFAna/Core/MathUtil.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <list>
using namespace osc;

static const int kPREM   = 0;
static const int kSTACEY = 1;

EarthModel::EarthModel(const char* which, double tol) 
{
  if      (std::string("PREM")  ==which) this->InitPREM();
  else if (std::string("STACEY")==which) this->InitStacey();
  // else if (std::string("ModelX")==which) this->InitModelX();
  else {
    std::cerr << __FILE__ << ":" << __LINE__ 
	      << " Model '" << which << "' is not supported." 
	      << " Currently only PREM is supported." << std::endl;
    abort();
  }
  this->MakeLayers(tol);
}

//......................................................................
///
/// Return the electron number density at radius r
///
/// \param r - radius in units of km
/// \returns electron densty in units of moles/cm^3
///
double EarthModel::Ne(double r) 
{
  return this->ZoverA(r)*this->Density(r);
}

//......................................................................
///
/// Earth density in g/cm^3 as a function of radius in km
///
double EarthModel::Density(double r) 
{
  switch (fModelID) {
  case kPREM:   return this->DensityPREM(r);
  case kSTACEY: return this->DensityStacey(r);
  default:      abort();
  }
  return 0.0;
}

//......................................................................
///
/// Initialize the Earth model to "PREM"
///
/// [*] Dziewonski and Anderson, Physics of the Earth and Planetary
/// Interiors, 25 (1981) 297-356. Table I.
///
void EarthModel::InitPREM()
{
  fModelID      = kPREM;
  fRouterCore   = 3480.0;
  fRearth       = 6371.0;

  fRregion.push_back(0.0);
  fRregion.push_back(1221.5); // Inner core
  fRregion.push_back(3480.0); // Outer core
  fRregion.push_back(3630.0); // D''
  fRregion.push_back(5701.0); // Lower mantle
  fRregion.push_back(6151.0); // Transition zone
  fRregion.push_back(6291.0); // Low velocity zone
  fRregion.push_back(6346.6); // LID
  fRregion.push_back(6368.0); // Crust
  fRregion.push_back(6371.0); // Ocean
}

//......................................................................
///
/// The PREM Earth density profile in units of g/cm^3 as a function of
/// radius r (km)
///
/// [*] Dziewonski and Anderson, Physics of the Earth and Planetary
/// Interiors, 25 (1981) 297-356. Table I.
///
/// \param r - radius in km
/// \returns density in g/cm^3
///
double EarthModel::DensityPREM(double r) 
{
  double x = r/fRearth;
  
  // The regions are checked in order from largest to smallest
  if (r>1221.5&&r<=3480.0) return 12.5815+x*(-1.2638+x*(-3.6426+x*(-5.5281)));
  if (r>3480.0&&r<=5701.0) return  7.9565+x*(-6.4761+x*(5.5283+x*(-3.0807)));
  if (r<=1221.5)           return 13.0885-8.8381*x*x;
  if (r>5771.0&&r<=5971.0) return 11.2494-8.0298*x;
  if (r>5971.0&&r<=6151.0) return  7.1089-3.8045*x;
  if (r>6151.0&&r<=6291.0) return  2.6910+0.6924*x;
  if (r>5701.0&&r<=5771.0) return  5.3197-1.4836*x;
  if (r>6291.0&&r<=6346.6) return  2.6910+0.6924*x;
  if (r>6346.6&&r<=6356.0) return  2.90;
  if (r>6356.0&&r<=6368.0) return  2.60;
  if (r>6368.0&&r<=6371.0) return  1.020;
  return 0.0;
}

//......................................................................
///
/// Initialize the Earth model to Stacey:
///
/// [*] F.D. Stacey, Physics of the Earth (Wiley, New York, 1969)
///
/// \param r - radius in km
/// \returns density in g/cm^3
///
void EarthModel::InitStacey()
{
  fModelID      = kSTACEY;
  fRouterCore   = 6371.0-2878.0;
  fRearth       = 6371.0;

  fRregion.push_back(0.0);
  fRregion.push_back(6371.0-5121.0); // Inner core
  fRregion.push_back(6371.0-2878.0); // Outer core
  fRregion.push_back(6371.0-650.0);  // Lower mantle
  fRregion.push_back(6371.0-350.0);  // Middle mantle
  fRregion.push_back(6371.0-15.0);   // Upper mantle
  fRregion.push_back(6371.0);        // Crust
}

//......................................................................
///
/// Earth density profile according to
///
/// [*] F.D. Stacey, Physics of the Earth (Wiley, New York, 1969)
///
double EarthModel::DensityStacey(double r)
{
  static double crustD[2] = {       0.0,
				    15.0};
  static double upperMantleD[6] = { 15.0,
				    60.0,
				    100.0,
				    200.0,
				    300.0, 
				    350.0};
  static double middleMantleD[6] = {350.0,
				    400.0,
				    413.0,
				    500.0,
				    600.0, 
				    650.0};
  static double lowerMantleD[14] = {650.0,
				    800.0,
				    984.0,
				    1000.0,
				    1200.0, 
				    1400.0,
				    1600.0,
				    1800.0,
				    2000.0,
				    2200.0, 
				    2400.0,
				    2600.0,
				    2800.0,
				    2878.0};
  static double outerCoreD[14] = {  2878.0,
				    3000.0,
				    3200.0,
				    3400.0,
				    3600.0, 
				    3800.0,
				    4000.0,
				    4200.0,
				    4400.0,
				    4600.0, 
				    4800.0,
				    4982.0,
				    5000.0,
				    5121.0};
  static double innerCoreD[8] ={    5121.0,
				    5200.0,
				    5400.0,
				    5600.0,
				    5800.0, 
                                    6000.0,
				    6200.0,
				    6371.0};
  static double crustRho[2] =       { 2.840,
				      2.840};
  static double upperMantleRho[6] = { 3.313,
				      3.332,
				      3.348,
				      3.387,
				      3.424,
				      3.441};
  static double middleMantleRho[6] = {3.700,
				      3.775,
				      3.795,
				      3.925,
				      4.075, 
				      4.150};
  static double lowerMantleRho[14] = {4.200,
				      4.380,
				      4.529,
				      4.538,
				      4.655,
				      4.768,
				      4.877,
				      4.983,
				      5.087,
				      5.188,
				      5.288,
				      5.387,
				      5.487,
				      5.527};
  static double outerCoreRho[14] = {  9.927,
				      10.121,
				      10.421,
				      10.697,
				      10.948, 
				      11.176,
				      11.383,
				      11.570,
				      11.737,
				      11.887, 
				      12.017,
				      12.121,
				      12.130,
				      12.197};
  static double innerCoreRho[7] = {   12.229,
				      12.301,
				      12.360,
				      12.405,
				      12.437, 
				      12.455,
				      12.460};
  // ============ END DENSITY TABLE ============
  
  if (r>fRearth) return 0.0;
  double       d = fRearth-r; // Table uses depth
  int i;
  
  if (d>=crustD[0] && d<crustD[1]) {
    for (i=0; i<1; ++i) {
      if (d>=crustD[i] && d<crustD[i+1]) {
	return crustRho[i];
      }
    }
  }
  if (d>=upperMantleD[0] && d<upperMantleD[5]) {
    for (i=0; i<5; ++i) {
      if (d>=upperMantleD[i] && d<upperMantleD[i+1]) {
	return upperMantleRho[i];
      }
    }
  }
  if (d>=middleMantleD[0] && d<middleMantleD[5]) {
    for (i=0; i<5; ++i) {
      if (d>=middleMantleD[i] && d<middleMantleD[i+1]) {
	return middleMantleRho[i];
      }
    }
  }
  if (d>=lowerMantleD[0] && d<lowerMantleD[13]) {
    for (i=0; i<13; ++i) {
      if (d>=lowerMantleD[i] && d<lowerMantleD[i+1]) {
	return lowerMantleRho[i];
      }
    }
  }
  if (d>=outerCoreD[0] && d<outerCoreD[13]) {
    for (i=0; i<13; ++i) {
      if (d>=outerCoreD[i] && d<outerCoreD[i+1]) {
	return outerCoreRho[i];
      }
    }
  }
  if (d>=innerCoreD[0] && d<innerCoreD[7]) {
    for (i=0; i<7; ++i) {
      if (d>=innerCoreD[i] && d<innerCoreD[i+1]) {
	return innerCoreRho[i];
      }
    }
  }
  return 0.0;
}

//......................................................................
///
/// Mean ratio of atomic number to atomic mass for different segments
/// of the Earth. These are taken from:
///
/// [*] Bahcall and Krastev, Phys. Rev. C56, p.2839 (1997)
///
/// \param r - radius in km
/// \returns Z/A at this radius
///
double EarthModel::ZoverA(double r) 
{
  if (r<fRouterCore) return 0.468;
  return 0.497;
}

//......................................................................

double EarthModel::AveNe(double r1, double r2, int nsample)
{
  double sum = 0.0;
  double r;
  for (int i=0; i<nsample; ++i) {
    r = r1 + (r2-r1)*((float)i+0.5)/(float)nsample;
    sum += this->Ne(r);
  }
  return sum/(float)nsample;
}

//......................................................................
///
/// Make layers of constant density keeping the electron density
/// profile to within a fractional tolerance
///
void EarthModel::MakeLayers(double tol) 
{
  unsigned int nsample = 20;
  for (unsigned int i=0; i<fRregion.size()-1; ++i) {
    // Add layers within each region until all layers are within
    // tolerance
    double rRegionLo = fRregion[i];
    double rRegionHi = fRregion[i+1];
    double r1 = rRegionLo;
    double r2 = rRegionHi;
    double r;
    double ne;
    while (1) {
      double ave  = this->AveNe(r1, r2, nsample);
      bool   isok = true;
      for (unsigned int j=0; j<nsample; ++j) {
	r = r1+(r2-r1)*((float)j+0.5)/(float)nsample;
	ne = this->Ne(r);
	if (fabs(ne-ave)/ne>tol) { isok = false; break; }
      }
      if (isok) {
	// Layer is within tolerance - accept it
	fRlo.push_back(r1);
	fRhi.push_back(r2);
	fNe. push_back(ave);
	if (r2>=rRegionHi) {
	  // Finished subdividing this region
	  break;
	}
	else {
	  // Ready next iteration of subdivision
	  r1 = r2;
	  r2 = rRegionHi;
	}
      } // if (isok) ...
      else {
	// Layer is not within tolerance - reduce its size
	r2 -= 1.0; // Step back some distance (km) and try again
	if (r2<=r1) r2 = r1+0.5;
      }
    } // making layers within region
  } // loop on regions
}

//......................................................................

void EarthModel::GetLayers(std::vector<double>& rlo,
			   std::vector<double>& rhi,
			   std::vector<double>& ne)
{
  rlo.resize(fRlo.size());
  rhi.resize(fRlo.size());
  ne. resize(fRlo.size());
  for (unsigned int i=0; i<fRlo.size(); ++i) {
    rlo[i] = fRlo[i];
    rhi[i] = fRhi[i];
    ne[i]  = fNe[i];
  }
}

//......................................................................
///
/// Find the list of distances a neutrino travels through each layer
/// of the Earth model
///
/// \param anu  - Altitude of neutrino production (km above sea level)
/// \param cosQ - cosine of neutrnio zenith angle (+1=down, -1=up)
/// \param adet - Altitude of detector (km above sea level)
///
void EarthModel::LineProfile(double anu,
			     double cosQ,
			     double adet,
			     std::list<double>& Ls,
			     std::list<double>& Ns)
{
  // A typical number density for the crust
  static const double crustNe = 2.6*0.5;
  
  // Locations are relative to a cicular Earth centered at (0,0)
  double rdet; // Radial location of detector (km)
  double rnu;  // Radial location of neutrino (km)
  double Lnu;  // Total flight distance of neutrino (km)
  double x1;   // Location of detector (km)
  double y1;   // Location of detector (km)
  double x2;   // Location of neutrino production (km)
  double y2;   // Location of neutrino production (km)
  
  rdet = fRearth+adet;
  rnu  = fRearth+anu;
  Lnu  = rdet*(sqrt(util::sqr(rnu/rdet)+cosQ*cosQ-1.0)-cosQ);

  x2 = 0.0;
  y2 = rdet;
  x1 = x2+sqrt(1.0-cosQ*cosQ)*Lnu;
  y1 = y2+cosQ*Lnu;
  
  unsigned int i;
  int n1, n2;
  double L;
  double xina,  yina;
  double xouta, youta;
  double xinb,  yinb;
  double xoutb, youtb;
  Ls.clear();
  Ns.clear();

  //
  // For down-going neutrinos, only consider matter in the top layer
  //
  if (cosQ>=0.0) {
    if (adet>0.0) {
      Ls.push_back(Lnu);
      Ns.push_back(0.0);
    }
    else {
      n1 = this->IntersectLineAndCircle(x1, y1, x2, y2, fRearth,
					&xouta, &youta, &xoutb, &youtb);
      L = util::pythag(x1-xoutb, y1-youtb);
      Ls.push_back(L);
      Ns.push_back(0.0);
      
      L = util::pythag(x2-xoutb, y2-youtb);
      Ls.push_back(L);
      Ns.push_back(crustNe); // Typical rock density
    }
    return;
  }
  
  //
  // For up-going neutrinos, find the intersections of the neutrino's
  // path with all layers of the Earth model.
  //
  for (i=0; i<fRhi.size()-1; ++i) {
    n1 = this->IntersectLineAndCircle(x1, y1, x2, y2, fRhi[i],
				      &xouta, &youta, &xoutb, &youtb);
    
    n2 = this->IntersectLineAndCircle(x1, y1, x2, y2, fRlo[i],
				      &xina, &yina, &xinb, &yinb);
    //
    // Create segments for each layer we intersect
    //
    if (n1==2 && n2<2) {
      //
      // Intersect outer radius, but not inner radius: 
      // Only one one segment to create.
      //
      L = util::pythag(xouta-xoutb, youta-youtb);
      
      Ls.push_back(L);
      Ns.push_back(fNe[i]);
    }
    else if (n1==2 && n2==2) {
      //
      // Path intersects both inner and outer radii: Two sections to
      // create, one on the way in, one on the way out.
      //
      L = util::pythag(xouta-xina, youta-yina);
      Ls.push_front(L);
      Ns.push_front(fNe[i]);
      
      L = util::pythag(xoutb-xinb, youtb-yinb);
      Ls.push_back(L);
      Ns.push_back(fNe[i]);
    }
  } // Loop on all but last layer
  
  //
  // Handle the last layer specially
  //
  i = fRhi.size()-1;
  n1 = this->IntersectLineAndCircle(x1, y1, x2, y2, fRhi[i],
				    &xouta, &youta, &xoutb, &youtb);
  n2 = this->IntersectLineAndCircle(x1, y1, x2, y2, fRlo[i],
				    &xina, &yina, &xinb, &yinb);
  if (n1==2 && n2<2) {
    //
    // Intersect outer radius, but not inner radius: 
    // One one segment to create
    //
    if (youtb>y2) {
      //
      // Don't over shoot the detector location
      //
      L = util::pythag(xouta-x2, youta-y2);
    }
    else {
      L = util::pythag(xouta-xoutb, youta-youtb);
    }
    Ls.push_back(L);
    Ns.push_back(fNe[i]);
  }
  else if (n1==2 && n2==2) {
    //
    // Path intersects both inner and outer radii: Two sections to
    // create, one on the way in, one on the way out.
    //
    L = util::pythag(xouta-xina, youta-yina);
    Ls.push_front(L);
    Ns.push_front(fNe[i]);
    
    if (youtb>y2) {
      //
      // Don't over shoot the detector location
      //
      L = util::pythag(x2-xinb, y2-yinb);
    }
    else {
      L = util::pythag(xoutb-xinb, youtb-yinb);
    }
    Ls.push_back(L);
    Ns.push_back(fNe[i]);
  }  
  //
  // Add segment to get neutrino from production to Earth's surface
  //
  L = util::pythag(xouta-x1, youta-y1);
  Ls.push_front(L);
  Ns.push_front(0.0); // Density of air~=0
  
  //
  // If required, add segment to get from sea level up to the
  // detector's position
  //
  if (youtb<y2) {
    L = util::pythag(xoutb-x2, youtb-y2);
    Ls.push_back(L);
    Ns.push_back(crustNe); // Density of of standard rock (mole/cm^3)
  }

  //
  // The following is a useful piece of debugging - check that we've
  // accounted for all pieces of the neutrino flight distance
  //
  if (false) {
    std::list<double>::iterator itr(Ls.begin());
    std::list<double>::iterator itrEnd(Ls.end());
    double ltot = 0.0;
    for (; itr!=itrEnd; ++itr) {
      ltot += (*itr);
    }
    std::cout << ltot << " " << Lnu << std::endl;
    if (fabs(ltot-Lnu)>1e-6) abort();
  }
}

//......................................................................

int EarthModel::IntersectLineAndCircle(double  x1, double  y1,
				       double  x2, double  y2,
				       double  r,
				       double* xa, double* ya,
				       double* xb, double* yb) 
{
  double dx    = x2-x1;
  double dy    = y2-y1;
  double drsqr = dx*dx + dy*dy;
  double D     = x1*y2 - x2*y1;
  double arg   = r*r*drsqr - D*D;
  if (arg<0.0) return 0;

  double sgndy = 1.0;
  if (dy<0.0) sgndy = -1.0;

  *xa = ( D*dy - sgndy * dx * sqrt(arg) ) / drsqr;
  *ya = (-D*dx - fabs(dy)   * sqrt(arg) ) / drsqr;
  *xb = ( D*dy + sgndy * dx * sqrt(arg) ) / drsqr;
  *yb = (-D*dx + fabs(dy)   * sqrt(arg) ) / drsqr;

  if (arg==0.0) return 1;
  return 2;
}

////////////////////////////////////////////////////////////////////////
