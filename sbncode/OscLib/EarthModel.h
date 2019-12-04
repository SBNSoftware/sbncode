////////////////////////////////////////////////////////////////////////
/// \file    EarthModel.h
/// \brief   Provide and interface to different Earth models
/// \author  messier@indiana.edu
/// \version $Id: EarthModel.h,v 1.1 2012-06-29 02:06:31 bckhouse Exp $
////////////////////////////////////////////////////////////////////////
#ifndef EARTHMODEL_H
#define EARTHMODEL_H
#include <list>
#include <vector>

namespace osc {
  class EarthModel 
  {
  public:
    /// Construct an Earth model
    ///
    /// \param which - Which model to build (eg. PREM, STACEY)
    /// \param tol - Fractional tolerance to use to construct layers.
    ///
    EarthModel(const char* which, double tol);

    /// Return the electron number density at radius r
    ///
    /// \param r - radius in units of km
    /// \returns electron densty in units of moles/cm^3
    ///
    double Ne(double r);

    /// Earth density in g/cm^3 as a function of radius in km
    ///
    /// \param r - radius in km
    /// \returns density in units of g/cm^3
    ///
    double Density(double r);

    /// Ratio of atomic Z to A inside the Earth
    ///
    /// \param r - radius in km
    /// \returns ratio of Z to A
    double ZoverA(double r);
    
    /// Return a list of shells of constant density based on the Earth
    /// model that represent the model to within the specified
    /// tolerance
    ///
    /// \param rlo - list of lower radial edges of layers
    /// \param rhi - list of upper radial edges of layers
    /// \param ne  - electron density (
    ///
    void GetLayers(std::vector<double>& rlo, 
		   std::vector<double>& rhi,
		   std::vector<double>& ne);

    /// Find the list of distances a neutrino travels through each layer
    /// of the Earth model
    ///
    /// \param anu  - Altitude of neutrino production (km above sea level)
    /// \param cosQ - cosine of neutrnio zenith angle (+1=down, -1=up)
    /// \param adet - Altitude of detector (km above sea level)
    ///
    void LineProfile(double prodL, double cosQ, double rdet,
		     std::list<double>& Ls,
		     std::list<double>& Ns);
  private:
    void   InitPREM();
    void   InitStacey();
    double DensityPREM(double r);
    double DensityStacey(double r);
    double AveNe(double r1, double r2, int nsample);
    void   MakeLayers(double tol);
    int    IntersectLineAndCircle(double  x1, double  y1,
				  double  x2, double  y2,
				  double  r,
				  double* xa, double* ya,
				  double* xb, double* yb);
  private:
    int                 fModelID;    ///< Which model has been selected?
    double              fRouterCore; ///< Radius of outer core in km
    double              fRearth;     ///< Radius of Earth in km
    std::vector<double> fRregion;    ///< Region boundaries in km

    std::vector<double> fRlo;    ///< Inner radius of ith layer (km)
    std::vector<double> fRhi;    ///< Outer radius of ith layer (km)
    std::vector<double> fNe;     ///< Electron density of ith layer (mole/cm^3)
  };
}
#endif
////////////////////////////////////////////////////////////////////////
