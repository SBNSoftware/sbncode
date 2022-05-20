#ifndef TPCGEOALG_H_SEEN
#define TPCGEOALG_H_SEEN


/////////////////////////////////////////////////////////////////////
// TPCGeoAlg.h
//
// Wrapper for easy access to TPC geometry
// Written by T Brooks (tbrooks@fnal.gov) for SBND, May 2019
// Ported to icaruscode by C Hilgenberg (chilgenb@fnal.gov), Oct 2021
/////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"

// c++
#include <vector>


namespace sbn{

  class TPCGeoAlg {
  public:

    TPCGeoAlg();

    ~TPCGeoAlg();

    // Getters
    double MinX() const;
    double MinY() const;
    double MinZ() const;
    double MaxX() const;
    double MaxY() const;
    double MaxZ() const;
    double CpaWidth() const;

    // Functions for applying fiducial volume cuts to total volume
    bool InFiducial(geo::Point_t point, double fiducial);
    bool InFiducial(geo::Point_t point, double fiducial, double fiducialTop);
    bool InFiducial(geo::Point_t point, double minXCut, double minYCut, double minZCut, 
                    double maxXCut, double maxYCut, double maxZCut);
    
    // Is point inside given TPC
    bool InsideTPC(geo::Point_t point, const geo::TPCGeo& tpc, double buffer=0.);

    // Determine which TPC a collection of hits is detected in (-1 if multiple)
    int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits);
    // Determine the drift direction for a collection of hits (-1, 0 or 1 assuming drift in X)
    int DriftDirectionFromHits(std::vector<art::Ptr<recob::Hit>> hits);
    // Work out the drift limits for a collection of hits
    std::pair<double, double> XLimitsFromHits(std::vector<art::Ptr<recob::Hit>> hits);

    double MinDistToWall(geo::Point_t point) const;

    // Determine if a true particle is ever inside the TPC volume
    bool InVolume(const simb::MCParticle& particle);
    // Determine if a true particle is contained inside the TPC volume
    bool IsContained(const simb::MCParticle& particle);
    // Determine if a true particle enters the TPC volume
    bool EntersVolume(const simb::MCParticle& particle);
    // Determine if a true particle crosses the TPC volume
    bool CrossesVolume(const simb::MCParticle& particle);
    // Determine if a true particle crosses either APA
    bool CrossesApa(const simb::MCParticle& particle);

    std::pair<TVector3, TVector3> CrossingPoints(const simb::MCParticle& particle);
    double TpcLength(const simb::MCParticle& particle);

  private:

    double fMinX;
    double fMinY;
    double fMinZ;
    double fMaxX;
    double fMaxY;
    double fMaxZ;
    double fCpaWidth;

    geo::GeometryCore const* fGeometryService;

  };

}

#endif


