#ifndef SBNCODE_SINGLEPHOTONANALYSIS_HELPER_MATH_H
#define SBNCODE_SINGLEPHOTONANALYSIS_HELPER_MATH_H

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/Shower.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Utilities/ensurePointer.h"

#include <numeric>
#include <vector>

namespace single_photon
{
  //templates
  // sort indices in descending order
  template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v){

      std::vector<size_t> idx(v.size());
      std::iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

      return idx;
    }

  // sort indices such that elements in v are in ascending order
  template <typename T>
     std::vector<size_t> sort_indexes_rev(const std::vector<T> &v) {

      std::vector<size_t> idx(v.size());
      std::iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

      return idx;
    }

  // check if two vectors have same elements (regardless of the order), and arrange their elements in order
  template<typename T>
    bool marks_compare_vec_nonsense(std::vector<T>& v1, std::vector<T>& v2)
    {
      std::sort(v1.begin(), v1.end());
      std::sort(v2.begin(), v2.end());
      return v1 == v2;
    }

//inline functions
  inline double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
    double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
    return wire;
  }


  /* dot product of wire_dir and shower direction vectors */
  inline double getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir){
    return wire_dir.Dot(shower_dir);
  }
  inline double degToRad(double deg){ return deg * M_PI/180; }
  inline double radToDeg(double rad){ return rad * 180/M_PI; }



  //-----------------HELPER FUNCTIONS -----------
  ////line between x1 and x2, point x0;
  double dist_line_point( std::vector<double>&X1, std::vector<double>& X2, std::vector<double>& point);

  double impact_paramater_shr(double x, double y, double z, art::Ptr<recob::Shower> & shr);

  // invariant mass of a particle that decays to two showers
  double  implied_invar_mass(double vx, double vy, double vz, art::Ptr<recob::Shower> & s1, double E1,  art::Ptr<recob::Shower> &s2, double E2);

  // invariant mass of two showers, calculated directly from shower directions
  double  invar_mass(art::Ptr<recob::Shower> & s1, double E1,  art::Ptr<recob::Shower> &s2, double E2);

  double getMedian(std::vector<double> thisvector);

  /* returns (generally) best median dEdx of all 3
   * planes, usually plane 2  */
  double getAmalgamateddEdx(
      double angle_wrt_plane0, 
      double angle_wrt_plane1, 
      double angle_wrt_plane2, 
      double median_plane0, 
      double median_plane1, 
      double median_plane2, 
      int plane0_nhits, 
      int plane1_nhits, 
      int plane2_nhits);

  /* returns the number of hits on the plane picked by function getAmalgamateddEdx */
  int getAmalgamateddEdxNHits(
      double amalgamateddEdx, 
      double median_plane0, 
      double median_plane1, 
      double median_plane2,
      int plane0_nhits, 
      int plane1_nhits, 
      int plane2_nhits);


  /**
   *@brief Calculates the four corners of a rectangle of given length and width around a cluster given the start point and axis direction
   *@param cluster_start - the start position of a cluster in CM
   *@param cluster_axis - calculated from the cluster end minus the cluster start
   *@param width - typically ~1cm
   *@param length - typically a few cm
   *
   * */
  std::vector<std::vector<double>> buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length);

}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_HELPER_MATH_H
