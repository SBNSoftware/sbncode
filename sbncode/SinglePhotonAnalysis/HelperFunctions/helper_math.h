#ifndef SBNCODE_SINGLEPHOTONANALYSIS_HELPER_MATH_H
#define SBNCODE_SINGLEPHOTONANALYSIS_HELPER_MATH_H

#include <numeric>

namespace single_photon
{
  //-----------------HELPER FUNCTIONS -----------
  ////line between x1 and x2, point x0;
  inline double dist_line_point( std::vector<double>&X1, std::vector<double>& X2, std::vector<double>& point){
    double x1 =X1.at(0);
    double y1 =X1.at(1);
    double z1 =X1.at(2);

    double x2 =X2.at(0);
    double y2 =X2.at(1);
    double z2 =X2.at(2);

    double x0 =point.at(0);
    double y0 =point.at(1);
    double z0 =point.at(2);

    double x10 = x1-x0;
    double y10 = y1-y0;
    double z10 = z1-z0;

    double x21 = x2-x1;
    double y21 = y2-y1;
    double z21 = z2-z1;

    double t = -(x10*x21+y10*y21+z10*z21)/fabs(x21*x21+y21*y21+z21*z21 );

    // right, but can be simplified
    double d2 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0)+(z2-z1)*(z1-z0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));

    return sqrt(d2);
  }

  // minimum distance between point and dead wire
  //    double distanceToNearestDeadWire(int plane, double Ypoint, double Zpoint,
  //                const geo::GeometryCore * geom,
  //                std::vector<std::pair<int,int>> & bad_channel_list_fixed_mcc9 ){
  //
  //             double min_dist = 999999;
  //            
  //             for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
  //                   int channel = bad_channel_list_fixed_mcc9[i].first;                                       
  //                   int is_ok = bad_channel_list_fixed_mcc9[i].second;       
  //                        if(is_ok>1)continue;
  //
  //                   auto wireids = geom->ChannelToWire(channel); //type of wireids: IDs of all the connected wires to the channel
  //                   auto result = geom->WireEndPoints(wireids[0]);
  //                        
  //                    //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 
  ////                    std::cout<<wireids[0].Plane<<" "<<result.start().X()<<std::endl;
  //                    if(plane != (int)wireids[0].Plane) continue;
  //
  //        // for the dead wires on the same plane 
  //                    std::vector<double> start = {0.0,result.start().Y(),result.start().Z()};
  //                    std::vector<double> end = {0.0,result.end().Y(),result.end().Z()};
  //                    std::vector<double> point = {0.0,Ypoint,Zpoint};
  //                    double dist = dist_line_point(start,end,point);
  //                    min_dist = std::min(dist,min_dist);
  //             }
  //                 
  //    return min_dist;
  //
  //    }

  //--------------- end of copying from bad_channel_matching.h

  inline double impact_paramater_shr(double x, double y, double z, art::Ptr<recob::Shower> & shr){

    std::vector<double> vert = {x,y,z}; 
    std::vector<double> start = {shr->ShowerStart().X(), shr->ShowerStart().Y(),shr->ShowerStart().Z()};
    std::vector<double> abit = {shr->ShowerStart().X() + shr->Direction().X(),  shr->ShowerStart().Y()+shr->Direction().Y(),  shr->ShowerStart().Z()+shr->Direction().Z()};

    return dist_line_point(start, abit, vert);

  }

  // invariant mass of a particle that decays to two showers
  inline double  implied_invar_mass(double vx, double vy, double vz, art::Ptr<recob::Shower> & s1, double E1,  art::Ptr<recob::Shower> &s2, double E2){

    double s1x = s1->ShowerStart().X()-vx;
    double s1y = s1->ShowerStart().Y()-vy;
    double s1z = s1->ShowerStart().Z()-vz;
    double norm1  = std::hypot(s1x,s1y,s1z);
    s1x = s1x/norm1; //unit vector pointing to shower start from point (vx, vy,vz)
    s1y = s1y/norm1;
    s1z = s1z/norm1;

    double s2x = s2->ShowerStart().X()-vx;
    double s2y = s2->ShowerStart().Y()-vy;
    double s2z = s2->ShowerStart().Z()-vz;
    double norm2  = std::hypot(s2x,s2y,s2z);
    s2x = s2x/norm2; // unit vector pointing to shower start from point (vx, vy, vz)
    s2y = s2y/norm2;
    s2z = s2z/norm2;

    return sqrt(2.0*E1*E2*(1.0-(s1x*s2x+s1y*s2y+s1z*s2z)));


  }

  // invariant mass of two showers, calculated directly from shower directions
  inline double  invar_mass(art::Ptr<recob::Shower> & s1, double E1,  art::Ptr<recob::Shower> &s2, double E2){

    double s1x = s1->Direction().X();
    double s1y = s1->Direction().Y();
    double s1z = s1->Direction().Z();

    double s2x = s2->Direction().X();
    double s2y = s2->Direction().Y();
    double s2z = s2->Direction().Z();

    return sqrt(2.0*E1*E2*(1.0-(s1x*s2x+s1y*s2y+s1z*s2z)));

  }



  // sort indices in descending order
  template <typename T>
    inline std::vector<size_t> sort_indexes(const std::vector<T> &v) {

      std::vector<size_t> idx(v.size());
      std::iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

      return idx;
    }

  // sort indices such that elements in v are in ascending order
  template <typename T>
    inline  std::vector<size_t> sort_indexes_rev(const std::vector<T> &v) {

      std::vector<size_t> idx(v.size());
      std::iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

      return idx;
    }


  // check if two vectors have same elements (regardless of the order), and arrange their elements in order
  template<typename T>
    inline bool marks_compare_vec_nonsense(std::vector<T>& v1, std::vector<T>& v2)
    {
      std::sort(v1.begin(), v1.end());
      std::sort(v2.begin(), v2.end());
      return v1 == v2;
    }



  inline double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
    double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
    return wire;
  }

  /* unit vector orthogonal to the  wire direction of plane -- usually named as wire_dir */
  inline TVector3 getWireVec(int plane){
    TVector3 wire_dir;
    if (plane == 0){
      wire_dir = {0., -sqrt(3) / 2., 1 / 2.};
    } else if (plane == 1){
      wire_dir = {0., sqrt(3) / 2., 1 / 2.};
    } else if (plane == 2) {
      wire_dir = {0., 0., 1.};
    }
    return wire_dir;
  }

  /* dot product of wire_dir and shower direction vectors */
  inline double getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir){
    return wire_dir.Dot(shower_dir);
  }

  /* returns angles between wire direction of plane  and shower_dir) 
   *  shower_dir needs to be unit vector */
  inline double getAnglewrtWires(TVector3 shower_dir,int plane){

    TVector3 wire_dir = getWireVec(plane);
    double cos_theta =  getCoswrtWires(shower_dir, wire_dir);

    double theta = acos(cos_theta);
    // return abs(theta);
    return abs(M_PI/2 - theta);

  }

  inline  double degToRad(double deg){ return deg * M_PI/180; }
  inline double radToDeg(double rad){ return rad * 180/M_PI; }

  inline double getMedian(std::vector<double> thisvector){
    size_t len = thisvector.size();
    if(len < 1) return NAN;

    std::sort(thisvector.begin(), thisvector.end());
    if(len % 2 != 0){//even - return average of two at median
      return 0.5*(thisvector[len/2]+thisvector[len/2+1]);
    }else{//odd - return the median
      return thisvector[len/2];
    }
  }

  /* returns (generally) best median dEdx of all 3
   * planes, usually plane 2  */
  inline double getAmalgamateddEdx(
      double angle_wrt_plane0, 
      double angle_wrt_plane1, 
      double angle_wrt_plane2, 
      double median_plane0, 
      double median_plane1, 
      double median_plane2, 
      int plane0_nhits, 
      int plane1_nhits, 
      int plane2_nhits){
    //if the shower is within 10 degrees of the wires on plane 2, consider planes 1 and 0
    if(angle_wrt_plane2< degToRad(10)){
      //if it's too close to the wires on either of the planes, then stick with plane 2
      if (angle_wrt_plane1> degToRad(20)|| angle_wrt_plane0>degToRad(20) ){
        //but if it's outside of the range on plane 1, choose that
        if(angle_wrt_plane1> angle_wrt_plane0){
          return median_plane1;
        } else{
          return median_plane0;
        }
      }
    }
    if (plane2_nhits< 2){
      if (plane1_nhits >=2 ){
        return median_plane1;
      } else if (plane0_nhits >=2 ){
        return median_plane0;
      }
    }
    return median_plane2;
  }

  /* returns the number of hits on the plane picked by function getAmalgamateddEdx */
  inline int getAmalgamateddEdxNHits(
      double amalgamateddEdx, 
      double median_plane0, 
      double median_plane1, 
      double median_plane2,
      int plane0_nhits, 
      int plane1_nhits, 
      int plane2_nhits){
    if (amalgamateddEdx == median_plane0){
      return plane0_nhits;
    }
    if (amalgamateddEdx == median_plane1){
      return plane1_nhits;
    }
    if (amalgamateddEdx == median_plane2){
      return plane2_nhits;
    }
    return -999;
  }


  /**
   *@brief Calculates the four corners of a rectangle of given length and width around a cluster given the start point and axis direction
   *@param cluster_start - the start position of a cluster in CM
   *@param cluster_axis - calculated from the cluster end minus the cluster start
   *@param width - typically ~1cm
   *@param length - typically a few cm
   *
   * */
  inline std::vector<std::vector<double>> buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length){
    std::vector<std::vector<double>> corners;

    //get the axis perpedicular to the cluster axis
    double perp_axis[2] = {-cluster_axis[1], cluster_axis[0]};

    //create a vector for each corner of the rectangle on the plane
    //c1 = bottom left corner
    std::vector<double> c1 = {cluster_start[0] + perp_axis[0] * width / 2,  cluster_start[1] + perp_axis[1] * width / 2};
    //c2 = top left corner
    std::vector<double> c2 = {c1[0] + cluster_axis[0] * length, c1[1] + cluster_axis[1] * length};
    //c3 = bottom right corner
    std::vector<double> c3 = {cluster_start[0] - perp_axis[0] * width / 2, cluster_start[1] - perp_axis[1] * width / 2};
    //c4 = top right corner
    std::vector<double> c4 ={c3[0] + cluster_axis[0] * length, c3[1] + cluster_axis[1] * length}; 

    //save each of the vectors
    corners.push_back(c1);
    corners.push_back(c2);
    corners.push_back(c4);
    corners.push_back(c3);
    return corners;
  }

}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_HELPER_MATH_H
