#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/fiducial_volume.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

namespace single_photon
{

void setTPCGeom(para_all& paras){

  //copy these from CAFMaker/CAFMaker_module.c
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
    }
    paras.fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: paras.fTPCVolumes) {
    paras.s_tpc_active_XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    paras.s_tpc_active_YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    paras.s_tpc_active_ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    paras.s_tpc_active_XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    paras.s_tpc_active_YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    paras.s_tpc_active_ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();
    if(g_is_verbose){
      std::cout<<""<<__FUNCTION__<<" || Active TPC info: X:("<<paras.s_tpc_active_XMin<<","<<paras.s_tpc_active_XMax<<")";
      std::cout<<" Y:("<<paras.s_tpc_active_YMin<<","<<paras.s_tpc_active_YMax<<")";
      std::cout<<" Z:("<<paras.s_tpc_active_ZMin<<","<<paras.s_tpc_active_ZMax<<")"<<std::endl;
    }
  }

}


/* inside TPC or not? */
bool isInTPCActive(std::vector<double> & vec, para_all& paras){
  if( vec.size() != 3){
    throw cet::exception("single_photon") << " The coordinate dimension is not 3!";
  }

  bool is_x = (vec[0] > paras.s_tpc_active_XMin && vec[0]< paras.s_tpc_active_XMax );
  bool is_y = (vec[1] > paras.s_tpc_active_YMin && vec[1]< paras.s_tpc_active_YMax );
  bool is_z = (vec[2] > paras.s_tpc_active_ZMin && vec[2]< paras.s_tpc_active_ZMax );
  bool inside = is_x&&is_y&&is_z;

  return inside;
}


/* returns minimum distance to the TPC active boundary; returns -999 if the point is not in TPC active volume */
double distToTPCActive(std::vector<double>&vec, para_all& paras){
  if(!isInTPCActive(vec, paras)) return -999;
  double min_x = std::min( fabs(vec[0] - paras.s_tpc_active_XMin) ,  fabs(vec[0] - paras.s_tpc_active_XMax));
  double min_y = std::min( fabs(vec[1] - paras.s_tpc_active_YMin) ,  fabs(vec[1] - paras.s_tpc_active_YMax));
  double min_z = std::min( fabs(vec[2] - paras.s_tpc_active_ZMin) ,  fabs(vec[2] - paras.s_tpc_active_ZMax));

  return ( (min_x<min_y) ? std::min(min_x,min_z)  : std::min(min_y,min_z)    );
}


/* returns minimum distance to the TPCActive boundary around the Cathode Plane Assemble; returns -999 if the point is not in TPC active volume */
double distToCPA(std::vector<double>&vec, para_all& paras){
  if(!isInTPCActive(vec, paras)) return -999;
  double dx = std::min( fabs(vec[0] - (-0.45)) ,  fabs(vec[0] - 0.45));

  return dx;
}




//    int isInSCB(std::vector<double> & vec){
//        return isInSCB(0.0,vec);
//    }

int distToSCB(double & dist, std::vector<double> &vec, para_all& paras){
  //CHECK!
  dist = distToTPCActive( vec, paras);
  int is_it_in = 0;
  if(isInTPCActive( vec, paras)) is_it_in++;
  return is_it_in;
  //NOT USE SCB YET, bring it back later!
  //
  //        //this one returns the distance to the boundary
  //        bool ans = false;
  //        double dist_yx = 999999;
  //
  //        int iseg = 0;
  //        if (!isInTPCActive(vec)){
  //            dist=-1;
  //            return 0; // is it in active volume?
  //        }
  //        double cut = 0.0;
  //
  //        TVector3 pt(&vec[0]);
  //
  //        Int_t z_idx = (Int_t)pt.Z()/100;  // Int_t: signed integer 4 bytes
  //
  //        Int_t z_idx_YX = z_idx;//YX-view effective z index, it is only different for z > 10m area, where we want to appliy 9m<z<10m YX boundary, still need to keep the original z_idx bc it's needed in ZX-view
  //        if (z_idx_YX==10) z_idx_YX-=1;
  //        double tbi = -10000.; //means "to be initialized"
  //
  //
  //        double ptX[6] = {0.+cut, tbi, m_SCB_YX_TOP_x2_array-cut, m_SCB_YX_BOT_x2_array-cut, tbi, 0.+cut};
  //        double ptY[6] = {m_SCB_YX_TOP_y1_array-cut, m_SCB_YX_TOP_y1_array-cut, tbi, tbi, m_SCB_YX_BOT_y1_array+cut, m_SCB_YX_BOT_y1_array+cut};
  //
  //        TGeoPolygon *polyXY = new TGeoPolygon(6);
  //
  //        ptX[1] = m_SCB_YX_TOP_x1_array[z_idx_YX+1];
  //        ptX[4] = m_SCB_YX_BOT_x1_array[z_idx_YX+1];
  //        ptY[2] = m_SCB_YX_TOP_y2_array[z_idx_YX+1];
  //        ptY[3] = m_SCB_YX_BOT_y2_array[z_idx_YX+1];
  //
  //        polyXY->SetXY(ptX,ptY);
  //        polyXY->FinishPolygon();
  //        double testpt[2] = {pt.X(), pt.Y()};
  //
  //        //cout << "is testpt ("<< pt.X()<<", "<<pt.Y()<<") contrained? "<<  polyXY->Contains(testpt)<<endl;
  //        //cout << "area ? " << polyXY->Area()<< endl;
  //
  //        //polyXY->Draw();    
  //
  //        Bool_t XY_contain = polyXY->Contains(testpt);
  //        dist_yx = polyXY->Safety(testpt, iseg); // Compute minimum distance from testpt to any segment.
  //
  //        if(0<z_idx && z_idx<10){
  //            double up_z = pt.Z()-s_tpc_active_z_low; // gonna bet, if it's middle enough to be up_z or down_z is smaller than the safefy in this z regime (1m,10m), it is safe to set up_z = z-0, down_z=1036.8-z 
  //            double down_z = s_tpc_active_z_high-pt.Z();
  //            double min_z =  std::min(up_z,down_z);
  //
  //            XY_contain ?  dist = std::min(dist_yx,min_z) : dist = -1  ;
  //
  //            delete polyXY;
  //            return  (XY_contain ? 1 : 0 );
  //        }
  //
  //        //up or down
  //        double top_y = s_tpc_active_y_high-pt.Y();
  //        double bottom_y = pt.Y()+s_tpc_active_y_high;
  //        double min_y = std::min(top_y, bottom_y);
  //
  //
  //        // if z_idx==0 or z_idx==10, they need xz view,  
  //        /// ZX view has Y dependence: Y sub-range from -116 to 116cm per 24cm
  //        Int_t y_idx = (pt.Y()+116.)/24;
  //        if (pt.Y()<-116. && pt.Y()>-116.5) y_idx = 0; //just the 0.5cm 
  //        if(y_idx<0 || y_idx>9) {
  //            dist = -1; 
  //            delete polyXY;  
  //            return 0;
  //        }
  //
  //        Bool_t ZX_contain = false;
  //        double arbit_out = 55555;
  //
  //        double d_dist = -1;
  //
  //        //upstream
  //        if(z_idx==0){
  //            double ptX_Up[5] = {0.+cut, m_SCB_ZX_Up_x1_array, m_SCB_ZX_Up_x2_array-cut, m_SCB_ZX_Up_x2_array-cut, 0+cut};
  //            double ptZ_Up[5] = {0.+cut,0.+cut,m_SCB_ZX_Up_z2_array, arbit_out, arbit_out};
  //
  //            TGeoPolygon *polyXZ_Up = new TGeoPolygon(5);
  //            polyXZ_Up->SetXY(ptX_Up,ptZ_Up);
  //            polyXZ_Up->FinishPolygon();
  //
  //            double testpt_Up[2] = {pt.X(), pt.Z()};
  //            ZX_contain = polyXZ_Up->Contains(testpt_Up);
  //            d_dist = polyXZ_Up->Safety(testpt_Up,iseg);
  //            delete polyXZ_Up;
  //   
  //        }
  //        if (z_idx==10){
  //            double ptX_Dw[5] = {0.+cut,m_SCB_ZX_Dw_x2_array-cut, m_SCB_ZX_Dw_x2_array-cut, tbi, 0.+cut};
  //            double ptZ_Dw[5] = {-arbit_out,-arbit_out,tbi,m_SCB_ZX_Dw_z1_array-cut, m_SCB_ZX_Dw_z1_array-cut};
  //
  //            ptX_Dw[3] = m_SCB_ZX_Dw_x1_array[y_idx+1];
  //            ptZ_Dw[2] = m_SCB_ZX_Dw_z2_array[y_idx+1];
  //
  //            TGeoPolygon *polyXZ_Dw = new TGeoPolygon(5);
  //            polyXZ_Dw->SetXY(ptX_Dw,ptZ_Dw);
  //            polyXZ_Dw->FinishPolygon();
  //
  //            double testpt_Dw[2] = {pt.X(), pt.Z()};
  //            ZX_contain = polyXZ_Dw->Contains(testpt_Dw);
  //            d_dist = polyXZ_Dw->Safety(testpt_Dw,iseg);
  //            delete polyXZ_Dw;
  //        }
  //
  //        delete polyXY;
  //
  //        ans = XY_contain && ZX_contain;
  //        ans ? dist = std::min(d_dist,min_y) :  dist=-1;
  //
  //        return (ans ? 1 : 0);
}


}
