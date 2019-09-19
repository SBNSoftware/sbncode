#include <iostream>
#include "AutoVetoTools.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace ana {
  namespace AutoVetoAnalysis {

      bool IsAV(const core::ProviderManager &manager, const TVector3 &point){
          const geo::GeometryCore* geom = manager.GetGeometryProvider();
          geo::CryostatGeo const& cryo0 = geom->Cryostat(0);
          geo::CryostatGeo const& cryo1 = geom->Cryostat(1);
          geo::TPCGeo const& tpc00 = cryo0.TPC(0);
          geo::TPCGeo const& tpc01 = cryo0.TPC(1);
          geo::TPCGeo const& tpc10 = cryo1.TPC(0);
          geo::TPCGeo const& tpc11 = cryo1.TPC(1);

          if(tpc00.ContainsPosition(point) ||
            tpc01.ContainsPosition(point) ||
            tpc10.ContainsPosition(point) ||
            tpc11.ContainsPosition(point) )
              return true; 

          return false;
      }

      bool IsFV(const core::ProviderManager &manager, const TVector3 &point, const TVector3 &fid1, const TVector3 &fid2) {
      
          const geo::GeometryCore* geom = manager.GetGeometryProvider();
          geo::CryostatGeo const& cryo0 = geom->Cryostat(0);
          geo::CryostatGeo const& cryo1 = geom->Cryostat(1);
          geo::TPCGeo const& tpc00 = cryo0.TPC(0);
          geo::TPCGeo const& tpc01 = cryo0.TPC(1);
          geo::TPCGeo const& tpc10 = cryo1.TPC(0);
          geo::TPCGeo const& tpc11 = cryo1.TPC(1);

          if(tpc00.ContainsPosition(point)
             && tpc00.InFiducialX(point.X(),fid1.X(),fid2.X())
             && tpc00.InFiducialY(point.Y(),fid1.Y(),fid2.Y())
             && tpc00.InFiducialZ(point.Z(),fid1.Z(),fid2.Z()) )
              return true;
          if(tpc01.ContainsPosition(point)
             /*tpc01.InFiducialX(point.X(),fid1.X(),fid2.X())
             */&& tpc01.InFiducialY(point.Y(),fid1.Y(),fid2.Y())
             && tpc01.InFiducialZ(point.Z(),fid1.Z(),fid2.Z()) )
              return true;
          if(tpc10.ContainsPosition(point)
             /*tpc10.InFiducialX(point.X(),fid1.X(),fid2.X())
             */&& tpc10.InFiducialY(point.Y(),fid1.Y(),fid2.Y())
             && tpc10.InFiducialZ(point.Z(),fid1.Z(),fid2.Z()) )
              return true;
          if(tpc11.ContainsPosition(point)
             && tpc11.InFiducialX(point.X(),fid2.X(),fid1.X())
             && tpc11.InFiducialY(point.Y(),fid1.Y(),fid2.Y())
             && tpc11.InFiducialZ(point.Z(),fid1.Z(),fid2.Z()) )
              return true;

          return false;
      }

      int MacToADReg(int mac) {

          if(mac>=107 && mac<=190) return 30; //top
          if(mac>=191 && mac<=204) return 31; //rim west
          if(mac>=205 && mac<=218) return 32; //rim east
          if(mac>=219 && mac<=224) return 33; //rim south
          if(mac>=225 && mac<=230) return 34; //rim north
          if(            mac<=12 ) return 40; //west side, south stack
          if(mac>=13  && mac<=24 ) return 41; //west side, center stack
          if(mac>=25  && mac<=36 ) return 42; //west side, north stack
          if(mac>=37  && mac<=48 ) return 43; //east side, south stack
          if(mac>=49  && mac<=60 ) return 44; //east side, center stack
          if(mac>=61  && mac<=72 ) return 45; //east side, north stack
          if(mac>=73  && mac<=84 ) return 46; //south
          if(mac>=85  && mac<=92 ) return 47; //north
          if(mac>=93 && mac<=106) return 50; //bottom

          return 0;
     } 

  }  // namespace AutoVetoAnalysis
}  // namespace ana

