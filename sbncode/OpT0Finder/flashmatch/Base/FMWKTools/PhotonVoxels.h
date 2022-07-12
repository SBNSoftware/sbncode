#ifndef PhotonVoxels_h 
#define PhotonVoxels_h 1

#include "TVector3.h"
#include <iostream>

namespace sim {


  class PhotonVoxel{
  public:
    PhotonVoxel(double xMin, 
		double xMax, 
		double yMin, 
		double yMax, 
		double zMin, 
		double zMax, 
		int N = 0) ;
    PhotonVoxel();

  private:
    double xVoxelMin;
    double xVoxelMax;
    double yVoxelMin;
    double yVoxelMax;
    double zVoxelMin;
    double zVoxelMax;

    int NPhotons;

  public:

    TVector3 GetLowerCorner() const;
    TVector3 GetUpperCorner() const;
    TVector3 GetCenter()      const;

  };


  class PhotonVoxelDef
  {
  public:
    PhotonVoxelDef(double xMin, 
		   double xMax, 
		   int xN, 
		   double yMin, 
		   double yMax, 
		   int yN, 
		   double zMin, 
		   double zMax, 
		   int z);
    PhotonVoxelDef();
    
  private:
    TVector3 fLowerCorner;
    TVector3 fUpperCorner;
    int      fxSteps;
    int      fySteps;
    int      fzSteps;

#ifndef __GCCXML__
  public:
    TVector3 GetRegionUpperCorner() const;
    TVector3 GetRegionLowerCorner() const;
    TVector3 GetSteps() const;


    TVector3 GetVoxelSize()const;

    int GetNVoxels() const;

    int GetVoxelID(TVector3) const;
    int GetVoxelID(double*)  const;
    int GetVoxelID(double x, double y, double z)  const;
    bool IsLegalVoxelID(int) const;

    PhotonVoxel      GetPhotonVoxel(int ID) const;
    std::vector<int> GetVoxelCoords(int ID) const;
    //PhotonVoxel      GetContainingVoxel(TVector3) const{};

    bool operator==(const PhotonVoxelDef &rhs) const;
    bool operator!=(const PhotonVoxelDef &rhs) const 
      { return ! ((*this)==rhs); }

#endif

  };
}

#endif
