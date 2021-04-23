////////////////////////////////////////////////////////////////////////
// \file PhotonVisibilityService.h
//
// \brief Service to report opdet visibility to different points in
//         the system
//
// \author bjpjones@mit.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef PHOTONVISIBILITYSERVICE_H
#define PHOTONVISIBILITYSERVICE_H


//#include "art/Framework/Services/Registry/ActivityRegistry.h"
//#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "PhotonLibrary.h"
#include "PhotonVoxels.h"
#include <cassert>

///General LArSoft Utilities
namespace phot{
  
  class PhotonVisibilityService {
  public:
   
    //PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    //PhotonVisibilityService(const std::string library="uboone_photon_library_v5.root");
    PhotonVisibilityService(const std::string library="PhotonLibrary-20180801.root");
    
    //void reconfigure(fhicl::ParameterSet const& p);
    
    //double GetQuenchingFactor(double dQdx);
    
    //double DistanceToOpDet(                 double* xyz, unsigned int OpDet );
    //double SolidAngleFactor(                double* xyz, unsigned int OpDet );
    float GetVisibility(double* xyz, unsigned int OpChannel ) const;
    float GetVisibility(double x, double y, double z, unsigned int OpChannel ) const;

    std::vector<std::vector<float> > GetVisibilityXY(double z) const;
    std::vector<std::vector<float> > GetVisibilityYZ(double x) const;
    std::vector<std::vector<float> > GetVisibilityZX(double y) const;

    float  Fraction2AbsoluteX(float frac) const;
    float  Fraction2AbsoluteY(float frac) const;
    float  Fraction2AbsoluteZ(float frac) const;

    inline float  GetMaxX() const { return fXmax; }
    inline float  GetMaxY() const { return fYmax; }
    inline float  GetMaxZ() const { return fZmax; }
    inline float  GetMinX() const { return fXmin; }
    inline float  GetMinY() const { return fYmin; }
    inline float  GetMinZ() const { return fZmin; }
    inline int    GetNX() const { return fNx; }
    inline int    GetNY() const { return fNy; }
    inline int    GetNZ() const { return fNz; }
    inline size_t GetNOpChannels() const { return fNOpDetChannels; }

    const std::vector<float>* GetAllVisibilities( double* xyz ) const;

    inline const std::vector<std::vector<float> >& GetLibraryData() const
    { if(!fTheLibrary) LoadLibrary(); return fTheLibrary->GetData(); }
    
    void LoadLibrary() const;
    void StoreLibrary();
    
    
    void StoreLightProd(    int  VoxID,  double  N );
    void RetrieveLightProd( int& VoxID,  double& N ) const;
    
    void SetLibraryEntry(   int VoxID, int OpChannel, float N);
    float GetLibraryEntry( int VoxID, int OpChannel) const;
    const std::vector<float>* GetLibraryEntries( int VoxID ) const;

    
    bool IsBuildJob() const { return fLibraryBuildJob; }
    bool UseParameterization() const {return fParameterization;}

    sim::PhotonVoxelDef GetVoxelDef() const {return fVoxelDef; }
    int NOpChannels() const { return fNOpDetChannels; }

    const std::string& GetLibraryFilename() { return fLibraryFile; }; // Allows one to check loaded filename

    static PhotonVisibilityService& GetME(std::string filename="PhotonLibrary-20180801.root")
    {
      // argument allows user to choose file loaded
      // if relative, searches for file in designated folder
      // if absolute, searches for file with given path
      // see LoadLibrary()
      if(!_me) _me = new PhotonVisibilityService(filename);
      return *_me;
    }
    
  private:

    static PhotonVisibilityService* _me;
    
    int    fCurrentVoxel;
    double fCurrentValue;
    
    float  fXmin, fXmax;
    float  fYmin, fYmax;
    float  fZmin, fZmax;
    int    fNx, fNy, fNz;
    size_t fNOpDetChannels;
    bool fUseCryoBoundary;
    
    bool                 fLibraryBuildJob;
    bool                 fDoNotLoadLibrary;
    bool                 fParameterization;
    std::string          fLibraryFile;      
    mutable PhotonLibrary* fTheLibrary;
    sim::PhotonVoxelDef  fVoxelDef;
    
    
  }; // class PhotonVisibilityService
} //namespace utils
//DECLARE_ART_SERVICE(phot::PhotonVisibilityService, LEGACY)
#endif // UTIL_DETECTOR_PROPERTIES_H
