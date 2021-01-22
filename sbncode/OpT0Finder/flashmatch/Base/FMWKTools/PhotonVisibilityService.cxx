////////////////////////////////////////////////////////////////////////
//
//  \file PhotonVisibilityService_service.cc
//
////////////////////////////////////////////////////////////////////////
//
//  Ben Jones, MIT 2012
//
//  This service reports the visibility of a particular point in
//  the detector to each OpDet.  This is used by the fast
//  optical simulation and by track-light association algorithms.
//
//  Visibility is defined as the fraction of isotropically produced
//  photons from a detector voxel which are expected to reach the
//  OpDet in question.
//
//  This information is lookup up from a previousely generated
//  optical library file, whose path is specified to this service.
//
//  Note that it is important that the voxelization schemes match
//  between the library and the service instance for sensible results.
//
//
// Framework includes

// LArSoft includes
#include "PhotonVisibilityService.h"
#include "PhotonLibrary.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
//#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "Geometry/Geometry.h"
//#include "Geometry/CryostatGeo.h"
//#include "Geometry/OpDetGeo.h"
#include <chrono>

using namespace std::chrono;
namespace phot{

  PhotonVisibilityService* PhotonVisibilityService::_me = nullptr;

  //--------------------------------------------------------------------
  /*
  PhotonVisibilityService::PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry &) :
    fCurrentVoxel(0),
    fCurrentValue(0.),
    fXmin(0.),
    fXmax(0.),
    fYmin(0.),
    fYmax(0.),
    fZmin(0.),
    fZmax(0.),
    fNx(0),
    fNy(0),
    fNz(0),
    fUseCryoBoundary(false),
    fLibraryBuildJob(false),
    fDoNotLoadLibrary(false),
    fParameterization(false),
    fTheLibrary(0)
  {
    this->reconfigure(pset);
    mf::LogInfo("PhotonVisibilityService")<<"PhotonVisbilityService initializing"<<std::endl;
  }
  */

  //--------------------------------------------------------------------
  //void PhotonVisibilityService::reconfigure(fhicl::ParameterSet const& p)
  PhotonVisibilityService::PhotonVisibilityService(std::string library) :
    fCurrentVoxel(0),
    fCurrentValue(0.),
    fXmin( -395.0 ),
    fXmax(  -45.0 ),
    fYmin( -215.2 ),
    fYmax(  174.8 ),
    fZmin( -995.0 ),
    fZmax(  965.0 ),
    fNx(70),
    fNy(78),
    fNz(392),
    fNOpDetChannels(180),
    fUseCryoBoundary(true),
    fLibraryBuildJob(false),
    fDoNotLoadLibrary(false),
    fParameterization(false),
    fLibraryFile(library),
    fTheLibrary(nullptr)
  {
    fVoxelDef = sim::PhotonVoxelDef(fXmin, fXmax, fNx, fYmin, fYmax, fNy, fZmin, fZmax, fNz);
    return;
  }

  std::vector<std::vector<float> >
  PhotonVisibilityService::GetVisibilityYZ(double x) const {
    std::vector<std::vector<float> > result(fNy,std::vector<float>(fNz,0.));
    if(x < fXmin || x > fXmax) return result;

    if(fTheLibrary == 0)
      LoadLibrary();

    for(int iy=0; iy<fNy; ++iy) {
      for(int iz=0; iz<fNz; ++iz) {
	int vox_id = iy*fNx + iz * (fNy + fNx);
	double vis_sum = 0.;
	for(auto const& vis_pmt : *(fTheLibrary->GetCounts(vox_id)))
	  vis_sum += ((double)(vis_pmt));
	result[iy][iz] = vis_sum;
      }
    }
    return result;
  }

  std::vector<std::vector<float> >
  PhotonVisibilityService::GetVisibilityZX(double y) const {
    std::vector<std::vector<float> > result(fNz,std::vector<float>(fNx,0.));
    if(y < fYmin || y > fYmax) return result;

    if(fTheLibrary == 0)
      LoadLibrary();

    for(int ix=0; ix<fNx; ++ix) {
      for(int iz=0; iz<fNz; ++iz) {
	int vox_id = ix + iz * (fNy + fNx);
	double vis_sum = 0.;
	for(auto const& vis_pmt : *(fTheLibrary->GetCounts(vox_id)))
	  vis_sum += ((double)(vis_pmt));
	result[iz][ix] = vis_sum;
      }
    }
    return result;
  }

  std::vector<std::vector<float> >
  PhotonVisibilityService::GetVisibilityXY(double z) const {
    std::vector<std::vector<float> > result(fNx,std::vector<float>(fNy,0.));
    if(z < fZmin || z > fZmax) return result;

    if(fTheLibrary == 0)
      LoadLibrary();

    for(int ix=0; ix<fNx; ++ix) {
      for(int iy=0; iy<fNy; ++iy) {
	int vox_id = ix + iy * fNx;
	double vis_sum = 0.;
	for(auto const& vis_pmt : *(fTheLibrary->GetCounts(vox_id)))
	  vis_sum += ((double)(vis_pmt));
	result[ix][iy] = vis_sum;
      }
    }
    return result;
  }

  float PhotonVisibilityService::Fraction2AbsoluteX(float frac) const
  { assert(frac>=0. && frac <=1.); return fXmin + (fXmax - fXmin) * frac; }
  float PhotonVisibilityService::Fraction2AbsoluteY(float frac) const
  { assert(frac>=0. && frac <=1.); return fYmin + (fYmax - fYmin) * frac; }
  float PhotonVisibilityService::Fraction2AbsoluteZ(float frac) const
  { assert(frac>=0. && frac <=1.); return fZmin + (fZmax - fZmin) * frac; }

  //--------------------------------------------------------------------
  void PhotonVisibilityService::LoadLibrary() const
  {
    // Don't do anything if the library has already been loaded.
    std::cout<<"Loading library..."<<std::endl;
    if(fTheLibrary == 0) {
      fTheLibrary = new PhotonLibrary();


      if((!fLibraryBuildJob)&&(!fDoNotLoadLibrary)) {

	std::string LibraryFileWithPath;
	if ( fLibraryFile.front()!='/' ) {
	  // relative path provided. use designated folder in larlite
	  LibraryFileWithPath = std::string(getenv("FMATCH_DATADIR"));
	  LibraryFileWithPath += "/";
	  LibraryFileWithPath += fLibraryFile;
	}
	else {
	  // absolute path
	  LibraryFileWithPath = fLibraryFile;
	}
	if(!fParameterization) {
	  std::cout << "PhotonVisibilityService Loading photon library from file "
		    << LibraryFileWithPath
		    << std::endl;
	  size_t NVoxels = GetVoxelDef().GetNVoxels();
	  fTheLibrary->LoadLibraryFromFile(LibraryFileWithPath, NVoxels);
	}
      }
      else {
	// building library, so creating an empty one
        //art::ServiceHandle<geo::Geometry> geom;

        size_t NVoxels = GetVoxelDef().GetNVoxels();
	std::cout << " Vis service running library build job.  Please ensure "
		  << " job contains LightSource, LArG4, SimPhotonCounter"<<std::endl;
	fTheLibrary->CreateEmptyLibrary(NVoxels, fNOpDetChannels);
      }
    }
  }

  //--------------------------------------------------------------------
  void PhotonVisibilityService::StoreLibrary()
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    if(fLibraryBuildJob )
      {
	std::cout<< " Vis service "
		 << " Storing Library entries to file..." <<std::endl;
	fTheLibrary->StoreLibraryToFile(fLibraryFile);
      }
  }


  //------------------------------------------------------

  // Eventually we will calculate the light quenching factor here
  /*
  double PhotonVisibilityService::GetQuenchingFactor(double)
  {
    // for now, no quenching
    return 1.0;

  }
  */

  //------------------------------------------------------

  // Get a vector of the relative visibilities of each OpDet
  //  in the event to a point xyz

  const std::vector<float>* PhotonVisibilityService::GetAllVisibilities(double * xyz) const
  {
    int VoxID = fVoxelDef.GetVoxelID(xyz);
    return GetLibraryEntries(VoxID);
  }


  //------------------------------------------------------

  // Get distance to optical detector OpDet
  /*
  double PhotonVisibilityService::DistanceToOpDet( double* xyz, unsigned int OpDet )
  {
    art::ServiceHandle<geo::Geometry> geom;
    return geom->OpDetGeoFromOpDet(OpDet).DistanceToPoint(xyz);
  }
  */

  //------------------------------------------------------


  // Get the solid angle reduction factor for planar optical detector OpDet
  /*
  double PhotonVisibilityService::SolidAngleFactor( double* xyz, unsigned int OpDet )
  {
    art::ServiceHandle<geo::Geometry> geom;
    return geom->OpDetGeoFromOpDet(OpDet).CosThetaFromNormal(xyz);
  }
  */

  //------------------------------------------------------

  float PhotonVisibilityService::GetVisibility(double * xyz, unsigned int OpChannel) const
  {
    int VoxID = fVoxelDef.GetVoxelID(xyz);
    return GetLibraryEntry(VoxID, OpChannel);
  }

  float PhotonVisibilityService::GetVisibility(double x, double y, double z, unsigned int OpChannel) const
  {
    int VoxID = fVoxelDef.GetVoxelID(x,y,z);
    return GetLibraryEntry(VoxID, OpChannel);
  }


  //------------------------------------------------------

  void PhotonVisibilityService::StoreLightProd(int VoxID, double N)
  {
    fCurrentVoxel = VoxID;
    fCurrentValue = N;
    std::cout << " PVS notes production of " << N << " photons at Vox " << VoxID<<std::endl;
  }


  //------------------------------------------------------


  void PhotonVisibilityService::RetrieveLightProd(int& VoxID, double& N) const
  {
    N     = fCurrentValue;
    VoxID = fCurrentVoxel;
  }

  //------------------------------------------------------

  void PhotonVisibilityService::SetLibraryEntry(int VoxID, int OpChannel, float N)
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    fTheLibrary->SetCount(VoxID,OpChannel, N);
    std::cout<< " PVS logging " << VoxID << " " << OpChannel<<std::endl;
  }

  //------------------------------------------------------



  const std::vector<float>* PhotonVisibilityService::GetLibraryEntries(int VoxID) const
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    return fTheLibrary->GetCounts(VoxID);
  }

  //------------------------------------------------------

  float PhotonVisibilityService::GetLibraryEntry(int VoxID, int Channel) const
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    return fTheLibrary->GetCount(VoxID, Channel);
  }

  //------------------------------------------------------
  /*
  int PhotonVisibilityService::NOpChannels()
  {
    if(fTheLibrary == 0)
      LoadLibrary();

    return fTheLibrary->NOpChannels();
  }
  */
} // namespace
