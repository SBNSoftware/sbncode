/*
  PhotonLibrary basically copied from

 */

#ifndef PHOTONLIBRARY_H
#define PHOTONLIBRARY_H

#include "TTree.h"
#include "PhotonVoxels.h"
#include <vector>
#include <string>

namespace phot{
  
  class PhotonLibrary
  {
  public:
    PhotonLibrary();
    ~PhotonLibrary();

    float GetCount(size_t Voxel, size_t OpChannel);
    void   SetCount(size_t Voxel, size_t OpChannel, float Count);
    
    const std::vector<float>* GetCounts(size_t Voxel) const;
    inline const std::vector<std::vector<float> >& GetData() const
    { return fLookupTable; }
    
    void StoreLibraryToFile(std::string LibraryFile);
    void LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels);
    void CreateEmptyLibrary(size_t NVoxels, size_t NChannels);
    

    int NOpChannels() const { return fNOpChannels; }
    int NVoxels() const { return fNVoxels; }
    
  private:
    // fLookupTable[Voxel]->at(OpChannel) = Count
    std::vector<std::vector<float> > fLookupTable;
    size_t fNOpChannels;
    size_t fNVoxels;
    
    static std::vector<float> EmptyChannelsList;
    
    static std::vector<float>* EmptyList()
      { EmptyChannelsList.clear(); return &EmptyChannelsList; }
  };

}

#endif
