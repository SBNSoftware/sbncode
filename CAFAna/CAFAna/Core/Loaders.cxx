#include "CAFAna/Core/Loaders.h"

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Utilities.h"

#include <cassert>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  Loaders::Loaders()
  {
  }

  //----------------------------------------------------------------------
  Loaders::~Loaders()
  {
    // for(auto it: fLoaders) delete it.second;
  }

  //----------------------------------------------------------------------
  void Loaders::SetLoaderPath(const std::string& path,
                              caf::Det_t det,
                              DataMC datamc,
                              DataSource src,
                              SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);
    assert(det == caf::kFARDET || swap == kNonSwap);
    assert(src == kBeam || swap == kNonSwap);

    const Key_t key(det, datamc, src, swap);

    // Clear out the old one if necessary
    DisableLoader(det, datamc, src, swap);

    fLoaderPaths[key] = path;
  }

  //----------------------------------------------------------------------
  void Loaders::SetLoaderFiles(const std::vector<std::string>& files,
                               caf::Det_t det,
                               DataMC datamc,
                               DataSource src,
                               SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);
    assert(det == caf::kFARDET || swap == kNonSwap);
    assert(src == kBeam || swap == kNonSwap);

    const Key_t key(det, datamc, src, swap);

    // Clear out the old one if necessary
    DisableLoader(det, datamc, src, swap);

    fLoaderFiles[key] = files;
  }

  //----------------------------------------------------------------------
  void Loaders::AddLoader(SpectrumLoaderBase* file,
                               caf::Det_t det,
                               DataMC datamc,
                               DataSource src,
                               SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);
    assert(det == caf::kFARDET || swap == kNonSwap);
    assert(src == kBeam || swap == kNonSwap);

    const Key_t key(det, datamc, src, swap);

    // Clear out the old one if necessary
    DisableLoader(det, datamc, src, swap);

    fLoaders[key] = file;
  }

  //----------------------------------------------------------------------
  void Loaders::DisableLoader(caf::Det_t det,
                              DataMC datamc,
                              DataSource src,
                              SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);
    assert(det == caf::kFARDET || swap == kNonSwap);
    assert(src == kBeam || swap == kNonSwap);

    const Key_t key(det, datamc, src, swap);

    // Clear out the current one if possible
    auto it = fLoaders.find(key);
    if(it != fLoaders.end()){
      delete it->second;
      fLoaders.erase(it);
    }

    fLoaderPaths.erase(key);
    fLoaderFiles.erase(key);
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase& Loaders::GetLoader(caf::Det_t det,
                                         DataMC datamc,
                                         DataSource src,
                                         SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);
    assert(det == caf::kFARDET || swap == kNonSwap);
    assert(src == kBeam || swap == kNonSwap);

    const Key_t key(det, datamc, src, swap);

    // Look up and return. Use fNull if no loader is set for this config
    auto itLoader = fLoaders.find(key);
    if(itLoader != fLoaders.end()) return *itLoader->second;

    auto itPath = fLoaderPaths.find(key);
    if(itPath != fLoaderPaths.end()){
      fLoaders[key] = new SpectrumLoader(itPath->second, src);
      return *fLoaders[key];
    }
    auto itFiles = fLoaderFiles.find(key);
    if(itFiles != fLoaderFiles.end()){
      fLoaders[key] = new SpectrumLoader(itFiles->second, src);
      return *fLoaders[key];
    }

    return fNull;
  }

  //----------------------------------------------------------------------
  void Loaders::Go()
  {
    for(auto it: fLoaders) it.second->Go();
  }
}
