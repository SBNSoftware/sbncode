#pragma once

#include "CAFAna/Core/SpectrumLoaderBase.h"

namespace caf{
  typedef int Det_t;
  const int kNEARDET = 1;
  const int kFARDET = 2;
}
//#include "StandardRecord/SRHeader.h" // For Det_t

#include <map>

namespace ana
{
  class SpectrumLoader;

  /// \brief Collection of SpectrumLoaders for many configurations
  class Loaders
  {
  public:
    enum DataMC{kData, kMC};
    enum SwappingConfig{kNonSwap, kNueSwap, kNuTauSwap};
    enum FluxType{kFHC, kRHC};

    /// No loaders initialized. Use \ref SetLoaderPath to configure
    Loaders();
    ~Loaders();

    /// Configure loader via wildcard \a path
    void SetLoaderPath(const std::string& path,
                       caf::Det_t det,
                       DataMC datamc,
                       DataSource src = kBeam,
                       SwappingConfig swap = kNonSwap);

    /// Configure loader via explicit file list
    void SetLoaderFiles(const std::vector<std::string>& files,
                        caf::Det_t det,
                        DataMC datamc,
                        DataSource src = kBeam,
                        SwappingConfig swap = kNonSwap);

    void AddLoader(SpectrumLoaderBase*,
                        caf::Det_t det,
                        DataMC datamc,
                        DataSource src = kBeam,
                        SwappingConfig swap = kNonSwap);

    void DisableLoader(caf::Det_t det,
                       DataMC datamc,
                       DataSource src = kBeam,
                       SwappingConfig swap = kNonSwap);

    /// Retrieve a specific loader
    SpectrumLoaderBase& GetLoader(caf::Det_t det,
                                  DataMC datamc,
                                  DataSource src = kBeam,
                                  SwappingConfig swap = kNonSwap);

    /// Call Go() on all the loaders
    void Go();

  protected:
    typedef std::tuple<caf::Det_t, DataMC, DataSource, SwappingConfig> Key_t;

    // Hold a list of paths that have been set
    std::map<Key_t, std::string> fLoaderPaths;
    std::map<Key_t, std::vector<std::string>> fLoaderFiles;
    // Only reify them when someone actually calls GetLoader()
    std::map<Key_t, SpectrumLoaderBase*> fLoaders;

    /// We give this back when a loader isn't set for some configuration
    NullLoader fNull;
  };
} // namespace
