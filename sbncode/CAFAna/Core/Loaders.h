#pragma once

#include "CAFAna/Core/SpectrumLoaderBase.h"

#include <map>

namespace ana
{
  class SpectrumLoader;

  /// \brief Collection of SpectrumLoaders for many configurations
  class Loaders
  {
  public:
    enum DataMC{kData, kMC};
    enum SwappingConfig{kNonSwap, kNueSwap, kNuTauSwap, kIntrinsic};

    /// No loaders initialized. Use \ref SetLoaderPath to configure
    Loaders();
    ~Loaders();

    /// Configure loader via wildcard \a path
    void SetLoaderPath(const std::string& path,
                       DataMC datamc,
                       SwappingConfig swap = kNonSwap);

    /// Configure loader via explicit file list
    void SetLoaderFiles(const std::vector<std::string>& files,
                        DataMC datamc,
                        SwappingConfig swap = kNonSwap);

    void AddLoader(SpectrumLoaderBase*,
                        DataMC datamc,
                        SwappingConfig swap = kNonSwap);

    void DisableLoader(DataMC datamc,
                       SwappingConfig swap = kNonSwap);

    /// Retrieve a specific loader
    SpectrumLoaderBase& GetLoader(DataMC datamc,
                                  SwappingConfig swap = kNonSwap);

    /// Call Go() on all the loaders
    void Go();

  protected:
    typedef std::tuple<DataMC, SwappingConfig> Key_t;

    // Hold a list of paths that have been set
    std::map<Key_t, std::string> fLoaderPaths;
    std::map<Key_t, std::vector<std::string>> fLoaderFiles;
    // Only reify them when someone actually calls GetLoader()
    std::map<Key_t, SpectrumLoaderBase*> fLoaders;

    /// We give this back when a loader isn't set for some configuration
    NullLoader fNull;
  };
} // namespace
