#pragma once

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Core/Spectrum.h"

class TDirectory;

namespace ana
{

  /// Standard interface to all decomposition techniques
  class IDecomp
  {
  public:
    virtual ~IDecomp() = default;
    virtual Spectrum NCComponent()       const = 0;
    virtual Spectrum NumuComponent()     const = 0;
    virtual Spectrum AntiNumuComponent() const = 0;
    virtual Spectrum NueComponent()      const = 0;
    virtual Spectrum AntiNueComponent()  const = 0;

    virtual void SaveTo(TDirectory* dir) const = 0;
  };
}
