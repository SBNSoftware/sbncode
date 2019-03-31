#pragma once

#include "CAFAna/Analysis/Surface.h"

#include <vector>

namespace ana
{
  class MedianSurface: public Surface
  {
  public:
    MedianSurface(const std::vector<Surface>& throws);

    void DrawEnsemble(TH2* fc, Color_t color = kGray);

    void SaveTo(TDirectory * dir) const;
    static std::unique_ptr<MedianSurface> LoadFrom(TDirectory * dir);
  protected:
    std::vector<Surface> fThrows;
  };
}
