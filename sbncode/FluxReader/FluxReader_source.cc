#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "FluxReader.h"

namespace fluxr {
  typedef art::Source<FluxReader> FluxReaderSource;
}

DEFINE_ART_INPUT_SOURCE(fluxr::FluxReaderSource)
