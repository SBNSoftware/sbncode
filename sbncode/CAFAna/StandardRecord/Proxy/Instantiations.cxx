// This file is the only way BasicTypesProxy.cxx gets compiled at all (the
// srproxy package doesn't include binaries).
#include "SRProxy/BasicTypesProxy.cxx"

#include "sbncode/StandardRecord/SREnums.h"

// But this also gives us an opportunity to instantiate the template for
// various sbncode-specific enums that would otherwise be missing symbols.
namespace caf
{
  template class Proxy<Det_t>;
  template class Proxy<Plane_t>;
  template class Proxy<Wall_t>;
  template class Proxy<MCType_t>;
  template class Proxy<generator_>;
  template class Proxy<interaction_mode_>;
  template class Proxy<genie_status_>;
  template class Proxy<g4_process_>;
}
