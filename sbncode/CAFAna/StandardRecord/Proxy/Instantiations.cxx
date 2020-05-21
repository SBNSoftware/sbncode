// This file is the only way BasicTypesProxy.cxx gets compiled at all (the
// srproxy package doesn't include binaries).
#include "SRProxy/BasicTypesProxy.cxx"

#include "sbnanalysis/core/Experiment.hh"

// But this also gives us an opportunity to instantiate the template for
// various sbndcode-specific enums that would otherwise be missing symbols.
namespace caf
{
  template class Proxy<Experiment>;
}
