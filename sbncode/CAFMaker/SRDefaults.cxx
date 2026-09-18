/**
 * @file   sbncode/CAFMaker/SRDefaults.cxx
 * @brief  Discovery of default values in Standard Record classes.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 6, 2026
 * @see    sbncode/CAFMaker/SRDefaults.h
 */

#include "sbncode/CAFMaker/SRDefaults.h"


// -----------------------------------------------------------------------------
namespace {
  
  template<typename... SRobjs>
  std::tuple<SRobjs...> staticInit
    (std::tuple<SRobjs...> objs)
    { (std::get<SRobjs>(objs).setDefault(), ...); return objs; }

} // local namespace


// -----------------------------------------------------------------------------
caf::SRDefaults::SupportedTypes const caf::SRDefaults::defaultObjs
  = staticInit(caf::SRDefaults::SupportedTypes{});

// -----------------------------------------------------------------------------
