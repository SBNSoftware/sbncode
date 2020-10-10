/**
 * \file PSetUtils.h
 *
 * \ingroup PSet
 * 
 * \brief Utility functions in Base/PSet
 *
 * @author Kazu - Nevis 2015
 */

/** \addtogroup PSet

    @{*/

#ifndef __PSET_UTILS_H__
#define __PSET_UTILS_H__

#include "PSet.h"

namespace flashmatch {

  /// Given a configuration string, format to create flashmatch::PSet
  //std::string FormatPSetString(std::string fname);
  /// Given a configuration file (full path), read & parse contents to create flashmatch::PSet
  std::string ConfigFile2String(std::string fname);
  /// Given a configuration file (full path), create and return flashmatch::PSet
  PSet CreatePSetFromFile(std::string fname,std::string cfg_name="cfg");

}

#endif
/** @} */ // end of doxygen group
