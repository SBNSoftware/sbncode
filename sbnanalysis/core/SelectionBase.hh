#ifndef __sbnanalysis_core_SelectionBase__
#define __sbnanalysis_core_SelectionBase__

/**
 * A generic processor that writes an sbnanalysis tree.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include "ProcessorBase.hh"

namespace core {

/**
 * \class core::SelectionBase
 * \brief Base class for event selections
 *
 * See core::ProcessorBase for more details.
 */
class SelectionBase : public ProcessorBase {
public:
  /** Constructor. */
  SelectionBase();

  /** Destructor. */
  virtual ~SelectionBase();
};

}  // namespace core

#endif  // __sbnanalysis_core_SelectionBase__

