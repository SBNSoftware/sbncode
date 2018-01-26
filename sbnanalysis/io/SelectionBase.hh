#ifndef __sbnanalysis_io_SelectionBase__
#define __sbnanalysis_io_SelectionBase__

/**
 * A generic processor that writes an sbnanalysis tree.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include "ProcessorBase.hh"

namespace io {

class SelectionBase : public ProcessorBase {
public:
  SelectionBase();
  virtual ~SelectionBase();
};

}  // namespace io

#endif  // __sbnanalysis_io_SelectionBase__

