#include "CAFAna/Core/MultiVar.h"

#include <algorithm>
#include <map>
#include <set>

namespace ana
{
  // explicitly instantiate the templates for the types we know we have
  template class GenericMultiVar<caf::SRSpillProxy>;
  template class GenericMultiVar<caf::SRSliceProxy>;

  // Stupid hack to avoid colliding with the IDs of actual Vars. Just count
  // down through negative numbers.
  template<class T> int GenericMultiVar<T>::fgNextID = -1;
}
