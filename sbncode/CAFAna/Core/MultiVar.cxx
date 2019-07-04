#include "CAFAna/Core/MultiVar.h"

#include <algorithm>
#include <map>
#include <set>

namespace ana
{
  // Stupid hack to avoid colliding with the IDs of actual Vars. Just count
  // down through negative numbers.
  int MultiVar::fgNextID = -1;
}
