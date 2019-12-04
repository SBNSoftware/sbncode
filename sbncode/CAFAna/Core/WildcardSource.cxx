#include "CAFAna/Core/WildcardSource.h"

#include "CAFAna/Core/Utilities.h"

#include "ifdh.h"

#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
namespace
{
  // No-one else needs to see this
  std::vector<std::string> sorted(const std::vector<std::string>& v)
  {
    std::vector<std::string> ret;
    ret = v;
    std::sort(ret.begin(), ret.end());
    return ret;
  }
}

namespace ana
{
  //----------------------------------------------------------------------
  WildcardSource::WildcardSource(const std::string& wildcard,
                                 int stride, int offset)
    // MRCCLoader needs all its lists in matching order, so sort here
    : FileListSource(sorted(WildcardOrLs(wildcard)), stride, offset)
  {
  }

  //----------------------------------------------------------------------
  WildcardSource::~WildcardSource()
  {
  }

  //----------------------------------------------------------------------
  std::vector<std::string> WildcardSource::
  WildcardOrLs(const std::string& wildcard) const
  {
    std::vector<std::string> ret = Wildcard(wildcard);

    struct stat ss;
    // If we found nothing, but it's because pnfs isn't mounted, try ifdh ls
    // instead.
    if(ret.empty() &&
       wildcard.find("/pnfs/") == 0 &&
       stat("/pnfs/", &ss) != 0){

      std::cout << "No files matching " << wildcard << " but that's probably because /pnfs is not mounted on the current node." << std::endl;

      // This doesn't work because "ifdh ls" doesn't actually understand
      // wildcard syntax properly.
      //      ifdh i;
      //      ret = i.ls(wildcard, 0);
    }

    return ret;
  }
}
