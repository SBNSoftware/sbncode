#pragma once

#include <string>
#include <unordered_set>

namespace flat
{
  class IBranchPolicy
  {
  public:
    virtual bool Include(const std::string&) const = 0;
  };

  /// Branch policy based on a list loaded from a text file
  class FileListBranchPolicy: public IBranchPolicy
  {
  public:
    FileListBranchPolicy(const std::string& fname);

    bool Include(const std::string& s) const override
    {
      return fIncluded.count(s);
    }
  protected:
    std::unordered_set<std::string> fIncluded;
  };
}
