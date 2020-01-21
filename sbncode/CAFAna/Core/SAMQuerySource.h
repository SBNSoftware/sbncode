#include "CAFAna/Core/WildcardSource.h"

namespace ana
{
  /// \brief File source based on a SAM query or dataset (definition)
  ///
  /// Locates the files on bluearc or pnfs (bluearc preferred).
  class SAMQuerySource: public FileListSource
  {
  public:
    /// \param query May be a SAM dataset name or a SAM query string
    SAMQuerySource(const std::string& query, int stride = -1, int offset = -1);
    virtual ~SAMQuerySource();
  protected:
    std::vector<std::string> LocationsForSAMQuery(const std::string& str,
                                                  int stride, int offset);
    /// Take filenames, return locations suitable for TFile::Open()
    std::vector<std::string> LocateSAMFiles(const std::vector<std::string>& fnames);

    bool RunningOnGrid() const;
    std::string EnsureDataset(const std::string& query) const;
    std::string EnsureSnapshot(const std::string& def) const;
  };
}
