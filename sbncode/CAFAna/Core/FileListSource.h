#pragma once

#include "CAFAna/Core/IFileSource.h"

#include <string>
#include <vector>

namespace ana
{
  /// Simple file source based on an explicit list provided by the user
  class FileListSource: public IFileSource
  {
  public:
    /// default offset and stride mean obey cmd-line options
    FileListSource(const std::vector<std::string>& files,
		   int stride = -1, int offset = -1);
    virtual ~FileListSource();

    virtual TFile* GetNextFile() override;
    int NFiles() const override {return fN;}
  protected:
    std::vector<std::string> fFileNames; ///< The list of files
    std::vector<std::string>::iterator fIt; ///< Iterator into \ref fFileNames
    int fStride;
    int fN; ///< Number of files that will actually be returned
    TFile* fFile; ///< The most-recently-returned file
    static bool fgGotTickets; ///< Have we renewed our tickets?
  };
}
