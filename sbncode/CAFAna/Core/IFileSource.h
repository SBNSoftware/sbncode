#pragma once

class TFile;

namespace ana
{
  /// \brief Interface class for accessing ROOT files in sequence
  ///
  /// Used internally by \ref SpectrumLoaderBase etc.
  class IFileSource
  {
  public:
    virtual ~IFileSource() {}
    /// \brief Returns the next file in sequence, ready for reading
    ///
    /// A null return means that the end of the sequence has been reached.
    /// DO NOT close or delete the file that is returned.
    virtual TFile* GetNextFile() = 0;

    /// May return -1 indicating the number of files is not known
    virtual int NFiles() const {return -1;}
  };
}
