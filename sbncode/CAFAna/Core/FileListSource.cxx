#include "CAFAna/Core/FileListSource.h"

#include "CAFAna/Core/Utilities.h"

#include "TFile.h"

#include <cassert>
#include <iostream>

namespace ana
{
  bool FileListSource::fgGotTickets = false;

  //----------------------------------------------------------------------
  FileListSource::FileListSource(const std::vector<std::string>& files,
				 int stride, int offset)
    : fFileNames(files.begin(), files.end()),
      fIt(fFileNames.begin()),
      fStride(stride),
      fFile(0)
  {
    if(fFileNames.empty()){
      fN = 0;
      return;
    }

    if(offset < 0){
      if(getenv("CAFANA_OFFSET"))
	offset = atoi(getenv("CAFANA_OFFSET"));

      offset = std::max(offset, 0);
    }

    if(fStride < 0){
      if(getenv("CAFANA_STRIDE"))
	fStride = atoi(getenv("CAFANA_STRIDE"));

      fStride = std::max(fStride, 1);
    }

    if(fStride > int(files.size())){
      std::cerr << "Warning: stride " << fStride
                << " is greater than the number of files: " << files.size()
                << ". This is strange and inefficient." << std::endl;
      // Needs to be nonzero otherwise some callers go off into SAM-query
      // land. Having a slightly misleading progress bar is a small price to
      // pay in this weird case.
      fN = 1;
      return;
    }

    // How many files will we process from the list with this offset and
    // stride?
    fN = (int(files.size())-offset-1)/fStride+1;

    for(const std::string& loc: files){
      if(loc.rfind("/pnfs/", 0) == 0){ // ie begins with
        if(!fgGotTickets){
          // No kerberos ticket means no point trying to voms-proxy-init. It
          // likely also means we're in a grid job, where that would be
          // counterproductive anyway.
          if(system("klist -5 -s || klist -s") != 0) fgGotTickets = true;
        }

	if(!fgGotTickets){
          // This comes from NovaGridUtils, v02.10 onwards.
          //          system("setup_fnal_security -b");

	  fgGotTickets = true;
	  break;
	}
      }
    }

    for(int i = 0; i < offset; ++i){
      ++fIt;
      if(fIt == fFileNames.end()) break;
    }
  }

  //----------------------------------------------------------------------
  FileListSource::~FileListSource()
  {
    delete fFile;
  }

  //----------------------------------------------------------------------
  TFile* FileListSource::GetNextFile()
  {
    // Tidy up the last file we gave, which the caller no longer needs
    delete fFile;
    fFile = 0;

    if(fIt == fFileNames.end()) return 0; // Ran out of files

    // If the file is on pnfs rewrite it to an xrootd address
    std::string loc = *fIt;
    // loc = pnfs2xrootd(loc); // no-op for non /pnfs locations

    fFile = TFile::Open(loc.c_str()); // This pattern allows xrootd
    assert(fFile);

    for(int i = 0; i < fStride; ++i){
      if(fIt == fFileNames.end()) break;
      ++fIt; // Move on to the next file, for the subsequent call
    }

    return fFile;
  }
}
