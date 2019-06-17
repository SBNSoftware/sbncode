#include "CAFAna/Core/SAMProjectSource.h"

#include "ifdh.h"

#include <unistd.h>

#ifdef DARWINBUILD
#include <libgen.h>
#endif

#include "TFile.h"

namespace ana
{
  bool SAMProjectSource::fgGotTickets = false;

  //----------------------------------------------------------------------
  SAMProjectSource::SAMProjectSource(const std::string& proj, int fileLimit)
    : fIFDH(new ifdh_ns::ifdh)
    , fFile(0)
    , fNFiles(fileLimit)
  {
    if(!fgGotTickets){
      // No kerberos ticket means no point trying to voms-proxy-init. It likely
      // also means we're in a grid job, where that would be counterproductive
      // anyway.
      if(system("klist -5 -s || klist -s") != 0) fgGotTickets = true;
    }

    if(!fgGotTickets){
      // This comes from NovaGridUtils, v02.10 onwards.
      system("setup_fnal_security -b");

      fgGotTickets = true;
    }

    // If X509_USER_PROXY isn't set, set it manually. Can help unconfuse IFDH
    // about what x509 it's supposed to be using.
    if(getenv("X509_USER_PROXY") == 0)
      setenv("X509_USER_PROXY", TString::Format("/tmp/x509up_u%d", getuid()).Data(), 0);

    fIFDH->set_debug("0"); // shut up

    fProjectURL = fIFDH->findProject(proj, "nova");

    fProcessID = fIFDH->establishProcess(fProjectURL, "CAFAna", "v0.9", getenv("HOSTNAME"), getenv("USER"), "nova", "", fileLimit);
  }

  //----------------------------------------------------------------------
  SAMProjectSource::~SAMProjectSource()
  {
    if(fFile){
      // Tidy up the final file
      const std::string fname = fFile->GetName();
      delete fFile;
      unlink(fname.c_str());
    }

    // End the process cleanly
    fIFDH->endProcess(fProjectURL, fProcessID);

    // certainly wrong for fileLimit case
    // status = fIFDH.endProject(fProjectURL);

    fIFDH->cleanup();
  }

  //----------------------------------------------------------------------
  TFile* SAMProjectSource::GetNextFile()
  {
    if(fFile){
      // Tidy up the previous file
      const std::string fname = fFile->GetName();
      delete fFile;
      fFile = 0;
      unlink(fname.c_str());

      // And let SAM know we're done with it
      fIFDH->updateFileStatus(fProjectURL, fProcessID, fname, "consumed");
    }

    std::string tmp;

    tmp = fIFDH->getNextFile(fProjectURL, fProcessID);
    const std::string uri(tmp);
    if(uri.empty()) return 0; // out of files

    tmp = fIFDH->fetchInput(uri);
    const std::string fname(tmp);
    assert(!fname.empty());

    // Let SAM know we got it OK
    fIFDH->updateFileStatus(fProjectURL, fProcessID, fname, "transferred");

    // Additional newlines because ifdh currently spams us with certificate
    // messages.
    if(fNFiles < 0) std::cout << std::endl << "Processing " << basename((char *)fname.c_str()) << std::endl << std::endl;

    fFile = new TFile(fname.c_str());
    return fFile;
  }
} // namespace

