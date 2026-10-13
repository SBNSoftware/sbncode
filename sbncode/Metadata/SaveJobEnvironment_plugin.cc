/**
 * @file   sbncode/Metadata/SaveJobEnvironment_plugin.cc
 * @brief  Producer module writing job environment information into output.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 16, 2025
 */

// SBN libraries
#include "sbncode/Metadata/JobEnvironmentInfoExtractor.h"
#include "sbnobj/Common/Metadata/JobEnvironmentInfo.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"

// ROOT libraries
#include "TFile.h"
#include "TDirectory.h"

// C/C++ libraries
#include <iterator> // std::back_inserter()
#include <memory> // std::make_unique()
#include <optional>
#include <utility> // std::move()
#include <vector>



// -----------------------------------------------------------------------------
namespace sbn { class SaveJobEnvironment; }
/**
 * @brief Writes information from the job execution environment into output.
 * 
 * This result-level module writes into the `art::Results` an object with
 * information about the environment where the job is being executed.
 * 
 * Information is extracted by the `sbn::JobEnvironmentInfoExtractor` algorithm.
 * 
 * To add a new repository to the source of metadata collected by this plugin,
 * see the instructions and explanations in
 * @ref SBNsourceMetadataSystem "the SBN metadata system documentation".
 * 
 * 
 * Output
 * -------
 * 
 * In _art_/ROOT output file, `art::Results` level:
 *  * `std::vector<sbn::JobEnvironmentInfo>` objects including all the ones
 *    found in the input files, and in addition, for this job:
 *      * a snapshot of the full environment (like in `getEnv()` C function),
 *        lexicographically sorted as for the object requirement.
 * 
 * The information from this job is also saved in the `TFileService` output
 * file, if that service is configured. Input information is currently not
 * written into the `TFileService` output. Also note that the information is
 * saved in the main directory of the file, as opposed to the usual subdirectory
 * named after the module. The name of the object is `EnvInfo`, and, if
 * `sbn::JobEnvironmentInfo` dictionary library is available to ROOT, one simple
 * way to see it is to open the ROOT file with the ROOT interpreter and execute
 * `EvtInfo->dump(std::cout);`.
 * 
 * 
 * Module configuration
 * ---------------------
 * 
 * * `Repositories` (list of strings; default: see `DefaultRepositories`):
 *     the list of repositories whose version reporting tool will be queried.
 *     Directly passed as `JobEnvironmentInfoExtractor::Config::repositories`
 *     (see) to the exctraction algorithm.
 * * `WriteToTFileService` (flag, optional): if set to `false`, the information
 *     will not be written into the `TFileService` output file; if it is set to
 *     `true`, it will be written into the output of file the service, which is
 *     required to be configured; if the parameter is omitted, the information
 *     will be written only if the service is available.
 * 
 * 
 * Service dependencies
 * ---------------------
 * 
 * * `TFileService` (required only if `WriteToTFileService` is set to `true`)
 *     for output into the `TFile` managed by that service.
 * 
 */
class sbn::SaveJobEnvironment: public art::ResultsProducer {
  
    public:
  
  /// List of repository tools to be loaded by default.
  static std::vector<std::string> const DefaultRepositories;
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<std::string> Repositories{
      Name{ "Repositories" },
      Comment{ "Names of the repositories with version reporting tools" },
      DefaultRepositories
      };
    
    fhicl::OptionalAtom<bool> WriteToTFileService{
      Name{ "WriteToTFileService" },
      Comment
        { "enables or disables writing to TFileService (default: if available)" }
      };
    
  }; // Config
  
  using Parameters = Table<Config>;
  
  
  /// Constructor.
  SaveJobEnvironment(Parameters const& params);
  
  /// Fetches the information.
  virtual void beginJob() override;
  
  /// Reads the information on each new input file.
  virtual void readResults(art::Results const& results) override;

  /// Writes the information at output closure time.
  virtual void writeResults(art::Results& results) override;
  
  /// Clears after writing: does nothing.
  virtual void clear() override;
  
  /// Write information to `TFileService` file if available.
  virtual void endJob() override;
  
  
    private:
  
  // --- BEGIN --- Configuration -----------------------------------------------
  
  /// Names of the repositories to query for version metadata.
  std::vector<std::string> const fRepositoryNames;
  
  bool const fWriteToTFileService; ///< Whether to write info to `TFileService`.
  
  // --- END ----- Configuration -----------------------------------------------
  
  
  sbn::JobEnvironmentInfoExtractor fInfoExtractor; ///< Extraction algorithm.
  
  ///< All information from the input to be written, in the order it was read.
  std::vector<sbn::JobEnvironmentInfo> fInputInfo;
  
  /// Information from this job to be written.
  sbn::JobEnvironmentInfo fJobInfo;
  
  
  /// Fetches and returns all the information.
  sbn::JobEnvironmentInfo fetchInformation();
  
  
  /// Writes the current information into `outDir`.
  void writeInformationToTFile(TDirectory& outDir) const;
  
  /// Returns a correctly configured information extractor algorithm.
  sbn::JobEnvironmentInfoExtractor makeInfoExtractor() const;
  
  
  /// Returns whether we should write to `TFileService`.
  static bool parseWriteToTFileService
    (std::optional<bool> writeToTFileService);
  
  
}; // sbn::SaveJobEnvironment


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
std::vector<std::string> const sbn::SaveJobEnvironment::DefaultRepositories {
  "sbncode", "sbndcode", "icaruscode"
};


// -----------------------------------------------------------------------------
sbn::SaveJobEnvironment::SaveJobEnvironment(Parameters const& params)
  : fRepositoryNames{ params().Repositories() }
  , fWriteToTFileService
      { parseWriteToTFileService(params().WriteToTFileService()) }
  , fInfoExtractor{ makeInfoExtractor() }
{
  
  produces<std::vector<sbn::JobEnvironmentInfo>>();
  
  if (fWriteToTFileService) {
    mf::LogInfo{ "SaveJobEnvironment" }
      << "Will also save information into TFileService output file.";
  }
  
} // sbn::SaveJobEnvironment::SaveJobEnvironment()


// -----------------------------------------------------------------------------
void sbn::SaveJobEnvironment::beginJob() {
  
  fJobInfo = fetchInformation();
  
}


// -----------------------------------------------------------------------------
void sbn::SaveJobEnvironment::readResults(art::Results const& results) {
  
  std::vector const infoHandles
    = results.getMany<std::vector<sbn::JobEnvironmentInfo>>();
  
  mf::LogDebug{ "SaveJobEnvironment" }
    << "Found " << infoHandles.size() << " job information products in input.";
  
  for (art::Handle<std::vector<sbn::JobEnvironmentInfo>> const& infoHandle
    : infoHandles)
  {
    if (!infoHandle.isValid()) continue;
    std::copy
      (infoHandle->begin(), infoHandle->end(), back_inserter(fInputInfo));
  } // for
  
} // sbn::SaveJobEnvironment::readResults()


// -----------------------------------------------------------------------------
void sbn::SaveJobEnvironment::writeResults(art::Results& results) {

  mf::LogDebug{ "SaveJobEnvironment" }
    << "Information saved into art/ROOT output file:\n\n" << fJobInfo;
  
  auto allInfo
    = std::make_unique<std::vector<sbn::JobEnvironmentInfo>>(fInputInfo);
  allInfo->push_back(fJobInfo);
  
  // one copy each call
  results.put(std::move(allInfo));
  
}


// -----------------------------------------------------------------------------
void sbn::SaveJobEnvironment::clear() {
}


// -----------------------------------------------------------------------------
void sbn::SaveJobEnvironment::endJob() {
  
  if (fWriteToTFileService) {
    
    writeInformationToTFile(art::ServiceHandle<art::TFileService>()->file());
    
  }
  
} // sbn::SaveJobEnvironment::endJob()


// -----------------------------------------------------------------------------
sbn::JobEnvironmentInfo sbn::SaveJobEnvironment::fetchInformation() {
  
  return fInfoExtractor.extract(moduleDescription());
  
}


// -----------------------------------------------------------------------------
void sbn::SaveJobEnvironment::writeInformationToTFile(TDirectory& outDir) const
{
  
  // assuming copy constructor will be used
  outDir.WriteObject(&fJobInfo, "EnvInfo");
  
} // sbn::SaveJobEnvironment::writeInformationToTFile()


// -----------------------------------------------------------------------------
bool sbn::SaveJobEnvironment::parseWriteToTFileService
  (std::optional<bool> writeToTFileService)
{
  // if is requested that we don't write it to TFileService, stop it here:
  if (!writeToTFileService.value_or(true)) return false;
  
  // so either there is no explicit request, or that request is to use it:
  // at this point, we do need to know if TFileService is available
  // and if it's there, we'll want to write into it
  try {
    art::ServiceHandle<art::TFileService>();
    return true;
  }
  catch(art::Exception& e) {
    if (e.categoryCode() != art::errors::ServiceNotFound)
      throw; // something else entirely is happening: propagate the exception
    
    // not available: do we need it?
    if (!writeToTFileService) return false; // no explicit request, so no
    
    // there is a request (and if we are here it was positive): complain!
    throw art::Exception{ art::errors::Configuration, "", e }
      << "SaveJobEnvironment explicitly requested saving into TFileService"
        " output, but TFileService is not configured.\n";
    
  }
  
} // sbn::SaveJobEnvironment::parseWriteToTFileService()


// -----------------------------------------------------------------------------
sbn::JobEnvironmentInfoExtractor
sbn::SaveJobEnvironment::makeInfoExtractor() const {
  
  JobEnvironmentInfoExtractor::Config config;
  config.repositories = fRepositoryNames;
  return JobEnvironmentInfoExtractor{ config };
  
} // sbn::SaveJobEnvironment::makeInfoExtractor()


// -----------------------------------------------------------------------------
DEFINE_ART_RESULTS_PLUGIN(sbn::SaveJobEnvironment)


// -----------------------------------------------------------------------------
