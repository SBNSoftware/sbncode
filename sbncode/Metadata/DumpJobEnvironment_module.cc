/**
 * @file   sbncode/Metadata/DumpJobEnvironment_module.cc
 * @brief  Producer module writing job environment information into output.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 16, 2025
 */

// local libraries
#include "sbnobj/Common/Metadata/JobEnvironmentInfo.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/OutputModule.h"
#include "art/Framework/Principal/ResultsPrincipal.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/ConfigurationTable.h" // fhicl::WrappedTable
#include "fhiclcpp/types/TableFragment.h"
#include "fhiclcpp/types/Atom.h"

// C++ standard libraries
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
namespace sbn { class DumpJobEnvironment; }
/**
 * @brief Output module dumping input versions to screen.
 * 
 * The output module can be added to any of the end paths of an _art_ job to
 * get a complete dump of the SBN job environment metadata stored into the input
 * file.
 * 
 * For example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * process_name: JobMeta
 * 
 * services.message.destinations.MetadataLog: {
 *   type:       file
 *   filename:  "JobEnvironment.log"
 *   append:     false
 *   threshold:  INFO
 *   categories: {
 *     DumpJobEnvironment: { limit: -1 }
 *     default:            { limit: 0 }
 *   }
 * }
 * 
 * outputs.metadataDumper: { module_type: "DumpJobEnvironment" }
 * 
 * physics.streams: [ metadataDumper ]
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will dump the metadata into a text file named `JobEnvironment.log`.
 * 
 */
class sbn::DumpJobEnvironment: public art::OutputModule {
    public:
  
  /// Module configuration.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::TableFragment<art::OutputModule::Config> OutputModuleConfig;
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the messagefacility output category to be used" },
      "DumpJobEnvironment"
      };
    
  }; // Config

  using Parameters
    = fhicl::WrappedTable<Config, art::OutputModule::Config::KeysToIgnore>;

  explicit DumpJobEnvironment(Parameters const& params);

    private:
  
  // --- BEGIN --- Configuration -----------------------------------------------
  
  std::string const fLogCategory; ///< Messagefacility category for the output.
  
  // --- END ----- Configuration -----------------------------------------------
  
  
  /// Dumps the information from an handle into the output stream.
  template <typename Stream>
  void dumpInformation(
    Stream& out,
    art::Handle<std::vector<sbn::JobEnvironmentInfo>> const& infoHandle
    ) const;
  
  
  
  void write(art::EventPrincipal&) override {}
  void writeRun(art::RunPrincipal&) override {}
  void writeSubRun(art::SubRunPrincipal&) override {}
  
  /// Reads and prints all the metadata data products.
  void readResults(art::ResultsPrincipal const& results) override;

}; // sbn::DumpJobEnvironment


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
sbn::DumpJobEnvironment::DumpJobEnvironment(Parameters const& params)
  : OutputModule{ params().OutputModuleConfig }
  , fLogCategory{ params().LogCategory() }
{}


// -----------------------------------------------------------------------------
void sbn::DumpJobEnvironment::readResults
  (art::ResultsPrincipal const& principal)
{
  if (!principal.size()) return;
  
  art::ModuleContext const moduleContext{ moduleDescription() };
  art::Results const& results = principal.makeResults(moduleContext);
  
  std::vector<art::Handle<std::vector<sbn::JobEnvironmentInfo>>> infoHandles
    = results.getMany<std::vector<sbn::JobEnvironmentInfo>>();
  
  mf::LogInfo out{ fLogCategory };
  out << "Found " << infoHandles.size() << " job information entries in input.";
  
  for (art::Handle<std::vector<sbn::JobEnvironmentInfo>> const& infoHandle
    : infoHandles
  ) {
    out << '\n' << std::string(80, '*') << '\n';
    dumpInformation(out, infoHandle);
  }
  
} // sbn::DumpJobEnvironment::readResults()


// -----------------------------------------------------------------------------
template <typename Stream>
void sbn::DumpJobEnvironment::dumpInformation(
  Stream& out,
  art::Handle<std::vector<sbn::JobEnvironmentInfo>> const& infoHandle
) const {
  
  art::Provenance const* provenance = infoHandle.provenance();
  if (provenance) {
    out << "Information from '" << provenance->inputTag().encode() << "'";
  }
  else {
    out << "Information with unknown provenance";
  }
  
  if (infoHandle.isValid()) {
    out << " from " << infoHandle->size() << " sources\n";
    for (sbn::JobEnvironmentInfo const& info: *infoHandle) {
      out << std::string(80, '=') << '\n' << info;
    }
  }
  else out << "\n[information not available]\n";
  
} // sbn::DumpJobEnvironment::dumpInformation()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::DumpJobEnvironment)


// -----------------------------------------------------------------------------
