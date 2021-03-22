////////////////////////////////////////////////////////////////////////
//
// Name:  FileCatalogMetadataSBN.h.  
//
// Purpose:  Art service adds microboone-specific per-job sam metadata.
//
//           FCL parameters:
//
//           FCLName        - FCL file name.
//           FCLVersion     - FCL file version.
//           ProjectName    - Project name.
//           ProjectStage   - Project stage.
//           ProjectVersion - Project version.
//
//           Above values will be added in internal metadata of artroot
//           output files whenever this service is included in job
//           configuration (service does not need to be called).  The
//           public interface of this service consists of accessors that
//           allow other code to discover above metadata parameters.
//
// Created:  28-Oct-2014,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#ifndef FILECATALOGMETADATASBN_H
#define FILECATALOGMETADATASBN_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace util {

  class FileCatalogMetadataSBN
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataSBN(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataSBN() = default;

    // Accessors.

    const std::string& GetFCLName() const {return fFCLName;}
    const std::string& GetProjectName() const {return fProjectName;}
    const std::string& GetProjectStage() const {return fProjectStage;}
    const std::string& GetProjectVersion() const {return fProjectVersion;}
    const std::string& GetProjectSoftware() const {return fProjectSoftware;}
    const std::string& GetProductionName() const {return fProductionName;}
    const std::string& GetProductionType() const {return fProductionType;}


  private:

    // Callbacks.

    void postBeginJob();

    // Data members.

    std::string fFCLName;
    std::string fProjectName;
    std::string fProjectStage;
    std::string fProjectVersion;
    std::string fProjectSoftware;
    std::string fProductionName; //Production parameter, do not use if not running a production
    std::string fProductionType; //Production parameter, do not use if not running a production
  };

} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataSBN, LEGACY)

#endif
