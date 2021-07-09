////////////////////////////////////////////////////////////////////////
//
// Name:  FileCatalogMetadataSBN.h.  
//
// Purpose:  Art service adds microboone-specific per-job sam metadata.
//
//           FCL parameters:
//
//           FCLName         - FCL file name.
//           ProjectName     - Project name.
//           ProjectStage    - Project stage.
//           ProjectVersion  - Project version.
//           ProjectSoftware - Project software information.
//           ProductionName  - Production campaign name.
//           ProductionType  - Production campaign type.
//           Merge           - Merge flag.
//                              1 - Set merge.merge = 1 and merge.merged = 0
//                              0 - Set merge.merge = 0 and merge.merged = 0
//                             -1 - Do not generate merge parameters.
//           Parameters      - Arbitrary (key, value) parameters.
//                             Specify in fcl file as sequence-of-sequence:
//                             [[key1, value1], [key2, value2],...]
//           POTModuleLabel  - POTSummary module label (default "generator").
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
    int GetMerge() const {return fMerge;}
    const std::vector<std::pair<std::string, std::string>>& GetParameters() const {return fParameters;}

  private:

    // Callbacks.

    void postBeginJob();
    void postEndSubRun(art::SubRun const& subrun);
    void preCloseOutputFile(const std::string& label);

    // Data members.

    std::string fExperiment;
    std::string fFCLName;
    std::string fProjectName;
    std::string fProjectStage;
    std::string fProjectVersion;
    std::string fProjectSoftware;
    std::string fProductionName; //Production parameter, do not use if not running a production
    std::string fProductionType; //Production parameter, do not use if not running a production
    int fMerge;
    std::vector<std::pair<std::string, std::string>> fParameters;
    std::string fPOTModuleLabel;
    double fTotPOT;
  };

} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataSBN, LEGACY)

#endif
