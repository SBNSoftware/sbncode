////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataSBN_service.cc.  
//
// Purpose:  Implementation for FileCatalogMetadataSBN.
//
// Created:  28-Oct-2014,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include "sbncode/Metadata/FileCatalogMetadataSBN.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "cetlib_except/exception.h"

//--------------------------------------------------------------------
// Constructor.

util::FileCatalogMetadataSBN::
FileCatalogMetadataSBN(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
{
  // Insist on configuring Experiment from the fcl file (ideally) or the
  // environment.
  const char* expt = getenv("EXPERIMENT");
  if(expt) fExperiment = pset.get<std::string>("Experiment", expt); else fExperiment = pset.get<std::string>("Experiment");
  std::transform(fExperiment.begin(), fExperiment.end(), fExperiment.begin(), [](unsigned char c){return std::tolower(c);});

  // Get parameters.

  fFCLName = pset.get<std::string>("FCLName");
  fProjectName = pset.get<std::string>("ProjectName");
  fProjectStage = pset.get<std::string>("ProjectStage");
  fProjectVersion = pset.get<std::string>("ProjectVersion");    
  fProjectSoftware = pset.get<std::string>("ProjectSoftware","");    
  fProductionName = pset.get<std::string>("ProductionName","");  //Leave as default value if not running a production   
  fProductionType = pset.get<std::string>("ProductionType",""); //Leave as default value if not running a production
  fMerge = pset.get<int>("Merge", -1);
  fParameters = pset.get<std::vector<std::pair<std::string, std::string>>>("Parameters", std::vector<std::pair<std::string, std::string>>());

  // Register for callbacks.

  reg.sPostBeginJob.watch(this, &FileCatalogMetadataSBN::postBeginJob);
}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via FileCatalogMetadata service.
void util::FileCatalogMetadataSBN::postBeginJob()
{
  // Get art metadata service.

  art::ServiceHandle<art::FileCatalogMetadata> mds;

  // Add metadata.

  mds->addMetadata("fcl.name", fFCLName);
  mds->addMetadata(fExperiment + "_project.name", fProjectName);
  mds->addMetadata(fExperiment + "_project.stage", fProjectStage);
  mds->addMetadata(fExperiment + "_project.version", fProjectVersion);
  mds->addMetadata(fExperiment + "_project.software", fProjectSoftware);
  mds->addMetadata("production.name", fProductionName);
  mds->addMetadata("production.type", fProductionType);
  std::ostringstream ostr;
  if(fMerge >= 0) {
    if(fMerge > 0)
      mds->addMetadata("merge.merge", "1");
    else
      mds->addMetadata("merge.merge", "0");
    mds->addMetadata("merge.merged", "0");
  }
  for(auto const& param : fParameters)
    mds->addMetadata(param.first, param.second);
}

DEFINE_ART_SERVICE(util::FileCatalogMetadataSBN)
