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

//--------------------------------------------------------------------
// Constructor.

util::FileCatalogMetadataSBN::
FileCatalogMetadataSBN(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
{
  // Get parameters.

  fFCLName = pset.get<std::string>("FCLName");
  fProjectName = pset.get<std::string>("ProjectName");
  fProjectStage = pset.get<std::string>("ProjectStage");
  fProjectVersion = pset.get<std::string>("ProjectVersion");    
  fProjectSoftware = pset.get<std::string>("ProjectSoftware","");    
  fProductionName = pset.get<std::string>("ProductionName","");  //Leave as default value if not running a production   
  fProductionType = pset.get<std::string>("ProductionType",""); //Leave as default value if not running a production


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
  mds->addMetadata("sbnd_project.name", fProjectName);
  mds->addMetadata("sbnd_project.stage", fProjectStage);
  mds->addMetadata("sbnd_project.version", fProjectVersion);
  mds->addMetadata("sbnd_project.software", fProjectSoftware);
  mds->addMetadata("production.name", fProductionName);
  mds->addMetadata("production.type", fProductionType);
}

DEFINE_ART_SERVICE(util::FileCatalogMetadataSBN)
