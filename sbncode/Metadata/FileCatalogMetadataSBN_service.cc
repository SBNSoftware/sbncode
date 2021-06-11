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
  // Get parameters.

  fFCLName = pset.get<std::string>("FCLName");
  fProjectName = pset.get<std::string>("ProjectName");
  fProjectStage = pset.get<std::string>("ProjectStage");
  fProjectVersion = pset.get<std::string>("ProjectVersion");    
  fProjectSoftware = pset.get<std::string>("ProjectSoftware","");    
  fProductionName = pset.get<std::string>("ProductionName","");  //Leave as default value if not running a production   
  fProductionType = pset.get<std::string>("ProductionType",""); //Leave as default value if not running a production
  fMerge = pset.get<int>("Merge", -1);
  fParameters = pset.get<std::vector<std::string> >("Parameters", std::vector<std::string>());

  // It doesn't make sense for parameter vector to have an odd number of elements.

  if(fParameters.size() % 2 != 0) {
    throw cet::exception("FileCatalogMetadataSBN") 
      << "Parameter vector has odd number of elements.\n";
  }

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
  std::ostringstream ostr;
  if(fMerge >= 0) {
    if(fMerge > 0)
      mds->addMetadata("merge.merge", "1");
    else
      mds->addMetadata("merge.merge", "0");
    mds->addMetadata("merge.merged", "0");
  }
  for(unsigned int i=0; i<fParameters.size(); i += 2)
    mds->addMetadata(fParameters[i], fParameters[i+1]);
}

DEFINE_ART_SERVICE(util::FileCatalogMetadataSBN)
