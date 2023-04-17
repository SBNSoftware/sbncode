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
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "larcoreobj/SummaryData/POTSummary.h"

//--------------------------------------------------------------------
// Constructor.

util::FileCatalogMetadataSBN::
FileCatalogMetadataSBN(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) :
  fTotPOT(0.)
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
  fPOTModuleLabel = pset.get<std::string>("POTModuleLabel", "generator");

  // Register for callbacks.

  reg.sPostBeginJob.watch(this, &FileCatalogMetadataSBN::postBeginJob);
  reg.sPostEndSubRun.watch(this, &FileCatalogMetadataSBN::postEndSubRun);
  reg.sPreCloseOutputFile.watch(this, &FileCatalogMetadataSBN::preCloseOutputFile);
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

//--------------------------------------------------------------------
// PostEndSubrun callback.
void util::FileCatalogMetadataSBN::postEndSubRun(art::SubRun const& sr)
{

  art::ServiceHandle<art::FileCatalogMetadata> mds;

  art::Handle< sumdata::POTSummary > potListHandle;
  if(sr.getByLabel(fPOTModuleLabel,potListHandle)){
    std::lock_guard lock(fMutex);
    fTotPOT+=potListHandle->totpot;
  }
}

// PreCloseOutputFile callback.
void util::FileCatalogMetadataSBN::preCloseOutputFile(const std::string& label)
{
  std::lock_guard lock(fMutex);

  art::ServiceHandle<art::FileCatalogMetadata> mds;

  if(fTotPOT > 0.) {
    std::ostringstream streamObj;
    streamObj << fTotPOT;
    std::string strPOT = streamObj.str();
    mds->addMetadata("mc.pot", strPOT);
  }
  fTotPOT = 0.;
}

DEFINE_ART_SERVICE(util::FileCatalogMetadataSBN)
