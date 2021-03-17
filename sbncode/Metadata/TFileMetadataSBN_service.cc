////////////////////////////////////////////////////////////////////////
// Name:  TFileMetadataSBN_service.cc.  
//
// Purpose:  generate SBN-specific sam metadata for root Tfiles (histogram or ntuple files).
//
// FCL parameters: GenerateTFileMetadata: This needs to be set to "true" in the fcl file
//				    to generate metadata (default value: false)
//	     	   dataTier: Currrently this needs to be parsed by the user
//		     	     for ntuples, dataTier = root-tuple; 
//		             for histos, dataTier = root-histogram
//		             (default value: root-tuple)
//	           fileFormat: This is currently specified by the user,
//			       the fileFormat for Tfiles is "root" (default value: root)	
//
// Other notes: 1. This service uses the ART's standard file_catalog_metadata service
//		to extract some of the common (common to both ART and TFile outputs)
//	        job-specific metadata parameters, so, it is important to call this
//  		service in your fcl file
//		stick this line in your "services" section of fcl file:
//		FileCatalogMetadata:  @local::art_file_catalog_mc
//	
//              2. When you call FileCatalogMetadata service in your fcl file, and if 
//		you have (art) root Output section in your fcl file, and if you do not  
//		have "dataTier" specified in that section, then this service will throw
//		an exception. To avoid this, either remove the entire root Output section
//		in your fcl file (and remove art stream output from your end_paths) or 
//		include appropriate dataTier information in the section.If you are only
//		running analysis job, best way is to not include any art root Output section.
//
//	        3. This service is exclusively written to work with production (in other
//		words for jobs submitted through grid). Some of the metadata parameters
//		(output TFileName, filesize, Project related details) are captured/updated
//		during and/or after the workflow. 
//	
//
// Created:  21-Feb-2018,  D. Brailsford
//  based on the SBND version by T. Junk which is based on the 
//  based on the MicroBooNE example by S. Gollapinni
//
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include "sbncode/Metadata/TFileMetadataSBN.h"
#include "sbncode/Metadata/FileCatalogMetadataSBN.h"

#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "art_root_io/RootDB/SQLErrMsg.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art/Utilities/OutputFileInfo.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTimeStamp.h"

using namespace std;

//--------------------------------------------------------------------

// Constructor.
util::TFileMetadataSBN::TFileMetadataSBN(fhicl::ParameterSet const& pset, 
					   art::ActivityRegistry& reg):
  fGenerateTFileMetadata{pset.get<bool>("GenerateTFileMetadata")},
  fJSONFileName{pset.get<std::string>("JSONFileName")},
  fFileStats{"", art::ServiceHandle<art::TriggerNamesService const>{}->getProcessName()}
{
  if (fGenerateTFileMetadata) {
    md.fdata_tier   = pset.get<std::string>("dataTier");
    md.ffile_format = pset.get<std::string>("fileFormat");

    reg.sPostBeginJob.watch(this, &TFileMetadataSBN::postBeginJob);
    reg.sPostOpenFile.watch(this, &TFileMetadataSBN::postOpenInputFile);
    reg.sPostCloseFile.watch(this, &TFileMetadataSBN::postCloseInputFile);
    reg.sPostProcessEvent.watch(this, &TFileMetadataSBN::postEvent);
    reg.sPostBeginSubRun.watch(this, &TFileMetadataSBN::postBeginSubRun);
  }
}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via TFileMetadata service.
void util::TFileMetadataSBN::postBeginJob()
{ 
  // get the start time  
  md.fstart_time = time(0); 
  
  // Get art metadata service and extract paramters from there
  art::ServiceHandle<art::FileCatalogMetadata> artmds;
  
  art::FileCatalogMetadata::collection_type artmd;
  artmds->getMetadata(artmd);
  
  for(auto const & d : artmd)
    mdmap[d.first] = d.second;
    
  std::map<std::string,std::string>::iterator it;
  
  // if a certain paramter/key is not found, assign an empty string value to it
  
  if ((it=mdmap.find("application.family"))!=mdmap.end()) std::get<0>(md.fapplication) = it->second;
  else std::get<0>(md.fapplication) = "\" \"";   
   
  if ((it=mdmap.find("art.process_name"))!=mdmap.end()) std::get<1>(md.fapplication) = it->second;
  else std::get<1>(md.fapplication) = "\" \"";  
  
  if ((it=mdmap.find("application.version"))!=mdmap.end()) std::get<2>(md.fapplication) = it->second;
  else  std::get<2>(md.fapplication) = "\" \"";   
  
  if ((it=mdmap.find("group"))!=mdmap.end()) md.fgroup = it->second;
  else md.fgroup = "\" \"";   
    
  if ((it=mdmap.find("file_type"))!=mdmap.end()) md.ffile_type = it->second;
  else  md.ffile_type = "\" \"";  
    
  if ((it=mdmap.find("art.run_type"))!=mdmap.end()) frunType = it->second;
  else frunType = "\" \"";         	     	

}


//--------------------------------------------------------------------        
// PostOpenFile callback.
void util::TFileMetadataSBN::postOpenInputFile(std::string const& fn)
{
  // save parent input files here
  // 08/06 DBrailsford: Only save the parent string if the string is filled.  The string still exists (with 0 characters) for generation stage files.  See redmine issue 20124
  if (fn.length() > 0) md.fParents.insert(fn);
  fFileStats.recordInputFile(fn);
}

//--------------------------------------------------------------------  	
// PostEvent callback.
void util::TFileMetadataSBN::postEvent(art::Event const& evt, art::ScheduleContext)
{
  art::RunNumber_t run = evt.run();
  art::SubRunNumber_t subrun = evt.subRun();
  art::EventNumber_t event = evt.event();
  art::SubRunID srid = evt.id().subRunID();
      
  // save run, subrun and runType information once every subrun    
  if (fSubRunNumbers.count(srid) == 0){
    fSubRunNumbers.insert(srid);
    md.fruns.push_back(make_tuple(run, subrun, frunType));   
  }
  
  // save the first event
  if (md.fevent_count == 0) md.ffirst_event = event;
  md.flast_event = event;
  // event counter
  ++md.fevent_count;
    
}

//--------------------------------------------------------------------  	
// PostSubRun callback.
void util::TFileMetadataSBN::postBeginSubRun(art::SubRun const& sr)
{
  art::RunNumber_t run = sr.run();
  art::SubRunNumber_t subrun = sr.subRun();
  art::SubRunID srid = sr.id();

  // save run, subrun and runType information once every subrun
  if (fSubRunNumbers.count(srid) == 0){
    fSubRunNumbers.insert(srid);
    md.fruns.push_back(make_tuple(run, subrun, frunType));
  }
}

//--------------------------------------------------------------------
// PostCloseFile callback.
void util::TFileMetadataSBN::postCloseInputFile()
{
  // get metadata from the FileCatalogMetadataSBN service, which is filled on its construction
        
  art::ServiceHandle<util::FileCatalogMetadataSBN> paramhandle; 
  md.fFCLName = paramhandle->GetFCLName();
  md.fProjectName = paramhandle->GetProjectName();
  md.fProjectStage = paramhandle->GetProjectStage();
  md.fProjectVersion = paramhandle->GetProjectVersion();
  md.fProjectSoftware = paramhandle->GetProjectSoftware();
  md.fProductionName = paramhandle->GetProductionName();
  md.fProductionType = paramhandle->GetProductionType();

  //update end time
  md.fend_time = time(0);

  // convert start and end times into time format: Year-Month-DayTHours:Minutes:Seconds
  char endbuf[80], startbuf[80];
  struct tm tstruct;
  tstruct = *localtime(&md.fend_time);
  strftime(endbuf,sizeof(endbuf),"%Y-%m-%dT%H:%M:%S",&tstruct);
  tstruct = *localtime(&md.fstart_time);
  strftime(startbuf,sizeof(startbuf),"%Y-%m-%dT%H:%M:%S",&tstruct);
  
  // open a json file and write everything from the struct md complying to the 
  // samweb json format. This json file holds the below information temporarily. 
  // If you submitted a grid job invoking this service, the information from 
  // this file is appended to a final json file and this file will be removed
  
  std::ofstream jsonfile;
  jsonfile.open(fJSONFileName);
  jsonfile<<"{\n  \"application\": {\n    \"family\": "<<std::get<0>(md.fapplication)<<",\n    \"name\": ";
  jsonfile<<std::get<1>(md.fapplication)<<",\n    \"version\": "<<std::get<2>(md.fapplication)<<"\n  },\n  ";
  jsonfile<<"\"data_tier\": \""<<md.fdata_tier<<"\",\n  ";
  jsonfile<<"\"event_count\": "<<md.fevent_count<<",\n  ";
  jsonfile<<"\"file_format\": \""<<md.ffile_format<<"\",\n  ";
  jsonfile<<"\"file_type\": "<<md.ffile_type<<",\n  ";
  jsonfile<<"\"first_event\": "<<md.ffirst_event<<",\n  ";
  jsonfile<<"\"group\": "<<md.fgroup<<",\n  ";
  jsonfile<<"\"last_event\": "<<md.flast_event<<",\n  ";
  //if (md.fdataTier != "generated"){
  unsigned int c=0;
  //08/06 DBrailsford.  Only create the parent json object if there are parent names in the set.
  if (md.fParents.size() > 0){
    jsonfile<<"\"parents\": [\n";
    for(auto parent : md.fParents) {
      std::cout<<"Parent " << c << ": " << parent << std::endl;
      c++;
      size_t n = parent.find_last_of('/');
      size_t f1 = (n == std::string::npos ? 0 : n+1);
      jsonfile<<"    {\n     \"file_name\": \""<<parent.substr(f1)<<"\"\n    }";
      if (md.fParents.size()==1 || c==md.fParents.size()) jsonfile<<"\n";
      else jsonfile<<",\n"; 
    }      
    jsonfile<<"  ],\n  "; 
    //}   
    c=0;
  }
  jsonfile<<"\"runs\": [\n";
  for(auto &t : md.fruns){
    c++;
    jsonfile<<"    [\n     "<<std::get<0>(t)<<",\n     "<<std::get<1>(t)<<",\n     "<<std::get<2>(t)<<"\n    ]";
    if (md.fruns.size()==1 || c==md.fruns.size()) jsonfile<<"\n";
    else jsonfile<<",\n"; 
  }
  jsonfile<<"  ],\n";          

  if (md.fFCLName!="") jsonfile << "\"fcl.name\": \"" << md.fFCLName << "\",\n";
  if (md.fProjectName!="") jsonfile << "\"sbnd_project.name\": \"" << md.fProjectName << "\",\n";
  if (md.fProjectStage!="") jsonfile << "\"sbnd_project.stage\": \"" << md.fProjectStage << "\",\n";
  if (md.fProjectVersion!="") jsonfile << "\"sbnd_project.version\": \"" << md.fProjectVersion << "\",\n";
  if (md.fProjectSoftware!="") jsonfile << "\"sbnd_project.software\": \"" << md.fProjectSoftware << "\",\n";
  if (md.fProductionName!="") jsonfile << "\"production.name\": \"" << md.fProductionName << "\",\n";
  if (md.fProductionType!="") jsonfile << "\"production.type\": \"" << md.fProductionType << "\",\n";

  // put these at the end because we know they'll be there and the last one needs to not have a comma
  jsonfile<<"\"start_time\": \""<<startbuf<<"\",\n";
  jsonfile<<"\"end_time\": \""<<endbuf<<"\"\n";
  
  
  jsonfile<<"}\n";
  jsonfile.close();  

  fFileStats.recordFileClose();
  //TODO figure out how to make the name identical to the TFile
  //std::string new_name = fRenamer.maybeRenameFile("myjson.json",fJSONFileName);
}

DEFINE_ART_SERVICE(util::TFileMetadataSBN)
