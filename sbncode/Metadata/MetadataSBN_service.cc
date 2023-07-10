////////////////////////////////////////////////////////////////////////
// Name:  MetadataSBN_service.cc.
//
// Purpose:  generate SBN-specific sam metadata for root Tfiles (histogram or ntuple files).
//
// FCL parameters: Experiment: Experiment name ("sbnd" or "icarus").
//	     	   JSONFileName: Name of generated .json file(s).
//                 dataTier: Data tier(s).
//		     	     for ntuples, dataTier = root-tuple;
//		             for histos, dataTier = root-histogram
//		             (default value: root-tuple)
//	           fileFormat: This is currently specified by the user,
//			       the fileFormat for Tfiles is "root" (default value: root)
//                 Merge: Merge flag.
//                         1 - Set merge.merge = 1 and merge.merged = 0
//                         0 - Set merge.merge = 0 and merge.merged = 0
//                        -1 - Do not generate merge parameters.
//                 POTModuleLabel  - POTSummary module label (default "generator").
//
//                 Parameters JSONFileName, dataTier, and fileFormat can be single
//                 stringss or sequences of strings.  In case of sequences of length
//                 greater than one, multiple json files will be generated.  Sequences
//                 must be equal length.
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
#include <string>
#include <vector>

#include "sbncode/Metadata/MetadataSBN.h"
#include "sbncode/Metadata/FileCatalogMetadataSBN.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
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
util::MetadataSBN::MetadataSBN(fhicl::ParameterSet const& pset,
			       art::ActivityRegistry& reg):
  fFileStats{"", art::ServiceHandle<art::TriggerNamesService const>{}->getProcessName()}
{
  // Insist on configuring Experiment from the fcl file (ideally) or the
  // environment.
  const char* expt = getenv("EXPERIMENT");
  if(expt) fExperiment = pset.get<std::string>("Experiment", expt); else fExperiment = pset.get<std::string>("Experiment");
  std::transform(fExperiment.begin(), fExperiment.end(), fExperiment.begin(), [](unsigned char c){return std::tolower(c);});

  // Get scalar/vector parameters.
  if(pset.is_key_to_atom("JSONFileName"))
    fJSONFileName.push_back(pset.get<std::string>("JSONFileName"));
  else if(pset.is_key_to_sequence("JSONFileName"))
    fJSONFileName = pset.get<std::vector<std::string> >("JSONFileName");
    
  if(pset.is_key_to_atom("dataTier"))
    fDataTier.push_back(pset.get<std::string>("dataTier"));
  else if(pset.is_key_to_sequence("dataTier"))
    fDataTier = pset.get<std::vector<std::string> >("dataTier");
  if(fDataTier.size() != fJSONFileName.size())
    throw cet::exception("MetadataSBN") << "FCL sequence size mismatch.\n";

  if(pset.is_key_to_atom("fileFormat"))
    fFileFormat.push_back(pset.get<std::string>("fileFormat"));
  else if(pset.is_key_to_sequence("fileFormat"))
    fFileFormat = pset.get<std::vector<std::string> >("fileFormat");
  if(fFileFormat.size() != fJSONFileName.size())
    throw cet::exception("MetadataSBN") << "FCL sequence size mismatch.\n";

  if(pset.has_key("Merge")) {
    if(pset.is_key_to_atom("Merge"))
      fMerge.push_back(pset.get<int>("Merge"));
    else if(pset.is_key_to_sequence("Merge"))
      fMerge = pset.get<std::vector<int> >("Merge");
  }
  fPOTModuleLabel = pset.get<std::string>("POTModuleLabel", "generator");

  reg.sPostBeginJob.watch(this, &MetadataSBN::postBeginJob);
  reg.sPostOpenFile.watch(this, &MetadataSBN::postOpenInputFile);
  reg.sPostCloseFile.watch(this, &MetadataSBN::postCloseInputFile);
  reg.sPostProcessEvent.watch(this, &MetadataSBN::postEvent);
  reg.sPostBeginSubRun.watch(this, &MetadataSBN::postBeginSubRun);
  reg.sPostEndSubRun.watch(this, &MetadataSBN::postEndSubRun);

  // get metadata from the FileCatalogMetadataSBN service, which is filled on its construction
  art::ServiceHandle<util::FileCatalogMetadataSBN> paramhandle;
  md.fFCLName = paramhandle->GetFCLName();
  md.fProjectName = paramhandle->GetProjectName();
  md.fProjectStage = paramhandle->GetProjectStage();
  md.fProjectVersion = paramhandle->GetProjectVersion();
  md.fProjectSoftware = paramhandle->GetProjectSoftware();
  md.fProductionName = paramhandle->GetProductionName();
  md.fProductionType = paramhandle->GetProductionType();
  md.merge = -1;
  md.fTotPOT = 0.;
}

/// Un-quote quoted strings
std::string UnQuoteString(std::string s)
{
  if(s.size() < 2 || s[0] != '\"' || s[s.size()-1] != '\"') return s;
  s.erase(0, 1);
  s.erase(s.size()-1, 1);
  return s;
}

void MaybeCopyFromMap(const std::map<std::string, std::string>& in,
                      const std::string& key,
                      std::string& out)
{
  const auto it = in.find(key);
  if(it == in.end()){
    out = "";
  }
  else{
    out = UnQuoteString(it->second);
  }
}

void MaybeCopyToMap(const std::string& in,
                    const std::string& key,
                    std::map<std::string, std::string>& out)
{
  if(!in.empty()) out[key] = UnQuoteString(in);
}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via Metadata service.
void util::MetadataSBN::postBeginJob()
{
  // get the start time
  md.fstart_time = time(0);

  // Get art metadata service and extract paramters from there
  art::ServiceHandle<art::FileCatalogMetadata> artmds;

  art::FileCatalogMetadata::collection_type artmd;
  artmds->getMetadata(artmd);

  std::map<std::string, std::string> mdmap;
  for(const auto& d: artmd)
    mdmap[d.first] = UnQuoteString(d.second);

  // if a certain paramter/key is not found, assign an empty string value to it
  MaybeCopyFromMap(mdmap, "application.family",  std::get<0>(md.fapplication));
  MaybeCopyFromMap(mdmap, "art.process_name",    std::get<1>(md.fapplication));
  MaybeCopyFromMap(mdmap, "application.version", std::get<2>(md.fapplication));
  MaybeCopyFromMap(mdmap, "group",        md.fgroup);
  MaybeCopyFromMap(mdmap, "file_type",    md.ffile_type);
  MaybeCopyFromMap(mdmap, "art.run_type", frunType);
}


//--------------------------------------------------------------------
// PostOpenFile callback.
void util::MetadataSBN::postOpenInputFile(std::string const& fn)
{
  // save parent input files here
  // 08/06 DBrailsford: Only save the parent string if the string is filled.  The string still exists (with 0 characters) for generation stage files.  See redmine issue 20124
  if (fn.length() > 0) md.fParents.insert(fn);
  fFileStats.recordInputFile(fn);
}

//--------------------------------------------------------------------
// PostEvent callback.
void util::MetadataSBN::postEvent(art::Event const& evt, art::ScheduleContext)
{
  std::lock_guard lock(fMutex);

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
void util::MetadataSBN::postBeginSubRun(art::SubRun const& sr)
{
  std::lock_guard lock(fMutex);

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
// PostEndSubRun callback.
void util::MetadataSBN::postEndSubRun(art::SubRun const& sr)
{
  std::lock_guard lock(fMutex);

  art::Handle< sumdata::POTSummary > potListHandle;
  double fTotPOT = 0;
  if(sr.getByLabel(fPOTModuleLabel,potListHandle)){
    fTotPOT+=potListHandle->totpot;
  }

  md.fTotPOT += fTotPOT;
}


//--------------------------------------------------------------------
std::string Escape(const std::string& s)
{
  // If it's formatted as a dict or list, trust it's already formatted
  if(s.size() >= 2 && ((s[0] == '{' && s.back() == '}') || (s[0] == '[' && s.back() == ']'))) return s;

  // otherwise quote it
  return "\""+s+"\"";
}

//--------------------------------------------------------------------
std::string util::MetadataSBN::GetParentsString() const
{
  if(md.fParents.empty()) return "";

  unsigned int c = 0;

  std::string ret = "[\n";
  for(auto parent: md.fParents) {
    std::cout<<"Parent " << c << ": " << parent << std::endl;
    c++;
    size_t n = parent.find_last_of('/');
    size_t f1 = (n == std::string::npos ? 0 : n+1);
    ret += "    {\n     \"file_name\": \"" + parent.substr(f1) + "\"\n    }";
    if(md.fParents.size() == 1 || c == md.fParents.size()) ret += "\n";
    else ret += ",\n";
  }

  ret += "  ]";
  return ret;
}

//--------------------------------------------------------------------
std::string util::MetadataSBN::GetRunsString() const
{
  unsigned int c = 0;

  std::string ret = "[\n";
  for(auto&t :md.fruns){
    c++;
    ret += "    [\n     " + std::to_string(std::get<0>(t)) + ",\n     " + std::to_string(std::get<1>(t)) + ",\n     \"" + std::get<2>(t) + "\"\n    ]";
    if(md.fruns.size() == 1 || c == md.fruns.size()) ret += "\n";
    else ret += ",\n";
  }
  ret += "  ]";
  return ret;
}

//--------------------------------------------------------------------
void util::MetadataSBN::GetMetadataMaps(std::map<std::string, std::string>& strs,
                                        std::map<std::string, int>& ints,
					std::map<std::string, double>& doubles,
                                        std::map<std::string, std::string>& objs)
{
  strs.clear(); ints.clear(); doubles.clear(); objs.clear();

  objs["application"] = "{\"family\": \""+std::get<0>(md.fapplication)+"\", \"name\": \""+std::get<1>(md.fapplication)+"\", \"version\": \""+std::get<2>(md.fapplication)+"\"}";

  if(!md.fParents.empty()) objs["parents"] = GetParentsString();
  if(!md.fruns.empty()) objs["runs"] = GetRunsString();

  // convert start and end times into time format: Year-Month-DayTHours:Minutes:Seconds
  char endbuf[80], startbuf[80];
  struct tm tstruct;
  tstruct = *localtime(&md.fend_time);
  strftime(endbuf,sizeof(endbuf),"%Y-%m-%dT%H:%M:%S",&tstruct);
  tstruct = *localtime(&md.fstart_time);
  strftime(startbuf,sizeof(startbuf),"%Y-%m-%dT%H:%M:%S",&tstruct);

  strs["start_time"] = startbuf;
  strs["end_time"] = endbuf;

  strs["data_tier"] = md.fdata_tier;
  ints["event_count"] = md.fevent_count;
  strs["file_format"] = md.ffile_format;
  ints["first_event"] = md.ffirst_event;
  ints["last_event"] = md.flast_event;

  const std::string proj = fExperiment+"_project";
  MaybeCopyToMap(md.fFCLName, "fcl.name", strs);
  MaybeCopyToMap(md.fProjectName, proj+".name", strs);
  MaybeCopyToMap(md.fProjectStage, proj+".stage", strs);
  MaybeCopyToMap(md.fProjectVersion, proj+".version", strs);
  MaybeCopyToMap(md.fProjectSoftware, proj+".software", strs);
  MaybeCopyToMap(md.fProductionName, "production.name", strs);
  MaybeCopyToMap(md.fProductionType, "production.type", strs);

  MaybeCopyToMap(md.fgroup, "group", strs);
  MaybeCopyToMap(md.ffile_type, "file_type", strs);

  if(md.merge >= 0) {
    ints["merge.merge"] = (md.merge==0 ? 0 : 1);
    ints["merge.merged"] = 0;
  }
  doubles["mc.pot"] = md.fTotPOT;  
}

//--------------------------------------------------------------------
// PostCloseFile callback.
void util::MetadataSBN::postCloseInputFile()
{
  //update end time
  md.fend_time = time(0);

  // Loop over files.

  for(unsigned int i=0; i<fJSONFileName.size(); ++i) {

    // Update per-file metadata.

    md.fdata_tier = fDataTier[i];
    md.ffile_format = fFileFormat[i];
    if(fMerge.size() > i)
      md.merge = fMerge[i];
    else
      md.merge = -1;


    std::map<std::string, std::string> strs;
    std::map<std::string, int> ints;
    std::map<std::string, double> doubles;
    std::map<std::string, std::string> objs;
    GetMetadataMaps(strs, ints, doubles, objs);

    // open a json file and write everything from the struct md complying to the
    // samweb json format. This json file holds the below information temporarily.
    // If you submitted a grid job invoking this service, the information from
    // this file is appended to a final json file and this file will be removed

    if(!fJSONFileName[i].empty()){
      std::ofstream jsonfile;
      jsonfile.open(fJSONFileName[i]);
      jsonfile << "{\n";

      bool once = true;
      for(auto& it: objs){
	if(!once) jsonfile << ",\n";
	once = false;
	jsonfile << "  \"" << it.first << "\": " << it.second;
      }
      for(auto& it: strs){
	// Have to escape string outputs
	jsonfile << ",\n  \"" << it.first << "\": \"" << it.second << "\"";
      }
      for(auto& it: ints){
	jsonfile << ",\n  \"" << it.first << "\": " << it.second;
      }
      for(auto& it: doubles){
	jsonfile << ",\n  \"" << it.first << "\": " << it.second;
      }

      jsonfile<<"\n}\n";
      jsonfile.close();
    }
  }

  fFileStats.recordFileClose();
  //TODO figure out how to make the name identical to the TFile
  //std::string new_name = fRenamer.maybeRenameFile("myjson.json",fJSONFileName);
}

DEFINE_ART_SERVICE(util::MetadataSBN)
