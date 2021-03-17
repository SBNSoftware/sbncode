////////////////////////////////////////////////////////////////////////
// Name: MetadataSBN.h
//
// A struct datatype to hold the metadata information as it is extracted
// from various places.
//
// Created: 21-Feb-2017,  D. Brailsford
//   Based on the DUNE version T. Junk which is based on the MicroBooNE
//   version by S. Gollapinni
//
////////////////////////////////////////////////////////////////////////
#ifndef TFILEMETADATASBN_H
#define TFILEMETADATASBN_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/fwd.h"
#include "art/Framework/IO/PostCloseFileRenamer.h"
#include "art/Framework/IO/FileStatsCollector.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "canvas/Persistency/Provenance/IDNumber.h"
#include "fhiclcpp/ParameterSet.h"

#include <time.h>
#include <fstream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace util{

  class MetadataSBN
  {
  public:
    MetadataSBN(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    struct metadata {
      std::tuple<std::string, std::string, std::string> fapplication;
      //no crc information yet
      //std::vector<std::string> fcrcs;
      std::string fdata_tier;
      time_t fend_time;
      unsigned int fevent_count=0;
      std::string ffile_format;
      std::string ffile_type;
      art::EventNumber_t ffirst_event=0;
      std::string fgroup;
      art::EventNumber_t flast_event=0;
      std::set<std::string> fParents;
      std::vector<std::tuple<art::RunNumber_t,art::SubRunNumber_t,std::string>> fruns;
      time_t fstart_time=0;
      std::string fFCLName;
      std::string fProjectName;
      std::string fProjectStage;
      std::string fProjectVersion;
      std::string fProjectSoftware;
      std::string fProductionName; //Production parameter, do not use if not running a production
      std::string fProductionType; //Production parameter, do not use if not running a production
    };

    metadata md;
    std::set<art::SubRunID> fSubRunNumbers;

    void GetMetadataMaps(std::map<std::string, std::string>& strs,
                         std::map<std::string, int>& ints,
                         std::map<std::string, std::string>& objs);

  private:

    // Callbacks.
    void postBeginJob();
    void postOpenInputFile(std::string const& fn);
    void postEvent(art::Event const& ev, art::ScheduleContext);
    void postBeginSubRun(art::SubRun const& subrun);
    void postCloseInputFile();

    std::string GetParentsString() const;
    std::string GetRunsString() const;

    std::map<std::string,std::string> mdmapStr;
    std::map<std::string, int> mdmapInt;
    std::map<std::string, std::string> mdmapObj;

    // Fcl parameters.
    std::string fExperiment;
    std::string frunType;
    std::string fJSONFileName;
    art::FileStatsCollector fFileStats;
    art::PostCloseFileRenamer fRenamer{fFileStats};
  }; // class MetadataSBN

} //namespace utils

DECLARE_ART_SERVICE(util::MetadataSBN, LEGACY)

#endif
