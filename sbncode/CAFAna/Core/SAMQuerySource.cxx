#include "CAFAna/Core/SAMQuerySource.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/Utilities.h"

#include "ifdh.h"

#include <cassert>
#include <iostream>
#include <map>
#include <sys/stat.h>
#include <unistd.h>

#include "TString.h"

namespace ana
{
  //----------------------------------------------------------------------
  SAMQuerySource::SAMQuerySource(const std::string& query,
                                 int stride, int offset)
    // Stride and offset already taken account of in the query
    : FileListSource(LocationsForSAMQuery(query, stride, offset), 1, 0)
  {
  }

  //----------------------------------------------------------------------
  SAMQuerySource::~SAMQuerySource()
  {
  }

  //----------------------------------------------------------------------
  bool SAMQuerySource::RunningOnGrid() const
  {
    return getenv("_CONDOR_SCRATCH_DIR") != 0;
  }

  //----------------------------------------------------------------------
  std::string SAMQuerySource::EnsureDataset(const std::string& query) const
  {
    const char* user = getenv("GRID_USER");
    assert(user);

    TString dset = TString::Format("%s_cafana_%s", user, query.c_str());
    // Sanitize various special characters that can appear in queries
    dset.ReplaceAll(" ", "_");
    dset.ReplaceAll("(", "_OPEN_");
    dset.ReplaceAll(")", "_CLOSE_");
    dset.ReplaceAll(":", "_COLON_");
    dset.ReplaceAll("'", "_SQUOTE_");
    dset.ReplaceAll("\"", "_DQUOTE_");

    std::cout << "Creating dataset " << dset << " for query " << query << std::endl;

    // I would be much much happier to do this in proper code, but I'm not sure
    // how, there's no samweb C++ API?
    system(TString::Format("samweb list-definitions --defname %s | grep %s || samweb create-definition %s %s",
                           dset.Data(), dset.Data(), dset.Data(), query.c_str()).Data());

    return dset.Data();
  }

  //----------------------------------------------------------------------
  std::string SAMQuerySource::EnsureSnapshot(const std::string& def) const
  {
    const char* user = getenv("GRID_USER");
    assert(user);
    const char* cluster = getenv("CLUSTER");
    assert(cluster);
    const char* process = getenv("PROCESS");
    assert(process);

    // Jobs in the same cluster should share the same snapshot of the dataset
    // so as not to hammer SAM with multiple requests for the same file list,
    // but so that the dataset snapshot is updated with every new submission.
    const std::string snap = TString::Format("%s_cafana_snap_%s_%s",
                                             user, def.c_str(), cluster).Data();

    // I'd love to do all this with a proper API, but samweb doesn't seem to
    // have a C++ one? So we get this stew of system() calls...

    // Use this name as an indication that someone is working on creating the
    // snapshot and every one else should stand by.
    const std::string snaplock = TString::Format("%s_cafana_snap_lock_%s_%s",
                                                 user, def.c_str(), cluster).Data();

    // Try to create the lock. Success means we have to create the snapshot,
    // failure means someone else is working on it. The content of the
    // definition (the nova.special) doesn't matter, except it has to be unique
    // between the jobs, because trying to create an exact duplicate of an
    // existing definition counts as success.
    std::cout << "Checking lock " << snaplock << std::endl;
    if(system(TString::Format("samweb create-definition %s nova.special %s",
                              snaplock.c_str(), process).Data()) == 0){
      // No one took the lock, it's up to us. Make the actual snapshot
      std::cout << "Snapshotting " << def << " as " << snap << std::endl;
      system(TString::Format("samweb take-snapshot %s | xargs samweb create-definition %s snapshot_id",
                             def.c_str(), snap.c_str()).Data());
    }
    else{
      // Lock already exists, just wait for the real snapshot to be created
      double period = 1;
      while(system(TString::Format("samweb list-definitions --defname %s | grep %s",
                                   snap.c_str(), snap.c_str()).Data()) != 0){
        sleep(int(period));
        period *= 1.5;
        if(period > 60*30){
          std::cout << "We've been waiting a real long time for " << snap << " to be created. I don't think it's happening." << std::endl;
          abort();
        }
      }
    }
    std::cout << "Will use " << snap << std::endl;

    return snap;
  }

  //----------------------------------------------------------------------
  std::vector<std::string> SAMQuerySource::
  LocationsForSAMQuery(const std::string& str, int stride, int offset)
  {
    TString query = str;

    // This is an actual query
    if(query.Contains(' ') && RunningOnGrid()){
      // On the grid we want to convert that to a dataset we can snapshot below
      query = EnsureDataset(query.Data());
    }

    // This is a dataset name
    if(!query.Contains(' ')){
      if(getenv("CAFANA_USE_SNAPSHOTS")){
	query = "dataset_def_name_newest_snapshot "+query;
      }
      else{
	// Take one snapshot between all the jobs and share that
	if(RunningOnGrid()) query = EnsureSnapshot(query.Data());

	query = "defname: "+query;
      }
    }

    if(stride > 1){
      query += TString::Format(" with stride %d", stride).Data();
      if(offset > 0){
        query += TString::Format(" offset %d", offset).Data();
      }
    }


    std::cout << "Looking up files matching '" << query << "' using SAM...\n";

    std::vector<std::string> files;

    ifdh i;
    i.set_debug("0"); // shut up
    try{
      files = i.translateConstraints(query.Data());
    }
    catch(ifdh_util_ns::WebAPIException& e){
      // I like my error message below better, since this could well be a
      // mistyped filename.
    }

    if(files.empty()){
      std::cerr << "\nCan't find any files matching '" << str
		<< "'. Aborting..." << std::endl;
      abort();
    }

    // IFDH likes to give back an empty string as the last response
    // https://cdcvs.fnal.gov/redmine/issues/6718
    if(!files.empty() && files.back().empty()){
      files.pop_back();
    }

    return LocateSAMFiles(files);
  }

  //----------------------------------------------------------------------
  std::vector<std::string> SAMQuerySource::
  LocateSAMFiles(const std::vector<std::string>& fnames)
  {
    std::vector<std::string> ret;

    // We're going to fill this map of locations for all the files
    std::map<std::string, std::vector<std::string>> locmap;

    Progress prog(TString::Format("Looking up locations of %ld files using SAM", fnames.size()).Data());

    ifdh i;
    i.set_debug("0"); // shut up

    // locateFiles() saves the roundtrip time of talking to the server about
    // every file individually, but it seems to bog down for large
    // queries. Split the query into chunks. Experimentally this is about the
    // sweet spot.
    const unsigned int kStep = 50;
    for(unsigned int fIdx = 0; fIdx < fnames.size(); fIdx += kStep){
      prog.SetProgress(double(fIdx)/fnames.size());

      // The files we're looking up right now. Careful not to run off the end
      // of the vector.
      const std::vector<std::string> fslice(fnames.begin()+fIdx, fIdx+kStep < fnames.size() ? fnames.begin()+fIdx+kStep : fnames.end());

      const auto locslice = i.locateFiles(fslice);

      locmap.insert(locslice.begin(), locslice.end());
    }

    prog.Done();


    // Now go through the map and pick our favourite location for each file,
    // and do some cleanup.
    for(auto it: locmap){
      const std::string& f = it.first;
      const std::vector<std::string>& locs = it.second;

      int best = 0;

      std::string resolved;
      for(TString loc: locs){
	// Never try to access bluearc locations from the grid
	if(!RunningOnGrid() && loc.BeginsWith("novadata:") && best < 3){
	  loc.Remove(0, 9);

          // Rewrite FNAL bluearc paths so users with matching directory
          // structures offsite can access their local copies.
          if(std::getenv("NOVA_ANA" )) loc.ReplaceAll("/nova/ana",  std::getenv("NOVA_ANA"));
          if(std::getenv("NOVA_APP" )) loc.ReplaceAll("/nova/app",  std::getenv("NOVA_APP"));
          if(std::getenv("NOVA_DATA")) loc.ReplaceAll("/nova/data", std::getenv("NOVA_DATA"));
          if(std::getenv("NOVA_PROD")) loc.ReplaceAll("/nova/prod", std::getenv("NOVA_PROD"));

          // Check if the file exists at that location. If not, maybe pnfs has
          // it.
          struct stat junk;
          if(stat((resolved+'/'+f).c_str(), &junk) == 0){
            best = 3; // Prefer bluearc locations
            resolved = loc;
          }
	}

        if(loc.BeginsWith("dcache:") && best < 2){
          // Otherwise, used xrootd. Prefer "dache:" to "enstore:" because
          // "dcache:" probably means /pnfs/nova/persistent/ so no chance of a
          // huge wait for tape.
          best = 2;

          // FileListSource does the actual conversion to xrootd
          loc.ReplaceAll("dcache:/pnfs/", "/pnfs/");
          // Strip the bit in brackets from the end
          if(loc.First('(') >= 0) loc.Resize(loc.First('('));
          resolved = loc;
        }

	if(loc.BeginsWith("enstore:") && best < 1){
	  best = 1;

	  loc.ReplaceAll("enstore:/pnfs/", "/pnfs/");
	  if(loc.First('(') >= 0) loc.Resize(loc.First('('));
	  resolved = loc;
	}

      } // end for loc

      if(best == 0 || resolved.empty()){
	std::cerr << "\nCouldn't find a usable location for " << f
		  << "'. Aborting..." << std::endl;
	abort();
      }

      ret.push_back((resolved+'/'+f));
    } // end for fIdx

    return ret;
  }
}
