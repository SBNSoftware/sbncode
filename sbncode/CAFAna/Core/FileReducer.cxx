#include "CAFAna/Core/FileReducer.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/Utilities.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include "StandardRecord/StandardRecord.h"

#include <cassert>
#include <iostream>

#include <fenv.h>

#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"


namespace ana
{
  //----------------------------------------------------------------------
  void ClearTrueParticles(caf::StandardRecord* sr)
  {
    sr->true_particles.clear();
    sr->ntrue_particles = 0;
  }

  //----------------------------------------------------------------------
  FileReducer::FileReducer(const std::string& wildcard,
                           const std::string& outfile)
    : SpectrumLoaderBase(wildcard),
      fOutfile(outfile),
      fSpillCut(nullptr), fSliceCut(nullptr),
      fCopyMetadata(false)
  {
  }

  //----------------------------------------------------------------------
  FileReducer::FileReducer(const std::vector<std::string>& fnames,
                           const std::string& outfile)
    : SpectrumLoaderBase(fnames),
      fOutfile(outfile),
      fSpillCut(nullptr),
      fSliceCut(nullptr),
      fCopyMetadata(false)
  {
  }

  //----------------------------------------------------------------------
  FileReducer::~FileReducer()
  {
    delete fSpillCut;
    delete fSliceCut;
  }

  //----------------------------------------------------------------------
  void FileReducer::AddSpillCut(const SpillCut& cut)
  {
    if(fSpillCut){
      *fSpillCut = *fSpillCut && cut;
    }
    else{
      fSpillCut = new SpillCut(cut);
    }
  }

  //----------------------------------------------------------------------
  void FileReducer::AddSliceCut(const SliceCut& cut)
  {
    if(fSliceCut){
      *fSliceCut = *fSliceCut && cut;
    }
    else{
      fSliceCut = new SliceCut(cut);
    }
  }

  //----------------------------------------------------------------------
  void FileReducer::SetEventList(const std::string& fname)
  {
    FILE* f = fopen(fname.c_str(), "r");
    assert(f);

    while(!feof(f)){
      int run, subrun, event;
      fscanf(f, "%d %d %d", &run, &subrun, &event);
      fEventList.emplace(run, subrun, event);
    }

    fclose(f);
  }

  //----------------------------------------------------------------------
  void FileReducer::Go()
  {
    //    FloatingExceptionOnNaN fpnan;
    
    // Don't want overflow to happen. Set to 1 petabyte: effectively infinite.
    TTree::SetMaxTreeSize(1e15);

    if(fGone){
      std::cerr << "Error: can only call Go() once on a FileReducer" << std::endl;
      abort();
    }
    fGone = true;

    const int Nfiles = NFiles();

    Progress* prog = 0;
    
    TFile fout(fOutfile.c_str(), "RECREATE");
    TTree* trOut = new TTree("recTree", "recTree");
    {
      //      FloatingExceptionOnNaN fpnan(false);
      caf::StandardRecord dummy;
      trOut->Branch("rec", &dummy);
    }

    TH1* hPOTOut = new TH1F("TotalPOT", "", 1, 0, 1);
    TH1* hEventsOut = new TH1F("TotalEvents", "", 1, 0, 1);

    std::vector<std::string> fnames;

    //    std::map<std::string, std::string> meta;
    //    std::set<std::string> meta_mask;

    caf::StandardRecord* oldsr = 0;

    long nRecSeen = 0;
    long nRecPassed = 0;

    int fileIdx = -1;
    while(TFile* f = GetNextFile()){
      ++fileIdx;

      assert(!f->IsZombie());

      TH1* hPOT = (TH1*)f->Get("TotalPOT");

      assert(hPOT);
      hPOTOut->Add(hPOT);

      fnames.push_back(f->GetName());

      /*
      if(fCopyMetadata){
        TDirectory* meta_dir = (TDirectory*)f->Get("metadata");
        assert(meta_dir);
        CombineMetadata(meta, GetCAFMetadata(meta_dir), meta_mask);
      }
      */

      TH1* hEvents = (TH1*)f->Get("TotalEvents");
      assert(hEvents);
      hEventsOut->Add(hEvents);

      TTree* recTree = (TTree*)f->Get("recTree");
      assert(recTree);

      // Use this one for assessing cuts etc
      caf::Proxy<caf::StandardRecord> srProxy(0, recTree, "rec", 0, 0);

      // And if they pass load into this one for writing out
      caf::StandardRecord* sr = 0;
      recTree->SetBranchAddress("rec", &sr);

      const int Nentries = recTree->GetEntries();
      for(int n = 0; n < Nentries; ++n){
        ++nRecSeen;
        recTree->LoadTree(n);

        // Apply EventList cut if it's been enabled
        if(!fEventList.empty() &&
           !fEventList.count(std::make_tuple(srProxy.hdr.run,
                                             srProxy.hdr.subrun,
                                             srProxy.hdr.evt))) continue;

        /// see if we want to omit the event
        if(!fSpillCut || (*fSpillCut)(&srProxy)){
          recTree->GetEntry(n);

          if(sr != oldsr){
            trOut->SetBranchAddress("rec", &sr);
            oldsr = sr;
          }

          std::vector<int> tocut;
          for(unsigned int i = 0; i < srProxy.slc.size(); ++i){
            if(fSliceCut && !(*fSliceCut)(&srProxy.slc[i])) tocut.push_back(i);
          }

          // Remove slices in reverse order so that the indices remain valid
          for(auto it = tocut.rbegin(); it != tocut.rend(); ++it){
            sr->slc.erase(sr->slc.begin() + *it);
          }

          // Apply any additional reduction steps
          for(const auto& f: fReductionFuncs) f(sr);
          // This is kind of problematic since the proxy and actual record
          // could be out of sync. Let's just disable this option for now.
          //          for(const auto & f: fReductionFuncsWithProxy) f(sr, &srProxy);

          ++nRecPassed;
          trOut->Fill();
        }

        if(Nfiles >= 0 && !prog) prog = new Progress(TString::Format("Filling from %d files matching '%s'", Nfiles, fWildcard.c_str()).Data());

        if(n%100 == 0 && Nfiles == 1 && prog)
          prog->SetProgress(double(n)/Nentries);
      } // end for n

      if(prog) prog->SetProgress((fileIdx+1.)/Nfiles);
    } // end while GetNextFile

    fout.cd();
    trOut->Write();
    hPOTOut->Write();
    hEventsOut->Write();

    //    UpdateMetadata(meta, meta_mask, fnames);
    //    WriteCAFMetadata(fout.mkdir("metadata"), meta);

    fout.Close();

    if(prog){
      prog->Done();
      delete prog;
    }

    std::cout << "Passed " << nRecPassed << " / " << nRecSeen << " records";
    std::cout << std::endl;
  }

  //----------------------------------------------------------------------
  /*
  void FileReducer::UpdateMetadata(std::map<std::string, std::string>& meta,
                                   const std::set<std::string>& mask,
                                   const std::vector<std::string>& fnames) const
  {
    for(const std::string& m: mask){
      std::cerr << "Warning: metadata parameter '" << m << "' differs between input files and has been dropped from the output." << std::endl;
      meta.erase(m);
    }

    // change caf -> decaf in the metadata field,
    // if we actually reduced anything
    // and if parents have data_tier
    if ( (!fReductionFuncs.empty() || !fReductionFuncsWithProxy.empty())
        && meta.find("data_tier") != meta.end())
    {
      std::string decaf_tier = meta["data_tier"];
      assert(decaf_tier.size() >= 3);
      assert(decaf_tier.substr(decaf_tier.size()-3,3) == "caf");
      // don't make 'decaf' into 'dedecaf', however
      if (decaf_tier.size() < 5 || decaf_tier.substr(decaf_tier.size()-5,5) != "decaf")
        decaf_tier.replace(decaf_tier.size()-3,3,"decaf");
      meta["data_tier"] = decaf_tier;
    }

    const char* rel = getenv("SRT_BASE_RELEASE");
    if(rel) meta["decaf.base_release"] = rel;

    std::string parents = "[";
    for(const std::string& f: fnames){
      parents += "{\"file_name\":\""+std::string(basename((char *)f.c_str()))+"\"},";
    }
    if(parents[parents.size()-1] == ',') parents.resize(parents.size()-1);
    parents += "]";

    meta["parents"] = parents;

    // if there's one more than one parent, this is a concat.
    // then we need "sumcaf", "sumdecaf", etc. as appropriate
    if (fnames.size() > 1 && meta.find("data_tier") != meta.end())
    {
      auto & tier = meta["data_tier"];
      // if it ends with 'caf' and doesn't start with 'sum',
      // it should now start with 'sum'
      if (tier.substr(tier.size()-3, 3) == "caf" && tier.substr(0, 3) != "sum")
        tier = "sum" + tier;
    }

    // Allow user to override any metadata
    for(auto it: fMetaMap) meta[it.first] = it.second;
  }
  */

} // namespace
