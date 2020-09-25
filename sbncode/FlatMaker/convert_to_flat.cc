#include "sbncode/StandardRecord/StandardRecord.h"

#include "sbncode/FlatMaker/FlatRecord.h"

//#include "sbncode/CAFAna/Core/Progress.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv)
{
  // Have to do it here since we didn't figure out how to statically link it
  // yet
  //  gSystem->Load("libsbnanalysis_Event.so");

  if(argc != 3){
    std::cout << "Usage: convert_to_flat input.events.root output.flat.root"
              << std::endl;
    return 1;
  }

  const std::string inname = argv[1];
  const std::string outname = argv[2];

  //Find out if "proposal" appears in input string
  // int is_proposal_flag = 0;
  // if (inname.find("Proposal") || inname.find("proposal")) {
  //   is_proposal_flag = 1;
  //   std::cout << "Setting proposal flag to 1: SBND baseline will be adjusted" << std::endl;
  // }

  TFile* fin = TFile::Open(inname.c_str());

  TTree* tr = (TTree*)fin->Get("sbnreco");
  if(!tr){
    std::cout << "Couldn't find tree 'sbnreco' in " << inname << std::endl;
    return 1;
  }

  caf::StandardRecord* event = 0;
  tr->SetBranchAddress("reco_events", &event);

  // LZ4 is the fastest format to decompress. I get 3x faster loading with
  // thiscompared to the default, and the files are only slightly larger.
  TFile fout(outname.c_str(), "RECREATE", "",
             ROOT::CompressionSettings(ROOT::kLZ4, 1));

  // First, copy the POT information over
  ((TTree*)fin->Get("sbnsubrun"))->CloneTree()->Write("sbnsubrun");

  TTree* trout = new TTree("sbnana", "sbnana");
  // On NOvA, had trouble with memory usage (because several trees are open at
  // once?). Setting the maximum buffer usage (per tree) to 3MB (10x less than
  // default) fixed it. But it doesn't seem necessary for now on SBN.
  //  trout->SetAutoFlush(-3*1000*1000);

  flat::FlatRecord* rec = new flat::FlatRecord("sbnana.", trout, 0);//policy);

  //  ana::Progress prog("Converting '"+inname+"' to '"+outname+"'");
  for(int i = 0; i < tr->GetEntries(); ++i){
    //    prog.SetProgress(double(i)/tr->GetEntries());

    tr->GetEntry(i);
    //    float bl = event->truth[0].neutrino.baseline;
    //    if (bl < 150. && is_proposal_flag) event->truth[0].neutrino.baseline = bl - 10;

    rec->Fill(*event);
    trout->Fill();

    // This causes us to have much larger baskets, which seems like it should
    // be more efficient for semi-random access of systs, but doesn't seem to
    // help much in practice. It appears the limit is capped at 256e6.
    if(i == 1000) rec->OptimizeBaskets(1000*1000*1000, 1.1, ""/*"d"*/);
  }
  //  prog.Done();

  trout->Write();
  delete rec;

  return 0;
}
