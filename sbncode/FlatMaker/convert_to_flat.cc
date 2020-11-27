#include "sbncode/StandardRecord/StandardRecord.h"

#include "sbncode/FlatMaker/FlatRecord.h"

#include "sbncode/CAFAna/Core/Progress.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv)
{
  if(argc != 3){
    std::cout << "Usage: convert_to_flat input.events.root output.flat.root"
              << std::endl;
    return 1;
  }

  const std::string inname = argv[1];
  const std::string outname = argv[2];

  //Find out if "proposal" appears in input string
  bool is_proposal_flag = false;
  if (inname.find("Proposal") || inname.find("proposal")) {
    is_proposal_flag = true;
    std::cout << "Setting proposal flag to true: SBND baseline will be adjusted" << std::endl;
  }

  TFile* fin = TFile::Open(inname.c_str());

  TTree* tr = (TTree*)fin->Get("recTree");
  if(!tr){
    std::cout << "Couldn't find tree 'recTree' in " << inname << std::endl;
    return 1;
  }

  caf::StandardRecord* event = 0;
  tr->SetBranchAddress("rec", &event);

  // LZ4 is the fastest format to decompress. I get 3x faster loading with
  // this compared to the default, and the files are only slightly larger.
  TFile fout(outname.c_str(), "RECREATE", "",
             ROOT::CompressionSettings(ROOT::kLZ4, 1));

  TTree* trout = new TTree("recTree", "recTree");
  // On NOvA, had trouble with memory usage (because several trees are open at
  // once?). Setting the maximum buffer usage (per tree) to 3MB (10x less than
  // default) fixed it. But it doesn't seem necessary for now on SBN.
  //  trout->SetAutoFlush(-3*1000*1000);

  flat::Flat<caf::StandardRecord> rec(trout, "rec", "", 0);//policy);

  ana::Progress prog("Converting '"+inname+"' to '"+outname+"'");
  for(int i = 0; i < tr->GetEntries(); ++i){
    prog.SetProgress(double(i)/tr->GetEntries());

    tr->GetEntry(i);
    if(is_proposal_flag && !event->mc.nu.empty()){
      const float bl = event->mc.nu[0].baseline;
      if(bl < 150) event->mc.nu[0].baseline -= 10;
    }

    rec.Clear();
    rec.Fill(*event);
    trout->Fill();
  }
  prog.Done();

  trout->Write();

  TH1* hPOT = (TH1*)fin->Get("TotalPOT");
  TH1* hEvts = (TH1*)fin->Get("TotalEvents");
  fout.cd();
  hPOT->Write("TotalPOT");
  hEvts->Write("TotalEvents");

  return 0;
}
