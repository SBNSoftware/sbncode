#include "sbnanalysis/core/Event.hh"

#include "sbncode/FlatMaker/FlatRecoEvent.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv)
{
  // Have to do it here since we didn't figure out how to statically link it
  // yet
  gSystem->Load("libsbnanalysis_Event.so");

  if(argc != 3){
    std::cout << "Usage: convert_to_flat input.events.root output.flat.root"
              << std::endl;
    return 1;
  }

  std::string inname = argv[1];
  std::string outname = argv[2];

  // Can manage by passing '' ''
  if(inname.empty()) inname = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/output_SBNOsc_NumuSelection_Modern_Uboone.root";
  if(outname.empty()) outname = "output_SBNOsc_NumuSelection_Modern_Uboone.flat.root";

  TFile* fin = TFile::Open(inname.c_str());

  TTree* tr = (TTree*)fin->Get("sbnreco");
  if(!tr){
    std::cout << "Couldn't find tree 'sbnreco' in " << inname << std::endl;
    return 1;
  }

  event::RecoEvent* event = 0;
  tr->SetBranchAddress("reco_events", &event);

  TFile fout(outname.c_str(), "RECREATE");

  // First, copy the POT information over
  ((TTree*)fin->Get("sbnsubrun"))->CloneTree()->Write("sbnsubrun");

  TTree* trout = new TTree("sbnana", "sbnana");
  // Have trouble with memory usage (because several trees are open at
  // once?). Set the maximum buffer usage (per tree) to 3MB (10x less than
  // default)
  trout->SetAutoFlush(-3*1000*1000);

  flat::FlatRecoEvent* rec = new flat::FlatRecoEvent("sbnana.", trout, 0);//policy);

  for(int i = 0; i < tr->GetEntries(); ++i){
    // TODO use CAFAna progress bar
    if(i%100 == 0) std::cout << i/100 << " / " << tr->GetEntries()/100 << std::endl;

    tr->GetEntry(i);

    rec->Fill(*event);
    trout->Fill();
  }

  trout->Write();
  delete rec;

  return 0;
}
